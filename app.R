library("shiny")
library("SemiCompRisks")
library("scales")
library("plyr")
library("dplyr")
library("ggplot2")
library("gridExtra")
library("grid")
library("ggpubr")
library("markdown")


# Magic numbers to slowly delete ------------------------------------------

maxt <- 240
rbarnum <- maxt + 30
maxtplot <- rbarnum + 5
cens <- 10^8
nxts <- 100

# Functions ---------------------------------------------------------------

# Simulate data based on random intercept model (strategy 1)
simulate_RI <- function(n, 
                        theta.true, 
                        alpha1.true, alpha2.true, alpha3.true,
                        beta1.true, beta2.true, beta3.true,
                        kappa1.true, kappa2.true, kappa3.true,
                        cens,
                        frailty = TRUE) {
  
  # Cens -> vector
  if (length(cens) != 2) {
    cens <- c(cens, cens + 1)
  }
  
  # Frailty parameters
  if (frailty) {
    log_gamma_i <- log(rgamma(n, shape = 1/theta.true, rate = 1/theta.true))  
  } else {
    log_gamma_i <- 0
  }
  
  
  # Simulate potential confounder
  x_c <- rnorm(n, mean = 0, sd = 1)
  a <- rep(0:1, each = n)
  x1 <- x2 <- x3 <- cbind(c(x_c, x_c), a, c(log_gamma_i, log_gamma_i))
  
  # Simulate, adding one to each set of parameters
  sim_dat <- simID(cluster = NULL, x1, x2, x3, 
                   alpha1.true = alpha1.true,
                   alpha2.true = alpha2.true,
                   alpha3.true = alpha3.true,
                   beta1.true = c(0, beta1.true, 1),
                   beta2.true = c(0, beta2.true, 1),
                   beta3.true = c(0, beta3.true, 1),
                   kappa1.true = kappa1.true,
                   kappa2.true = kappa2.true,
                   kappa3.true = kappa3.true,
                   theta.true, SigmaV.true = NULL, cens)
  
  # Rename and add treatment column
  colnames(sim_dat) <- c("R", "delta_R", "T", "delta_T")
  sim_dat$ID <- c(1:n, 1:n)
  sim_dat$z <- factor(rep(0:1, each = n))
  
  # Make factor version of Z for nicer plotting
  sim_dat$World <- factor(sim_dat$z, levels = c("1", "0"), 
                          labels = c("Treated (z=1)", "Control (z=0)"))
  
  # Add version of R for plotting (where Rbar has a value)
  sim_dat$Rtoplot <- sim_dat$R
  sim_dat$Rtoplot[sim_dat$delta_R == 0] <- rbarnum * 2
  
  # Return created data
  return(sim_dat)
}

# Get principal state based on values of t, T0, and T1
get_state <- function(T0, T1, at_time) {
  if ((T0 > at_time) & (T1 > at_time)) {
    state_t <- "AA"
  } else if ((T0 > at_time) & (T1 < at_time)) {
    state_t <- "TK"
  } else if ((T0 < at_time) & (T1 > at_time)) {
    state_t <- "CK"
  } else if ((T0 < at_time) & (T1 < at_time)) {
    state_t <- "DD"
  } else {
    state_t <- "undef"
  }
  return(state_t)
}

# Calculate percentage of R observations which are truncated by T
get_pctRbar <- function(dat, rbarnum){
  
  # Assumes any censoring/truncation is due to T
  pctRbar <- group_by(dat, World) %>% 
    summarise(mean(delta_R == 0, na.rm = TRUE))
  
  # Add in variable containing the numerical "value" of Rbar (for plotting only)
  pctRbar$rbarnum <- rbarnum
  
  # Return
  return(pctRbar)  
}

# Loop over data points and evaluate principal state membership
# and individual-level causal effects at each time point
calc_states_diffs <- function(dat, maxt) {
  
  # Grid of points to evaluate
  xts <- seq(0, maxt, length.out = nxts)

  # Shells
  n <- nrow(dat)/2
  states <- diffs <- matrix(NA, nrow = n, ncol = length(xts))

  # Process  
  T0 <- dat$T[dat$z == 0]
  T1 <- dat$T[dat$z == 1]
  R0 <- dat$R[dat$z == 0]
  R1 <- dat$R[dat$z == 1]
  delta_R0 <- dat$delta_R[dat$z == 0]
  delta_R1 <- dat$delta_R[dat$z == 1]
  R0_eff <- ifelse(delta_R0 == 1, R0, maxt + 1)
  R1_eff <- ifelse(delta_R1 == 1, R1, maxt + 1)
  alive0 <- t(sapply(T0,  FUN = function(x) {(x > xts)*1}))
  alive1 <- t(sapply(T1,  FUN = function(x) {(x > xts)*1}))
  rehosp0 <- 1 - t(sapply(R0_eff,  FUN = function(x) {(x > xts)*1}))
  rehosp1 <- 1 - t(sapply(R1_eff,  FUN = function(x) {(x > xts)*1}))

  # Calculate principal state and individual-level causal effect
  diffs  <- (rehosp1 - rehosp0) * (alive0 * alive1)
  for (id in 1:n) {
  states[id, ] <- sapply(xts, FUN = get_state, T0 = T0[id], T1 = T1[id])
  }
  
  # Return list
  return(list(states = states, diffs = diffs))
}

# Constructing composition plot
make_cplot <- function(dat, maxt) {
  
  # Get state matrix
  states <- calc_states_diffs(dat, maxt)[["states"]]
  
  # Grid of points to evaluate
  xts <- seq(0, maxt, length.out = nxts)
  
  # Preprocessing
  comp_dat0 <- data.frame(time = xts)  
  pstates <- c("AA", "CK", "TK", "DD")
  for (pstate in pstates) {
    comp_dat0[[pstate]] <- colMeans(states == pstate)
  }
  
  # Data frame containing proportion in each state at each t
  comp_dat <- data.frame(Time = rep(xts, times = length(pstates)),
                         State = rep(pstates, each = length(xts)))
  
  comp_dat$Proportion <- c(comp_dat0[["AA"]], comp_dat0[["CK"]],
                           comp_dat0[["TK"]], comp_dat0[["DD"]])
  comp_dat$State <- factor(comp_dat$State, levels = rev(pstates), 
                           ordered = TRUE)

  # Make the plot
  cplot <- ggplot(comp_dat, aes(x = Time, y = Proportion, fill = State)) + 
    geom_area(alpha = 0.6, color = "black") + 
    scale_x_continuous(breaks = seq(0, maxt, by = 30)) +
    ggtitle("Principal state composition over time") + 
    theme(legend.position = "bottom", legend.key.size = unit(0.6,"line"))
  
  # Return
  return(cplot)
}

# Constructing causal effect (alpha(t,t)) plot
make_aplot <- function(dat, maxt) {
  
  # Grid of points to evaluate
  xts <- seq(0, maxt, length.out = nxts)
  
  # Get state and difference matrices
  res <- calc_states_diffs(dat, maxt)

  # Subset to always-alive group
  diffs_AA <- res$diffs
  diffs_AA[res$states != "AA"] <- NA
  countAA <- colSums(res$states == "AA")
  newmaxt <- ifelse(any(countAA == 0), xts[min(which(countAA == 0))], maxt)
  
  # Make data frame of causal effects
  ce_dat <- data.frame(Time = xts,
                       alpha_tt = colMeans(diffs_AA, na.rm = TRUE))
  
  # Make plot
  aplot <- ggplot(ce_dat, aes(x = Time, y = alpha_tt)) +
    geom_smooth(method = "loess", se = FALSE) + 
    ylab(expression(alpha(r,t)~with~r==t)) + 
    scale_x_continuous(breaks = seq(0, newmaxt, by = 30)) +
    ggtitle("Causal effect on rehospitalization among Always Alive at t") +
    labs(subtitle = expression(alpha(r,t) == P(R[1]<r~'|'~T[0] > t, T[1] > t)-
                                 P(R[0]<r~'|'~T[0] > t, T[1] > t)))
  
  # Return
  return(aplot)
  
}



# Plotting functions ------------------------------------------------------

# Nothing like a function stolen from StackOverflow!
# Function puts a common legend under (or next to) multiple ggplot objects
grid_arrange_shared_legend <- function(..., nrow = 1, ncol = length(list(...)), position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  gl <- c(gl, nrow = nrow, ncol = ncol)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  
}

# Produce R plot object
make_rplot <- function(dat, showRbar) {

  # Calculate percentage of R values truncated by T
  pctRbar <- group_by(dat, World) %>%
    summarise(mean(delta_R == 0, na.rm = TRUE))
  pctRbar$rbarnum <- rbarnum

  # Construct r plot
  rplot <- ggplot(dat, aes(Rtoplot, fill = World, colour = World)) +
    geom_density(alpha = 0.2, trim = maxt) +
    ylab("Density") +
    ggtitle("Rehospitalization")
  if (showRbar) {
    rplot <- rplot  +
      scale_x_continuous("Rehospitalization time",
                         limits = c(0, maxtplot),
                         breaks = c(seq(0, maxt, by = 30), rbarnum),
                         labels = c(seq(0, maxt, by = 30), expression(bold(bar(R))))) +
      annotate("segment", color = hue_pal()(2)[which.max(pctRbar[[2]])], 
               size = 1.5,
               x = rbarnum, xend = rbarnum, 
               y = 0, yend = max(pctRbar[[2]])) +
      annotate("segment", color = hue_pal()(2)[which.min(pctRbar[[2]])], 
               size = 1.5,
               x = rbarnum, xend = rbarnum, 
               y = 0, yend = min(pctRbar[[2]]))
  } else {
    rplot <- rplot + 
      scale_x_continuous("Rehospitalization time",
                         limits = c(0, maxtplot),
                         breaks = seq(0, maxt, by = 30),
                         labels = seq(0, maxt, by = 30))
  }
  
  return(rplot)
}

# Make T plot object
make_tplot <- function(dat) {
  tplot <- ggplot(dat, aes(T, fill = World, colour = World)) +
    geom_density(alpha = 0.2) +
    ylab("Density") + 
    ggtitle("Death") +
    scale_x_continuous("Death time",
                       limits = c(0, maxt),
                       breaks = seq(0, maxt, by = 30),
                       labels = seq(0, maxt, by = 30))
  return(tplot)
}

# Make combined R and T plot object
make_rt_comb_plot <- function(dat, showRbar) {
  rplot <- make_rplot(dat, showRbar)
  tplot <- make_tplot(dat)
  suppressWarnings(grid_arrange_shared_legend(rplot, tplot, nrow = 1, ncol = 2))
}

# ui and server functions -------------------------------------------------


# Define UI for random distribution app ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Semicompeting Risks Data Simulations"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      withMathJax(),
      
      # Show Rbar values in R dist plot
      h5("\\(P(R_z = \\bar{\\mathbb{R}})\\)"),
      checkboxInput("showRbar", label = "Show", value = TRUE),
      
      # Simulate data with subject specific frailties?
      h5("Simulate subject-specific frailty (if not, \\(\\gamma_i = 1 \\ \\forall \\ i\\))"),
      checkboxInput("frailty", label = "Include frailty", value = TRUE),
      
      # Input: Slider for the number of observations to generate ----
      sliderInput("n",
                  "Sample size (smaller = faster):",
                  value = 10000,
                  min = 1000,
                  max = 50000,
                  step = 1000),
      sliderInput("theta.true",
                  "True \\(\\theta\\):",
                  value = 1,
                  min = 0.01,
                  max = 1),
      sliderInput("alpha1.true",
                  "True \\(\\alpha_1\\):",
                  value = 0.8,
                  min = 0.01,
                  max = 2),
      sliderInput("alpha2.true",
                  "True \\(\\alpha_2\\):",
                  value = 1.2,
                  min = 0.01,
                  max = 2),
      sliderInput("alpha3.true",
                  "True \\(\\alpha_3\\):",
                  value = 0.15,
                  min = 1.2,
                  max = 2),
      sliderInput("beta1.true",
                  "True \\(\\beta_1\\):",
                  value = -0.6,
                  min = -1,
                  max = 1, 
                  step = 0.05),
      sliderInput("beta2.true",
                  "True \\(\\beta_2\\):",
                  value = -0.8,
                  min = -1,
                  max = 1, 
                  step = 0.05),
      sliderInput("beta3.true",
                  "True \\(\\beta_3\\):",
                  value = -1,
                  min = -1,
                  max = 1, 
                  step = 0.05),
      sliderInput("kappa1.true",
                  "True \\(\\kappa_1\\):",
                  value = 0.01,
                  min = 0.01,
                  max = 2),
      sliderInput("kappa2.true",
                  "True \\(\\kappa_2\\):",
                  value = 0.1,
                  min = 0.01,
                  max = 2),
      sliderInput("kappa3.true",
                  "True \\(\\kappa_3\\):",
                  value = 0.1,
                  min = 0.01,
                  max = 2)
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: plots
      tabsetPanel(type = "tabs",
                  tabPanel("Notation", uiOutput("notation")),
                  tabPanel("Event distributions", plotOutput("eventPlots")),
                  tabPanel("Principal State Composition", plotOutput("compPlot")),
                  tabPanel("Causal Effect", plotOutput("causalEffectPlot"))
      )
      
    )
  )
)

# Define server logic for random distribution app ----
server <- function(input, output) {
  
  # Simulate new data when parameter inputs change
  dat <- reactive({
    simulate_RI(input$n, input$theta.true, 
                input$alpha1.true, input$alpha2.true, input$alpha3.true,
                input$beta1.true, input$beta2.true, input$beta3.true,
                input$kappa1.true, input$kappa2.true, input$kappa3.true,
                cens, frailty = input$frailty)
  })
  
  # Write model notation
  output$notation <- renderUI({
    includeMarkdown("Strategy1models.md")
  })
  
  # Show plots of nonterminal and terminal events ----
  output$eventPlots <- renderPlot({
    make_rt_comb_plot(dat(), input$showRbar)
  })
  
  # Show composition of population by principal state over time ----
  output$compPlot <- renderPlot({
    make_cplot(dat(), maxt)
  })
  
  # Graph causal effect over time
  output$causalEffectPlot <- renderPlot({
    make_aplot(dat(), maxt)
  })
 
  # Testing tab for printing stuff
  output$garbage <- renderText({
    testing  
  })
}

# Create Shiny app ----
shinyApp(ui, server)