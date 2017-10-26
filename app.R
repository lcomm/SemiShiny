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
                        cens) {
  
  # Cens -> vector
  if (length(cens) != 2) {
    cens <- c(cens, cens + 1)
  }
  
  # Simulate potential confounder
  # (Not currently doing anything with the confounder x_c)
  x_c <- rnorm(n, mean = 0, sd = 1)
  a <- rep(0:1, each = n)
  
  # Frailty - draw and include as if it's a normal covariate with coef 1
  if (theta.true > 0) {
    log_gamma_true <- log(rgamma(n, 1/theta.true, 1/theta.true))
  }
  if (theta.true == 0) {
    log_gamma_true <- rep(0, n)
  }
  
  # Design matrices 
  x1 <- x2 <- x3 <- cbind(a, c(log_gamma_true, log_gamma_true))
  ID <- c(1:n, 1:n)
  
  # Simulate, adding one to each set of parameters
  # theta.true set to 0 here because we have already included the frailties in X
  sim_dat <- simID(cluster = NULL, x1, x2, x3, 
                   alpha1.true = alpha1.true,
                   alpha2.true = alpha2.true,
                   alpha3.true = alpha3.true,
                   beta1.true = c(beta1.true, 1),
                   beta2.true = c(beta2.true, 1),
                   beta3.true = c(beta3.true, 1),
                   kappa1.true = kappa1.true,
                   kappa2.true = kappa2.true,
                   kappa3.true = kappa3.true,
                   theta.true = 0, 
                   SigmaV.true = NULL, cens)
  
  # Rename and add treatment column
  colnames(sim_dat) <- c("R", "delta_R", "T", "delta_T")
  sim_dat$ID <- ID
  sim_dat$z <- as.factor(a)
  
  # Make factor version of Z for nicer plotting
  sim_dat$World <- factor(sim_dat$z, levels = c("1", "0"), 
                          labels = c("Treated (z=1)", "Control (z=0)"))
  
  # Add version of R for plotting (where Rbar has a value)
  sim_dat$Rtoplot <- sim_dat$R
  sim_dat$Rtoplot[sim_dat$delta_R == 0] <- rbarnum * 2
  
  # Return created data
  return(sim_dat)
}


# Make a matrix of whether elements of a vector are greater than times in a 
# time sequence vector
# Time are along the columns and rows correspond to elements of vec
# Rows are time-varying T/F indicators for survival past t_j
check_vector_gt_timeseq <- function(vec, timeseq) {
  return(matrix((vec > rep(timeseq, each = length(vec))), 
                nrow = length(vec), ncol = length(timeseq)))
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


# Calculate matrix of principal states over all time points
calc_states <- function(dat, xts) {
  
  # Shells
  n <- nrow(dat)/2
  states <- matrix(NA, nrow = n, ncol = length(xts))

  # Process  
  T0 <- dat$T[dat$z == 0]
  T1 <- dat$T[dat$z == 1]
  
  # Matrix of T/F for person i still alive at time j
  alive0 <- check_vector_gt_timeseq(vec = T0, timeseq = xts)
  alive1 <- check_vector_gt_timeseq(vec = T1, timeseq = xts)
  
  # Translate into survival principal states
  states <- matrix(NA, nrow = n, ncol = length(xts))
  states[(alive0 == 1) & (alive1 == 1)] <- "AA"
  states[(alive0 == 1) & (alive1 == 0)] <- "TK"
  states[(alive0 == 0) & (alive1 == 1)] <- "CK"
  states[(alive0 == 0) & (alive1 == 0)] <- "DD"
  
  # Return state matrix
  return(states)
}


# Calculate difference in readmission incidence by each t in a vector of times
calc_readmission_diffs <- function(dat, xts) { 
  
  # Make shell
  diffs <- matrix(NA, nrow = nrow(dat)/2, ncol = length(xts))
  
  # Process data set
  R0 <- dat$R[dat$z == 0]
  R1 <- dat$R[dat$z == 1]
  delta_R0 <- dat$delta_R[dat$z == 0]
  delta_R1 <- dat$delta_R[dat$z == 1]
  R0_eff <- ifelse(delta_R0 == 1, R0, max(xts) + 1)
  R1_eff <- ifelse(delta_R1 == 1, R1, max(xts) + 1)

  # Matrix of T/F whether readmission occurred for person i by time j
  readmitted0 <- 1 - check_vector_gt_timeseq(vec = R0_eff, timeseq = xts)
  readmitted1 <- 1 - check_vector_gt_timeseq(vec = R1_eff, timeseq = xts)
  
  # Calculate hospitalization incidence difference
  # Will later perform calculations only on AA
  diffs <- readmitted1 - readmitted0
  
  # Return list
  return(diffs)
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


# Produce R (rehospitalization) plot object
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


# Make T (death time) plot object
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


# Constructing composition plot
make_cplot <- function(dat, maxt) {
  
  # Grid of points to evaluate
  xts <- seq(0, maxt, length.out = nxts)
  
  # Get state matrix
  states <- calc_states(dat, xts)
  
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
    theme(legend.position = "bottom")
  
  # Return
  return(cplot)
}


# Constructing causal effect (alpha(t,t)) plot
make_aplot <- function(dat, maxt, yrange) {
  
  # Grid of points to evaluate
  xts <- seq(0, maxt, length.out = nxts)
  
  # Get state matrix
  states <- calc_states(dat, xts)
  
  # Get difference in cumulative incidence of readmission by t
  diffs <- calc_readmission_diffs(dat, xts)
  diffs[states != "AA"] <- NA
  
  # Subset to always-alive group
  countAA <- colSums(states == "AA")
  newmaxt <- ifelse(any(countAA == 0), xts[min(which(countAA == 0))], maxt)
  
  # Make data frame of causal effects
  ce_dat <- data.frame(Time = xts, alpha_tt = colMeans(diffs, na.rm = TRUE),
                       NumAlive = countAA)
  
  # Make plot
  aplot <- ggplot(ce_dat, aes(x = Time, y = alpha_tt, color = NumAlive)) +
    geom_line(stat = "smooth", method = "loess", size = 1.2, alpha = 0.2) + 
    geom_point(size = 1.5) +
    ylab(expression(alpha(r,t)~with~r==t)) + 
    ylim(yrange[1], yrange[2]) + 
    scale_x_continuous(breaks = seq(0, newmaxt, by = 30)) +
    scale_color_gradientn(colors = hue_pal()(7), 
                          limits = c(0, max(ce_dat$NumAlive)),
                          guide = guide_colorbar("Number Always-Alive")) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    ggtitle("Causal effect on rehospitalization among Always Alive at t") +
    labs(subtitle = expression(alpha(r,t) == P(R[1]<r~'|'~T[0] > t, T[1] > t)-
                                 P(R[0]<r~'|'~T[0] > t, T[1] > t)))
  
  # Return
  return(aplot)
}



# ui and server functions -------------------------------------------------

# Define UI ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("Semicompeting Risks Data Simulations"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      withMathJax(),
      
      ### Plot parameters ----
      h3("Plot Parameters"),
      
      # Show Rbar values in R dist plot
      h5("\\(P(R_z = \\bar{\\mathbb{R}})\\)"),
      checkboxInput("showRbar", label = "Show", value = TRUE),
      
      # Control range for causal effect graph
      sliderInput("yrange",
                  "Causal effect Y axis limits:",
                  value = c(-0.3, 0.3),
                  min = -0.5,
                  max = 0.5, 
                  step = 0.05),
      
      ### Simulation parameters ----
      h3("Data Generation Parameters"),
      sliderInput("n",
                  "Sample size (smaller = app runs faster):",
                  value = 10000,
                  min = 5000,
                  max = 50000,
                  step = 5000),
      sliderInput("theta.true",
                  "True \\(\\theta\\) (0 = no induced correlation):",
                  value = 0.5,
                  min = 0,
                  max = 2,
                  step = 0.1),
      sliderInput("alpha1.true",
                  "True \\(\\alpha_1\\) (bigger = more clumping):",
                  value = 0.8,
                  min = 0.01,
                  max = 4),
      sliderInput("alpha2.true",
                  "True \\(\\alpha_2\\) (bigger = more clumping):",
                  value = 1.2,
                  min = 0.1,
                  max = 2),
      sliderInput("alpha3.true",
                  "True \\(\\alpha_3\\) (bigger = more clumping):",
                  value = 1.2,
                  min = 0.1,
                  max = 2),
      sliderInput("beta1.true",
                  "True \\(\\beta_1\\) (positive = treatment accelerates):",
                  value = -0.6,
                  min = -2,
                  max = 2, 
                  step = 0.2),
      sliderInput("beta2.true",
                  "True \\(\\beta_2\\) (positive = treatment accelerates):",
                  value = -0.8,
                  min = -2,
                  max = 2, 
                  step = 0.2),
      sliderInput("beta3.true",
                  "True \\(\\beta_3\\) (positive = treatment accelerates):",
                  value = -1,
                  min = -2,
                  max = 2, 
                  step = 0.2),
      sliderInput("logkappa1.true",
                  "True \\(\\log(\\kappa_1) = \\beta_{01}\\) (smaller = smaller rate):",
                  value = -2,
                  min = -8.5,
                  max = 0, 
                  step = 0.1),
      sliderInput("logkappa2.true",
                  "True \\(\\log(\\kappa_2) = \\beta_{02}\\) (smaller = smaller rate):",
                  value = -6.5,
                  min = -8.5,
                  max = 0, 
                  step = 0.1),
      sliderInput("logkappa3.true",
                  "True \\(\\log(\\kappa_3) = \\beta_{03}\\) (smaller = smaller rate):",
                  value = -5,
                  min = -8.5,
                  max = 0, 
                  step = 0.1)
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: dgp description, density, composition, and effect plots
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
                exp(input$logkappa1.true), 
                exp(input$logkappa2.true), 
                exp(input$logkappa3.true),
                cens)
  })
  
  # Write model notation
  output$notation <- renderUI({
    includeMarkdown("Strategy1models.md")
  })
  
  # Show plots of nonterminal and terminal events
  output$eventPlots <- renderPlot({
    make_rt_comb_plot(dat(), input$showRbar)
  })
  
  # Show composition of population by principal state over time
  output$compPlot <- renderPlot({
    make_cplot(dat(), maxt)
  })
  
  # Graph causal effect over time
  output$causalEffectPlot <- renderPlot({
    make_aplot(dat(), maxt, yrange = input$yrange)
  })
 
}


# Create Shiny app ----
shinyApp(ui, server)
