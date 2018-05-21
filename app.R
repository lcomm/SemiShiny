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
                        sigma.true, 
                        alpha1.true, alpha2.true, alpha3.true,
                        alpha4.true, alpha5.true, alpha6.true,
                        kappa1.true, kappa2.true, kappa3.true,
                        kappa4.true, kappa5.true, kappa6.true,
                        cens) {
  
  # Cens -> vector
  if (length(cens) != 2) {
    cens <- c(cens, cens + 1)
  }
  
  # Simulate potential confounder
  # (Not currently doing anything with the confounder x_c)
  x_c <- rnorm(n, mean = 0, sd = 1)
  a <- rep(0:1, each = n)
  
  # Frailty - draw 
  if (sigma.true > 0) {
    gamma.true <- rgamma(n, 1/sigma.true, 1/sigma.true)
  }
  if (sigma.true == 0) {
    gamma.true <- rep(1, n)
  }
  
  # Design matrices (shell placeholder covariate for now)
  x1 <- x2 <- x3 <- cbind(log(gamma.true), 0)
  ID <- c(1:n)
  
  # Simulate data under control and treatment
  # theta.true (sigma.true in our notation) is set to 0 here 
  # (we have already included frailties in the linear design matrix)
  sim_control <- simID(cluster = NULL, 
                       x1 = x1, x2 = x2, x3 = x3, 
                       alpha1.true = alpha1.true,
                       alpha2.true = alpha2.true,
                       alpha3.true = alpha3.true,
                       beta1.true = c(1, 0),
                       beta2.true = c(1, 0),
                       beta3.true = c(1, 0),
                       kappa1.true = kappa1.true,
                       kappa2.true = kappa2.true,
                       kappa3.true = kappa3.true,
                       theta.true = 0, 
                       SigmaV.true = NULL, 
                       cens = cens)
  
  sim_treated <- simID(cluster = NULL, 
                       x1 = x1, x2 = x2, x3 = x3, 
                       alpha1.true = alpha4.true,
                       alpha2.true = alpha5.true,
                       alpha3.true = alpha6.true,
                       beta1.true = c(1, 0),
                       beta2.true = c(1, 0),
                       beta3.true = c(1, 0),
                       kappa1.true = kappa4.true,
                       kappa2.true = kappa5.true,
                       kappa3.true = kappa6.true,
                       theta.true = 0, 
                       SigmaV.true = NULL, 
                       cens = cens)
  
  # Rename and add treatment column
  sim_dat <- rbind(sim_control, sim_treated)
  colnames(sim_dat) <- c("R", "delta_R", "T", "delta_T")
  sim_dat$ID <- rep(ID, 2)
  sim_dat$z <- rep(0:1, each = n)
  
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
    ggtitle("Hospital readmission")
  if (showRbar) {
    rplot <- rplot  +
      scale_x_continuous("Readmission time",
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
      scale_x_continuous("Readmission time",
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


# Constructing TV-SACE causal effect (Q(t,t)) plot
make_tvsace_plot <- function(dat, maxt, yrange) {
  
  # Grid of points to evaluate
  xts <- seq(0, maxt, length.out = nxts)
  
  # Get state matrix
  states <- calc_states(dat, xts)
  
  # Get difference in cumulative incidence of readmission by t
  # (Time-varying survivor average causal effect)
  diffs <- calc_readmission_diffs(dat, xts)
  diffs[states != "AA"] <- NA
  
  # Subset to always-alive group
  countAA <- colSums(states == "AA")
  newmaxt <- ifelse(any(countAA == 0), xts[min(which(countAA == 0))], maxt)
  
  # Make data frame of causal effects
  ce_dat <- data.frame(Time = xts, alpha_tt = colMeans(diffs, na.rm = TRUE),
                       NumAlive = countAA)
  
  # Make plot
  tvsaceplot <- ggplot(ce_dat, aes(x = Time, y = alpha_tt, color = NumAlive)) +
    geom_line(stat = "smooth", method = "loess", size = 1.2, alpha = 0.2) + 
    geom_point(size = 1.5) +
    ylab(expression(Q(r,t)~with~r==t)) + 
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
  return(tvsaceplot)
}


# ui and server functions -------------------------------------------------

# Define UI ----
ui <- fluidPage(
  theme = shinythemes::shinytheme("flatly"),
  
  # App title ----
  titlePanel("Semicompeting Risks Data Simulations"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      withMathJax(),
      
      ### Simulation parameters ----
      h3("Data Generation Parameters"),
      sliderInput("n",
                  "Sample size (smaller = app runs faster):",
                  value = 10000,
                  min = 5000,
                  max = 50000,
                  step = 5000),
      sliderInput("sigma.true",
                  "True \\(\\sigma\\), variance of frailties (0 = no heterogeneity):",
                  value = 0.5,
                  min = 0,
                  max = 2,
                  step = 0.1),
      sliderInput("alpha1.true",
                  "Model 1 control shape \\(\\alpha_1\\) (bigger = more clumping):",
                  value = 0.8,
                  min = 0.01,
                  max = 4),
      sliderInput("alpha2.true",
                  "Model 2 control shape \\(\\alpha_2\\) (bigger = more clumping):",
                  value = 1.2,
                  min = 0.1,
                  max = 2),
      sliderInput("alpha3.true",
                  "Model 3 control shape \\(\\alpha_3\\) (bigger = more clumping):",
                  value = 1.2,
                  min = 0.1,
                  max = 2),
      sliderInput("alpha4.true",
                  "Model 1 treated shape \\(\\alpha_4\\) (bigger = more clumping):",
                  value = 0.8,
                  min = 0.01,
                  max = 4),
      sliderInput("alpha5.true",
                  "Model 2 treated shape \\(\\alpha_5\\) (bigger = more clumping):",
                  value = 1.2,
                  min = 0.1,
                  max = 2),
      sliderInput("alpha6.true",
                  "Model 3 treated shape \\(\\alpha_6\\) (bigger = more clumping):",
                  value = 1.2,
                  min = 0.1,
                  max = 2),
      sliderInput("kappa1.true",
                  "Model 1 control reference hazard \\(\\kappa_1\\) (bigger = faster):",
                  value = 0.2,
                  min = 0.15,
                  max = 7.5, 
                  step = 0.05),
      sliderInput("kappa2.true",
                  "Model 2 control reference hazard \\(\\kappa_2\\) (bigger = faster):",
                  value = 0.1,
                  min = 0.15,
                  max = 7.5, 
                  step = 0.05),
      sliderInput("kappa3.true",
                  "Model 3 control reference hazard \\(\\kappa_3\\) (bigger = faster):",
                  value = 0.35,
                  min = 0.15,
                  max = 7, 
                  step = 0.05),
      sliderInput("kappa4.true",
                  "Model 1 treated reference hazard \\(\\kappa_4\\) (bigger = faster):",
                  value = 0.5,
                  min = 0.15,
                  max = 7.5, 
                  step = 0.05),
      sliderInput("kappa5.true",
                  "Model 2 treated reference hazard \\(\\kappa_5\\) (bigger = faster):",
                  value = 0.7,
                  min = 0.15,
                  max = 7.5, 
                  step = 0.05),
      sliderInput("kappa6.true",
                  "Model 3 treated reference hazard \\(\\kappa_6\\) (bigger = faster):",
                  value = 0.2,
                  min = 0.15,
                  max = 7, 
                  step = 0.05),
      
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
                  step = 0.05)
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: dgp description, density, composition, and effect plots
      tabsetPanel(type = "tabs",
                  tabPanel("Notation", uiOutput("notation")),
                  tabPanel("Event distributions", plotOutput("eventPlots")),
                  tabPanel("Principal State Composition", plotOutput("compPlot")),
                  tabPanel("Time-varying SACE", plotOutput("tvsacePlot"))
      )
      
    )
  )
)

# Define server logic for random distribution app ----
server <- function(input, output) {
  
  # Simulate new data when parameter inputs change
  dat <- reactive({
    simulate_RI(input$n, input$sigma.true, 
                input$alpha1.true, input$alpha2.true, input$alpha3.true,
                input$alpha4.true, input$alpha5.true, input$alpha6.true,
                input$kappa1.true, input$kappa2.true, input$kappa3.true,
                input$kappa4.true, input$kappa5.true, input$kappa6.true,
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
  
  # Graph TV-SACE effect over time
  output$tvsacePlot <- renderPlot({
    make_tvsace_plot(dat(), maxt, yrange = input$yrange)
  })
 
}


# Create Shiny app ----
shinyApp(ui, server)
