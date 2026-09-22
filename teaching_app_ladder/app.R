## WODIN-style teaching app for the trachoma "SIS ladder" model
## (partial immunity + irreversible disease progression, no age structure).
##
## Local prototype built with Shiny (no Node.js/Docker/odin.api required).
## Sibling of ../teaching_app (the basic 2-compartment SIS model) - see
## that folder's app.R for the simpler starting point this was built from.
##
## Run with:  shiny::runApp("teaching_app_ladder", launch.browser = TRUE)

library(shiny)
library(ggplot2)
library(markdown)  # required by includeMarkdown() - explicit so rsconnect bundles it

## The 'ragg' graphics device DLL on this machine is incompatible with the
## installed R version (built-under-4.3.3 packages vs. R 4.0.5 runtime) -
## fall back to the standard png/Cairo device shiny uses when ragg is absent.
options(shiny.useragg = FALSE)

source(file.path("R", "ladder_model.R"))

`%+%` <- function(a, b) paste0(a, b)

ui <- fluidPage(
  titlePanel("Trachoma teaching app: SIS ladder model (partial immunity)"),
  tabsetPanel(
    tabPanel(
      "Simulator",
      sidebarLayout(
        sidebarPanel(
          h4("Transmission"),
          sliderInput("beta", "Transmission intensity, \u03b2 (per year)",
                      min = 0, max = 10, value = 3, step = 0.1),
          sliderInput("lExponent", "Infectiousness decay across rungs, l-exponent",
                      min = 0, max = 2, value = 0.05277, step = 0.01),
          hr(),
          h4("Recovery & partial immunity"),
          numericInput("n_rungs", "Number of infection-history rungs", value = 5, min = 1, max = 20),
          sliderInput("nu0", "Recovery rate at rung 1, \u03bd0 (per year)",
                      min = 0.1, max = 10, value = 0.885594, step = 0.05),
          sliderInput("nu1", "Asymptotic recovery rate at high rungs, \u03bd1 (per year)",
                      min = 0.1, max = 10, value = 4.208717, step = 0.05),
          sliderInput("nuExponent", "Speed of \u03bd0\u2192\u03bd1 transition, nu-exponent",
                      min = 0, max = 3, value = 0.751195, step = 0.01),
          hr(),
          h4("Irreversible disease progression"),
          sliderInput("r_TS", "Scarring accumulation rate, r_TS (per year)",
                      min = 0, max = 1, value = 0.10, step = 0.01),
          sliderInput("r_TT", "Progression TS \u2192 TT, r_TT (per year)",
                      min = 0, max = 1, value = 0.02, step = 0.01),
          sliderInput("alpha", "TT resolution rate (e.g. surgery), \u03b1 (per year)",
                      min = 0, max = 1, value = 0.04, step = 0.01),
          hr(),
          h4("Demography & population"),
          numericInput("life_expectancy", "Life expectancy (years)", value = 55, min = 5),
          numericInput("N", "Population size, N", value = 1000, min = 10),
          sliderInput("I0_prop", "Initial prevalence (in rung 1, I_1)",
                      min = 0, max = 1, value = 0.05, step = 0.01),
          hr(),
          h4("Mass Drug Administration (MDA)"),
          sliderInput("mda_coverage", "Coverage per round",
                      min = 0, max = 1, value = 0, step = 0.05),
          numericInput("mda_start", "Year of first round", value = 20, min = 0),
          numericInput("mda_interval", "Years between rounds", value = 1, min = 0.1),
          hr(),
          numericInput("t_end", "Simulation horizon (years)", value = 80, min = 1),
          hr(),
          uiOutput("r0_display")
        ),
        mainPanel(
          plotOutput("prevalence_plot", height = "280px"),
          plotOutput("compartment_plot", height = "260px"),
          plotOutput("rung_plot", height = "280px"),
          plotOutput("disease_plot", height = "260px")
        )
      )
    ),
    tabPanel(
      "About the model",
      withMathJax(includeMarkdown("explanation.md"))
    )
  )
)

server <- function(input, output, session) {

  sim <- reactive({
    run_ladder_model(
      beta = input$beta, nu0 = input$nu0, nu1 = input$nu1,
      nuExponent = input$nuExponent, lExponent = input$lExponent,
      n_rungs = input$n_rungs,
      r_TS = input$r_TS, r_TT = input$r_TT, alpha = input$alpha,
      life_expectancy = input$life_expectancy,
      N = input$N, I0_prop = input$I0_prop, t_end = input$t_end,
      mda_coverage = input$mda_coverage, mda_interval = input$mda_interval,
      mda_start = input$mda_start
    )
  })

  output$r0_display <- renderUI({
    df <- sim()
    tagList(
      strong(sprintf("R0 (invasion into a naive population) = %.2f",
                      attr(df, "R0_invasion"))),
      br(),
      em("Note: once partial immunity is established, effective "
         %+% "transmission depends on the mix of rungs, not just this number.")
    )
  })

  output$prevalence_plot <- renderPlot({
    df <- sim()
    mda_times <- attr(df, "mda_times")
    p <- ggplot(df, aes(time, prevalence_infection)) +
      geom_line(colour = "firebrick", linewidth = 1) +
      labs(x = "Time (years)", y = "Prevalence (sum of I rungs)/N",
           title = "Infection prevalence over time") +
      ylim(0, 1) +
      theme_minimal(base_size = 14)
    if (length(mda_times) > 0) {
      p <- p + geom_vline(xintercept = mda_times, linetype = "dashed",
                           colour = "steelblue", alpha = 0.6)
    }
    p
  })

  output$compartment_plot <- renderPlot({
    df <- sim()
    vars <- c("S_total", "I_total")
    df_long <- reshape(df[, c("time", vars)], direction = "long",
                        varying = list(vars), v.names = "count",
                        timevar = "compartment", times = vars)
    df_long$compartment <- factor(df_long$compartment, levels = vars)
    ggplot(df_long, aes(time, count, colour = compartment)) +
      geom_line(linewidth = 1) +
      labs(x = "Time (years)", y = "People",
           title = "Total susceptible vs total infected", colour = NULL) +
      theme_minimal(base_size = 14)
  })

  output$rung_plot <- renderPlot({
    df <- sim()
    n <- attr(df, "n_rungs")
    i_cols <- sprintf("I_%d", 1:n)
    df_long <- reshape(df[, c("time", i_cols)], direction = "long",
                        varying = list(i_cols), v.names = "count",
                        timevar = "rung", times = i_cols)
    df_long$rung <- factor(df_long$rung, levels = i_cols)
    df_long$proportion <- df_long$count / input$N
    ggplot(df_long, aes(time, proportion, fill = rung)) +
      geom_area(position = "stack") +
      scale_fill_manual(values = grDevices::hcl.colors(n, palette = "Dynamic"), name = NULL) +
      labs(x = "Time (years)", y = "Proportion of N",
           title = "Composition of infected rungs (I_1 = first infection ... I_n = most repeat infections)") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "bottom")
  })

  output$disease_plot <- renderPlot({
    df <- sim()
    vars <- c("TS", "TT")
    df_long <- reshape(df[, c("time", vars)], direction = "long",
                        varying = list(vars), v.names = "count",
                        timevar = "compartment", times = vars)
    df_long$compartment <- factor(df_long$compartment, levels = vars)
    ggplot(df_long, aes(time, count, colour = compartment)) +
      geom_line(linewidth = 1) +
      labs(x = "Time (years)", y = "People",
           title = "Irreversible disease burden: scarring (TS) and trichiasis (TT)",
           colour = NULL) +
      theme_minimal(base_size = 14)
  })
}

shinyApp(ui, server)
