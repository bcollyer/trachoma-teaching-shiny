## WODIN-style teaching app for a simple SIS trachoma model.
##
## Local prototype built with Shiny (no Node.js/Docker/odin.api required).
## Mimics the layout of the real mrc-ide/wodin teaching apps
## (https://epimodels.dide.ic.ac.uk/): a "Simulator" tab with parameter
## sliders + plots, and an "About the model" explanation tab.
##
## Run with:  shiny::runApp("teaching_app", launch.browser = TRUE)

library(shiny)
library(ggplot2)
library(markdown)  # required by includeMarkdown() - explicit so rsconnect bundles it

## The 'ragg' graphics device DLL on this machine is incompatible with the
## installed R version (built-under-4.3.3 packages vs. R 4.0.5 runtime) -
## fall back to the standard png/Cairo device shiny uses when ragg is absent.
options(shiny.useragg = FALSE)

source(file.path("R", "sis_model.R"))

ui <- fluidPage(
  titlePanel("Trachoma teaching app: simple SIS model"),
  tabsetPanel(
    tabPanel(
      "Simulator",
      sidebarLayout(
        sidebarPanel(
          h4("Transmission"),
          sliderInput("beta", "Transmission rate, \u03b2 (per year)",
                      min = 0, max = 10, value = 4, step = 0.1),
          sliderInput("gamma", "Recovery rate, \u03b3 (per year)",
                      min = 0.1, max = 10, value = 2, step = 0.1),
          sliderInput("I0_prop", "Initial prevalence",
                      min = 0, max = 1, value = 0.3, step = 0.01),
          numericInput("N", "Population size, N", value = 1000, min = 10),
          hr(),
          h4("Mass Drug Administration (MDA)"),
          sliderInput("mda_coverage", "Coverage per round",
                      min = 0, max = 1, value = 0, step = 0.05),
          numericInput("mda_start", "Year of first round", value = 2, min = 0),
          numericInput("mda_interval", "Years between rounds", value = 1, min = 0.1),
          hr(),
          numericInput("t_end", "Simulation horizon (years)", value = 20, min = 1),
          hr(),
          uiOutput("r0_display")
        ),
        mainPanel(
          plotOutput("prevalence_plot", height = "350px"),
          plotOutput("compartment_plot", height = "300px")
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
    run_sis_model(
      beta = input$beta,
      gamma = input$gamma,
      N = input$N,
      I0_prop = input$I0_prop,
      t_end = input$t_end,
      mda_coverage = input$mda_coverage,
      mda_interval = input$mda_interval,
      mda_start = input$mda_start
    )
  })

  output$r0_display <- renderUI({
    df <- sim()
    r0 <- attr(df, "R0")
    eq <- attr(df, "equilibrium_prevalence")
    tagList(
      strong(sprintf("R0 = %.2f", r0)),
      br(),
      sprintf("Equilibrium prevalence (no MDA): %.1f%%", 100 * eq)
    )
  })

  output$prevalence_plot <- renderPlot({
    df <- sim()
    mda_times <- attr(df, "mda_times")
    p <- ggplot(df, aes(time, prevalence)) +
      geom_line(colour = "firebrick", linewidth = 1) +
      labs(x = "Time (years)", y = "Prevalence (I/N)",
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
    df_long <- reshape(df[, c("time", "S", "I")], direction = "long",
                        varying = list(c("S", "I")), v.names = "count",
                        timevar = "compartment", times = c("S", "I"))
    ggplot(df_long, aes(time, count, colour = compartment)) +
      geom_line(linewidth = 1) +
      labs(x = "Time (years)", y = "People", title = "S and I compartments",
           colour = NULL) +
      theme_minimal(base_size = 14)
  })
}

shinyApp(ui, server)
