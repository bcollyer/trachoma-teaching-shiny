## WODIN-style teaching app for the FULL trachoma model: age structure
## (18 groups, from trachoma_first_principles_300_history.R) combined with
## the infection-history "ladder" (partial immunity, from teaching_app_ladder).
##
## Local prototype built with Shiny (no Node.js/Docker/odin.api required).
## Reuses the real model code in ../trachoma_first_principles_300_history.R
## via R/model_lib.R rather than reimplementing it.
##
## Run with:  shiny::runApp("teaching_app_age_ladder", launch.browser = TRUE)

library(shiny)
library(ggplot2)
library(markdown)  # required by includeMarkdown() - explicit so rsconnect bundles it

## The 'ragg' graphics device DLL on this machine is incompatible with the
## installed R version (built-under-4.3.3 packages vs. R 4.0.5 runtime) -
## fall back to the standard png/Cairo device shiny uses when ragg is absent.
options(shiny.useragg = FALSE)

source(file.path("R", "model_lib.R"))

ui <- fluidPage(
  titlePanel("Trachoma teaching app: age-structured ladder model"),
  tabsetPanel(
    tabPanel(
      "Simulator",
      sidebarLayout(
        sidebarPanel(
          h4("Transmission"),
          sliderInput("transRate", "Transmission intensity, \u03b2 (per year)",
                      min = 0, max = 10, value = 3, step = 0.1),
          sliderInput("epsilon", "Age-assortative mixing, \u03b5",
                      min = 0, max = 1, value = 0.5, step = 0.05),
          hr(),
          h4("Recovery & partial immunity"),
          numericInput("m_hist", "Number of infection-history rungs",
                       value = 5, min = 1, max = 20),
          sliderInput("nu0", "Recovery rate at rung 1, \u03bd0 (per year)",
                      min = 0.1, max = 10, value = 0.885594, step = 0.05),
          sliderInput("nu1", "Asymptotic recovery rate at high rungs, \u03bd1 (per year)",
                      min = 0.1, max = 10, value = 4.208717, step = 0.05),
          sliderInput("nuExponent", "Speed of \u03bd0\u2192\u03bd1 transition, nu-exponent",
                      min = 0, max = 3, value = 0.751195, step = 0.01),
          sliderInput("lExponent", "Infectiousness decay across rungs, l-exponent",
                      min = 0, max = 2, value = 0.05277, step = 0.01),
          hr(),
          h4("Irreversible disease progression"),
          sliderInput("k_TS", "Scarring accumulation rate, k_TS (per year)",
                      min = 0, max = 1, value = 0.08, step = 0.01),
          sliderInput("k_TT", "Progression TS \u2192 TT, k_TT (per year)",
                      min = 0, max = 1, value = 0.02, step = 0.01),
          hr(),
          h4("Demography & simulation"),
          numericInput("alpha", "Birth rate, \u03b1", value = 0.04, min = 0, step = 0.005),
          numericInput("years", "Simulation horizon (years)", value = 100, min = 20, max = 300),
          numericInput("dt_out", "Plot time resolution (years per point)",
                       value = 0.1, min = 0.01, max = 5, step = 0.01),
          hr(),
          h4("Mass Drug Administration (MDA)"),
          sliderInput("mda_coverage", "Coverage per round",
                      min = 0, max = 1, value = 0, step = 0.05),
          numericInput("mda_start", "Year of first round", value = 40, min = 0),
          numericInput("mda_interval", "Years between rounds", value = 1, min = 0.1),
          sliderInput("drug_efficacy", "Drug efficacy", min = 0, max = 1, value = 1, step = 0.05),
          hr(),
          uiOutput("r0_display")
        ),
        mainPanel(
          plotOutput("timeseries_plot", height = "300px"),
          plotOutput("age_hist_infection", height = "300px"),
          plotOutput("age_hist_disease", height = "300px")
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
    run_age_ladder_model(
      transRate = input$transRate, epsilon = input$epsilon,
      nu0 = input$nu0, nu1 = input$nu1, nuExponent = input$nuExponent,
      lExponent = input$lExponent, alpha = input$alpha,
      k_TS = input$k_TS, k_TT = input$k_TT, m_hist = input$m_hist,
      years = input$years, dt_out = input$dt_out,
      mda_coverage = input$mda_coverage, mda_start = input$mda_start,
      mda_interval = input$mda_interval, drug_efficacy = input$drug_efficacy
    )
  })

  output$r0_display <- renderUI({
    r0 <- sim()$R0
    tagList(
      strong(sprintf("R0 (next-generation matrix) = %.2f", r0)),
      br(),
      em("Computed at the disease-free equilibrium age distribution.")
    )
  })

  output$timeseries_plot <- renderPlot({
    res <- sim()
    p <- ggplot(res$ts_df, aes(time, prevalence_infection)) +
      geom_line(colour = "firebrick", linewidth = 1) +
      labs(x = "Time (years)", y = "Overall infection prevalence",
           title = "Infection prevalence over time (all ages & rungs combined)") +
      ylim(0, NA) +
      theme_minimal(base_size = 14)
    if (length(res$mda_times) > 0) {
      p <- p + geom_vline(xintercept = res$mda_times, linetype = "dashed",
                           colour = "steelblue", alpha = 0.6)
    }
    p
  })

  output$age_hist_infection <- renderPlot({
    df <- sim()$age_df
    ggplot(df, aes(age_label, prevalence_infection)) +
      geom_col(fill = "#c0392b") +
      labs(x = "Age band (years)", y = "Prevalence",
           title = "Equilibrium infection prevalence by age") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })

  output$age_hist_disease <- renderPlot({
    df <- sim()$age_df
    df_long <- reshape(df[, c("age_label", "prevalence_TS", "prevalence_TT")],
                        direction = "long",
                        varying = list(c("prevalence_TS", "prevalence_TT")),
                        v.names = "prevalence",
                        timevar = "sequela", times = c("TS", "TT"))
    ggplot(df_long, aes(age_label, prevalence, fill = sequela)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c(TS = "#2E7D32", TT = "#E65100"), name = NULL) +
      labs(x = "Age band (years)", y = "Prevalence",
           title = "Equilibrium disease burden (TS, TT) by age") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top")
  })
}

shinyApp(ui, server)
