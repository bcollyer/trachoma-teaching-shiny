# Trachoma Teaching Apps (Shiny)

Interactive R Shiny apps for teaching trachoma (*Chlamydia trachomatis*)
transmission dynamics, mirroring (in local Shiny form) the WODIN-style
teaching apps in the `wodin-trachoma-site` companion repository.

Each folder is an independent, self-contained Shiny app (each can be
deployed separately, e.g. to https://www.shinyapps.io/):

- **`teaching_app/`** - the simplest model: a 2-compartment SIS model
  (S, I), no age structure, no partial immunity. MDA included via
  scheduled deSolve events.
- **`teaching_app_ladder/`** - adds an infection-history "ladder"
  (partial immunity via rung-dependent recovery/infectiousness) and
  irreversible disease progression (scarring, trichiasis), still no age
  structure.
- **`teaching_app_age_ladder/`** - the full model: adds the real
  18-age-group structure (age-assortative mixing, ageing, births,
  age-specific mortality) on top of the ladder mechanics, reusing the
  model code from `trachoma_first_principles_300_history.R` in the main
  modelling repository (vendored into `teaching_app_age_ladder/R/` so this
  app is self-contained for deployment).

## Running locally

From this folder:

```r
shiny::runApp("teaching_app", launch.browser = TRUE)
shiny::runApp("teaching_app_ladder", launch.browser = TRUE)
shiny::runApp("teaching_app_age_ladder", launch.browser = TRUE)
```

Each app has two tabs: **Simulator** (parameter sliders + plots) and
**About the model** (equations, diagrams, and biological rationale).

## Deploying to shinyapps.io

```r
install.packages("rsconnect")
rsconnect::setAccountInfo(name = "...", token = "...", secret = "...")  # from your shinyapps.io account dashboard
rsconnect::deployApp("teaching_app")
rsconnect::deployApp("teaching_app_ladder")
rsconnect::deployApp("teaching_app_age_ladder")
```

Each `deployApp()` call publishes that folder as its own app with its own
public URL.
