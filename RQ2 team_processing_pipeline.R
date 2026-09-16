# pipeline for RQ2 
# ============================ 1. [in team code] ============================
  ## load helper function to extract estimates for SDMA ##
  source(here::here("team_processing_helpers.R"))

# ============================ 2. [in team code] ============================
## extract estimates ##
teamXX_RQage_estimates <- bind_rows(
  
  # Main age effect on looking-time level (averaged across trials)
  extract_sdma_estimate(
    model = model_age,
    team_id = "team_XX",
    model_id = "model_age",
    estimand = "age_main_effect_on_looking_time",
    term = "age_days_c",                  # must match Parameter name in model_parameters()
    predictor = "age_days",
    contrast = NA_character_,
    model_role = "primary"
  ),
  
  # Age × trial interaction (age effect on habituation trajectory)
  extract_sdma_estimate(
    model = model_age,
    team_id = "team_XX",
    model_id = "model_age",
    estimand = "age_x_trial_effect_on_habituation",
    term = "age_days_c:trial_number_c",   # must match Parameter name
    predictor = "age_days",
    contrast = "interaction_with_trial",
    model_role = "primary"
  )
)
# ============================ 3. [in team code] ============================
# team level meta data and saving 

## add team specific analytic decisions
teamXX_RQage_estimates <- teamXX_RQage_estimates |>
  mutate(
    outcome_scale = "log(looking_time + 1)",
    age_scale = "age_days centered within sample",
    trial_scale = "trial_number centered within sample",
    covariates = NA_character_,  # or list them if you have any
    random_effects = "random intercepts for participant, trial_type, lab",
    analytic_sample = "complete cases used by fitted model",
    notes = "Standardized using model_parameters(..., standardize = 'refit')"
  )

## save per-team file
dir.create(
  here::here("RQage_sdma_estimates"),
  recursive = TRUE,
  showWarnings = FALSE
)

write_csv(
  teamXX_RQage_estimates,
  here::here(
    "RQage_sdma_estimates",
    "teamXX_RQage_estimates.csv"
  )
)

# ============================ 5. [in separate script] ============================ 
## combine team estimates ##

library(here)
library(readr)
library(dplyr)
library(purrr)

sdma_files <- list.files(
  here::here("RQage_sdma_estimates"),
  pattern = "_RQage_estimates\\.csv$",
  full.names = TRUE
)

sdma_all_teams_age <- sdma_files |>
  map_dfr(read_csv, show_col_types = FALSE)

write_csv(
  sdma_all_teams_age,
  here::here(
    "RQage_sdma_estimates",
    "sdma_all_teams_age.csv"
  )
)