#Power simulation for LMER model, code adapted from Jungerius (https://cjungerius.github.io/powersim/)

if (!require(pacman)) install.packages("pacman")

pacman::p_load(tidyverse, lmerTest, broom.mixed, styler)


# create dataset function
my_sim_data <- function(n_A, n_B,
                        beta_0, beta_g, beta_i, beta_gi,
                        tau_0, sigma,
                        n_trials_per_instance = 2,
                        seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
 #overall N 
  n_participants <- n_A + n_B
  
  # make participants and assign groups
  base <- tibble::tibble(
    participant = sprintf("P%04d", seq_len(n_participants)),
    group = sample(rep(c("A","B"), times = c(n_A, n_B)))
  )
  
  # expand to instances (T1,T2) and trials per instance
  dat <- base |>
    tidyr::expand_grid(instance = c("T1","T2")) |>
    tidyr::uncount(weights = n_trials_per_instance, .id = "trial") |>
    # ensure stable term names
    dplyr::mutate(
      group = factor(group, levels = c("A","B")),
      instance = factor(instance, levels = c("T1","T2"))
    ) |>
    # participant random intercept
    dplyr::group_by(participant) |>
    dplyr::mutate(u0 = rnorm(1, 0, tau_0)) |>
    dplyr::ungroup() |>
    # residual per observation
    dplyr::mutate(eps = rnorm(dplyr::n(), 0, sigma)) |>
    # generate score using treatment coding A/T1 as reference
    dplyr::mutate(
      g  = as.numeric(group == "B"),
      i  = as.numeric(instance == "T2"),
      score = beta_0 + beta_g*g + beta_i*i + beta_gi*g*i + u0 + eps
    )
  
  dat
}

single_run <- function(filename = NULL, ..., seed_offset = 0L) {
  dat_sim <- my_sim_data(..., seed = seed_offset)
  
  # run lmer and capture warnings
  ww <- ""
  suppressMessages(suppressWarnings(
    mod_sim <- withCallingHandlers({
      lmerTest::lmer(
        score ~ group * instance + (1 | participant),
        dat_sim, REML = FALSE
      )
    },
    warning = function(w) { ww <<- w$message })
  ))
  
  # tidy results
  sim_results <- broom.mixed::tidy(mod_sim) |>
    dplyr::mutate(warnings = ww)
  
  params <- list(...)
  for (name in names(params)) sim_results[[name]] <- params[[name]]
  
  if (!is.null(filename)) {
    append <- file.exists(filename)
    readr::write_csv(sim_results, filename, append = append)
  }
  sim_results
}

set.seed(123)

nreps <- 1000
#parameters based on the UK Biobank results
params <- tibble::tibble(rep = 1:nreps) |>
  dplyr::mutate(
    n_A = 96, 
    n_B = 36, 
    n_trials_per_instance = 1,
    beta_0 = 6.31, 
    beta_g = -0.02, 
    beta_i = -0.23, 
    beta_gi = 0.31,
    tau_0 = 1.74, 
    sigma = 1.64
  ) |>
  dplyr::select(-rep)

# run simulation
sims <- purrr::pmap_dfr(
  .l = dplyr::mutate(params ,seed_offset = 1000 + dplyr::row_number()),
  .f = single_run,
  filename = "sim_power.csv"
)

alpha <- 0.05
#output of the simulation
sims |>
  dplyr::filter(effect == "fixed", term == "groupB:instanceT2") |>
  dplyr::summarize(
    mean_estimate = mean(estimate),
    mean_se       = mean(std.error),
    power         = mean(p.value < alpha),
    .groups = "drop"
  )

