library(tidyverse) # Efficient data manipulation and plotting


Rcpp::sourceCpp("models.cpp") # Compile C++ functions implementing the various models

# Integrate ODE system for a given time sequence
solveModel <- function(pars, ic, tseq, ...) {
  deSolve::ode(func = pars$func, y = ic, parms = pars, times = tseq, ...) %>%
    organizeResults(pars)
}

# Integrate ODE system until its fixed point is reached
solveSteady <- function(pars, ic, ...) {
  rootSolve::runsteady(y = ic, time = c(0, Inf), func = pars$func,
                       parms = pars, ...)$y %>%
    enframe() %>%
    pivot_wider() %>%
    mutate(time = Inf, .before = 1) %>%
    organizeResults(pars)
}

# Put ODE solution results in tidy tibble
organizeResults <- function(sol, pars) {
  as_tibble(as.data.frame(sol)) %>%
    rename_with(~paste0("n_", 1:pars$S), 1 + 1:pars$S) %>%
    rename_with(~paste0("m_", 1:pars$S), 1 + pars$S + 1:pars$S) %>%
    pivot_longer(cols = !time, names_to = "variable", values_to = "v") %>%
    separate(col = variable, into = c("type", "species"), sep = "_", convert = TRUE) %>%
    pivot_wider(names_from = "type", values_from = "v") %>%
    mutate(species = as_factor(species), Sinit = pars$S, w = pars$w, theta = pars$theta,
           indVar = pars$indVar, sigma = pars$sigma[species], h2 = pars$h2[species],
           replicate = pars$replicate, rseedEps = pars$rseedEps)
}

# Generate random intraspecific standard deviations
randSigmas <- function(S, indVar) {
  if (indVar == "low") {
    runif(S, 0.003, 0.009)
  } else  if (indVar == "medium") {
    runif(S, 0.01, 0.03)
  } else {
    runif(S, 0.05, 0.1)
  }
}

# Generate higher-order interaction coefficients for intraHOI model
intraHOIcoefs <- function(S) {
  eps <- array(runif(S^3, 0, 0.04), c(S, S, S))
  for (i in 1:S) eps[i,i,] <- runif(S, 0.05, 0.1)
  eps
}

# Generate higher-order interaction coefficients for HOI model
HOIcoefs <- function(S) {
  array(runif(S^3, -0.01, 0.01), c(S, S, S))
}

# Assign HOI tensors for random HOI models, based on model name
assignRandomHOI <- function(model, rseedEps, S) {
  if (model == "intraHOI") {
    set.seed(rseedEps)
    intraHOIcoefs(S)
  } else if (model == "HOI") {
    set.seed(rseedEps)
    HOIcoefs(S)
  } else {
    NA_real_
  }
}

# Plot species abundances through time
plot_abundance <- function(dat) {
  ggplot(dat) +
    geom_line(aes(x = time, y = n, colour = as_factor(species))) +
    scale_y_continuous(name = "abundance", limits = c(0, NA)) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "none")
}

# Plot species trait values through time
plot_trait <- function(dat) {
  ggplot(dat) +
    geom_ribbon(aes(x = time, ymin = m - sigma, ymax = m + sigma,
                    fill = as_factor(species)), alpha = 0.15) +
    geom_line(aes(x = time, y = m, colour = as_factor(species))) +
    ylab("trait value") +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "none")
}

# Make a plot of the trait distributions at some given moment
plot_snapshot <- function(sol, bfun, limits = c(NA, NA, NA), res = 501) {
  # Normal distribution if the trait z is within 4 sigma of the mean, and NA otherwise:
  truncnorm <- function(z, m, s) if_else(abs(z - m) < 4 * s, dnorm(z, m, s), NA_real_)
  if (is.na(limits[1])) limits[1] <- min(sol$m - 4 * sol$sigma)
  if (is.na(limits[2])) limits[2] <- max(sol$m + 4 * sol$sigma)
  traitaxis <- seq(limits[1], limits[2], l = res)
  traits <- crossing(sol, trait = traitaxis) %>%
    mutate(density = n * truncnorm(trait, m, sigma))
  landscape <- tibble(trait = traitaxis) %>% # for plotting intrinsic rates
    mutate(r = bfun(trait) * if_else(
      is.na(limits[3]), max(traits$density, na.rm = TRUE), limits[3] * 0.98)) %>%
    mutate(r = if_else(r < 0, NA_real_, r)) # where b(z) < 0: set to NA
  ggplot(traits) +
    geom_line(aes(x = trait, y = density, colour = species), na.rm = TRUE) +
    geom_ribbon(aes(x = trait, ymin = 0, ymax = density, fill = species), alpha = 0.15) +
    geom_line(data = landscape, aes(x = trait, y = r), linetype = "dashed",
              colour = "darkred", alpha = 0.5, na.rm = TRUE) +
    scale_x_continuous(name = "trait value", limits = c(limits[1], limits[2])) +
    scale_y_continuous(name = "density", limits = c(0, limits[3])) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "none")
}

# Plot time series of abundances, time series of trait values, and
# snapshot of the trait distributions at time = moment
plot_all <- function(dat, bfun, moment = 0, limits = c(NA, NA, NA), res = 501) {
  cowplot::plot_grid(plot_abundance(dat),
                     plot_trait(dat),
                     plot_snapshot(filter(dat, time == moment), bfun, limits, res),
                     ncol = 1, align = "hv")
}


# Create table of parameters, solve equations, and save results
crossing( # Fully factorial combination of:
  S = 40, # Number of species
  theta = 0.5, # Width of intrinsic growth function
  replicate = 1:1, # Replicate number
  w = c(0.2), # Competition width
  indVar = c("medium"), # Degree of individual trait variation
  model = c("pairwise", "hier", "evoHOI", "hierHOI", "intraHOI", "HOI") # Which model?
) %>%
  mutate(func = case_when( # New column containing the actual function to be integrated
    model == "pairwise" ~ list(pairwise),
    model == "evoHOI" ~ list(evoHOI),
    model == "evoNP" ~ list(evoHOI_nopairwise),
    model == "hier" ~ list(hier),
    model == "hierHOI" ~ list(hierHOI),
    TRUE ~ list(randHOI)
  )) %>%
  mutate(sigma = map2(S, indVar, randSigmas)) %>% # Intraspecific trait std devs
  mutate(h2 = map(S, ~runif(.x, 0.1, 1.5))) %>% # Heritabilities (slightly randomized)
  mutate(kappa = if_else(model %in% c("evoHOI", "hierHOI"), 1, NA_real_)) %>%
  mutate(z0 = if_else(model %in% c("hierHOI"), 0.7, NA_real_)) %>%
  mutate(W = if_else(model %in% c("hier", "hierHOI"), 0.1, NA_real_)) %>%
  mutate(rseedEps = sample(100000:999999, nrow(.)),
         rseedEps = if_else(model %in% c("intraHOI", "HOI"), rseedEps, NA_integer_)) %>%
  mutate(pars = pmap(., list)) %>% # Merge all model parameters into one list
  # Initial conditions: initial densities = 1; initial trait means from [-theta, theta]:
  mutate(ic = pmap(list(S, theta, model), function(S, theta, model) {
    if (model %in% c("hier", "hierHOI")) {
      c(rep(1, times = S), runif(S, min = 0, max = 1.5))
    } else {
      c(rep(1, times = S), runif(S, min = -theta, max = theta))
    }
  })) %>%
  select(pars, ic, model) %>% # Keep only the necessary columns
  # Assign HOI tensors if needed, and solve equations:
  mutate(sol = map2(pars, ic, function(pars, ic) {
    pars$eps <- assignRandomHOI(pars$model, pars$rseedEps, pars$S) # Assign random HOIs
    solveModel(pars, ic, tseq = seq(0, 1e14, l = 1001)) # Integrate eqs
  } )) %>%
  unnest(sol) -> # Expand out the column "sol" containing the solution
  sol

sol %>%
  filter(n > 1e-6, time == max(time), model == "hierHOI") %>%
  print() %>%
  plot_snapshot(function(x) 1 - exp(-x), c(-0.05, NA, NA))

# sol %>%
#   filter(n > 1e-6, time %in% c(0, 1e14)) %>%
#   select(-pars, -ic, -replicate) %>%
#   write_rds("data_fig1.rds", compress = "xz")
