library(tidyverse) # Efficient data manipulation and plotting


clargs <- commandArgs(trailingOnly = TRUE)
if (length(clargs) > 0) {
  datFile <- clargs[1]
} else {
  datFile <- "dat_0001.rds"
}


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
           W = pars$W, z0 = pars$z0, indVar = pars$indVar, sigma = pars$sigma[species],
           h2 = pars$h2[species], replicate = pars$replicate, rseedEps = pars$rseedEps)
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


read_rds(datFile) %>%
  mutate(func = case_when(
    model == "pairwise" ~ list(pairwise),
    model == "evoHOI" ~ list(evoHOI),
    model == "evoNP" ~ list(evoHOI_nopairwise),
    model == "hier" ~ list(hier),
    model == "hierHOI" ~ list(hierHOI),
    TRUE ~ list(randHOI)
  )) %>%
  mutate(pars = map2(pars, func, ~{ pars[[1]]$func = func[[1]]; pars[[1]] })) %>%
  mutate(sol = map2(pars, ic, function(pars, ic) {
    pars$eps <- assignRandomHOI(pars$model, pars$rseedEps, pars$S) # Assign random HOIs
    solveSteady(pars, ic, stol = 1e-15) # Integrate equations
  } )) %>%
  unnest(sol) %>%
  select(-pars, -time) %>%
  write_rds(datFile)
