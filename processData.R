library(tidyverse)

invSimpson <- function(n) {
  p <- n / sum(n)
  1 / sum(p^2)
}

traitCV <- function(m) {
  traitdiff <- diff(sort(m))
  sd(traitdiff) / mean(traitdiff)
}

nullCVs <- function(S, replicates) {
  runif(S * replicates) %>%
    matrix(replicates, S) %>%
    apply(1, traitCV)
}

traitPval <- function(m, replicates) {
  mtrans <- (m - min(m)) / (max(m) - min(m)) # Squeeze traits into [0, 1] range
  nulls <- nullCVs(length(mtrans), replicates)
  observed <- traitCV(mtrans)
  length(nulls[nulls < observed]) / replicates
}

intraHOIcoefs <- function(S) {
  eps <- array(runif(S^3, 0, 0.04), c(S, S, S))
  for (i in 1:S) eps[i,i,] <- runif(S, 0.05, 0.1)
  eps
}

HOIcoefs <- function(S) {
  array(runif(S^3, -0.01, 0.01), c(S, S, S))
}

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

sigmoid <- function(x) pnorm(x * sqrt(2)) # Hierarchical competition kernel

epsilonEvoHOI <- function(m, v, w) {
  S <- length(m)
  eps <- array(0, c(S, S, S))
  for (i in 1:S) for (j in 1:S) for (k in 1:S) {
    dm_ij = (m[i] - m[j])^2
    dm_jk = (m[j] - m[k])^2
    dm_ki = (m[k] - m[i])^2
    num = 8*(dm_jk*v[i] + dm_ki*v[j] + dm_ij*v[k]) + 2*w^2*(dm_ij + dm_jk + dm_ki)
    denom = 16*(v[i]*v[j] + v[j]*v[k] + v[k]*v[i]) + 8*w^2*(v[i] + v[j] + v[k]) + 3*w^4
    eps[i,j,k] <- 2*w^2*exp(-num/denom) / sqrt(sqrt(pi)*w*denom)
  }
  eps
}

jac <- function(data, model, w, Sinit = 40) {
  spList <- as.numeric(data$species) # Indices of surviving species
  S <- length(spList)
  rseed <- data$rseedEps[1]
  kappa <- data$kappa[1]
  z0 <- data$z0[1]
  W <- data$W[1]
  m <- data$m
  dm <- outer(m, m, FUN = `-`)
  ddm <- outer(dm, m, FUN = `-`)
  v <- data$sigma^2
  sv <- outer(v, v, FUN = `+`)
  ssv <- outer(sv, v, FUN = `+`)
  if (model %in% c("hier", "hierHOI")) {
    alpha <- sigmoid(dm / sqrt(2*sv + w^2))
  } else {
    alpha <- w*exp(-dm^2 / (2*sv + w^2)) / sqrt(2*sv + w^2)
  }
  if (model == "evoHOI") {
    eps <- kappa * epsilonEvoHOI(m, v, w)
    # eps <- w / (2*ssv + w^2) * exp(-ddm^2 / (2*ssv + w^2))
  } else if (model == "hierHOI") {
    eps <- kappa * sigmoid((z0 + ddm) / sqrt(2*ssv + W^2))
  } else if (model %in% c("intraHOI", "HOI")) {
    eps <- assignRandomHOI(model, rseed, Sinit)[spList, spList, spList]
  } else {
    eps <- array(0, c(S, S, S))
  }
  alphaHOI <- matrix(0, S, S)
  for (i in 1:S) alphaHOI[i,] <- as.numeric((eps[i,,] + t(eps[i,,])) %*% data$n)
  J <- (-data$n) * (alpha + alphaHOI)
  eigen(J, only.values = TRUE)$values %>% abs() %>% log() %>% mean() %>% exp()
}


read_rds("data.rds") %>%
  filter(n > 1e-6) %>%
  nest(data = !model & !indVar & !w & !replicate) %>%
  mutate(diversity = map_dbl(data, ~invSimpson(.x$n)),
         richness = map_int(data, ~length(unique(.x$species))),
         clustering = map_dbl(data, ~traitPval(.x$m, 1000)),
         robustness = pmap_dbl(list(data, model, w), jac)) %>%
  mutate(indVar = fct_relevel(indVar, "low", "medium", "high"),
         model = fct_relevel(model, "pairwise", "evoHOI", "hier", "hierHOI",
                             "intraHOI", "HOI"),
         w = paste0("omega == ", w)) %>%
  select(-data) %>%
  write_rds("data_processed.rds", compress = "xz")
