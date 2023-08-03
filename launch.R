library(tidyverse)


randSigmas <- function(S, indVar) {
  if (indVar == "low") {
    runif(S, 0.003, 0.009)
  } else if (indVar == "medium") {
    runif(S, 0.01, 0.03)
  } else {
    runif(S, 0.05, 0.1)
  }
}

generate_script <- function(jobname, projectID, time, memory, runstr) {
  str_c("#!/bin/bash\n",
        "#SBATCH -J ", jobname, "\n",
        "#SBATCH -A ", projectID, "\n",
        "#SBATCH -t ", time, "\n",
        "#SBATCH --mem=", memory, "\n",
        "#SBATCH -n 1\n",
        "# \n",
        "# Run single task in foreground\n",
        "module add R/4.0.3-nsc1-gcc-7.3.0\n",
        "cd ", Sys.getenv("PWD"), "\n",
        "Rscript ", runstr, "\n",
        "# \n",
        "# Script ends here")
}

bash_script <- function(jobname, projectID, time, memory, runstr, launch) {
  write(generate_script(jobname, projectID, time, memory, runstr), jobname)
  if (launch) {
    system(paste("sbatch", jobname))
    Sys.sleep(0.1)
  }
}


# Write temporary input data files
crossing(
  S = 40,
  theta = 0.5,
  replicate = 1:40,
  w = c(0.1, 0.2, 0.3),
  indVar = c("low", "medium", "high"),
  model = c("pairwise", "hier", "evoHOI", "hierHOI", "intraHOI", "HOI")
) %>%
  mutate(sigma = map2(S, indVar, randSigmas)) %>%
  mutate(h2 = map(S, ~runif(.x, 0.1, 0.15))) %>%
  mutate(kappa = if_else(model %in% c("evoHOI", "hierHOI"), 1, NA_real_)) %>%
  mutate(z0 = if_else(model %in% c("hierHOI"), 0.7, NA_real_)) %>%
  mutate(W = if_else(model %in% c("hier", "hierHOI"), 0.1, NA_real_)) %>%
  mutate(rseedEps = sample(100000:999999, nrow(.)),
         rseedEps = if_else(model %in% c("intraHOI", "HOI"), rseedEps, NA_integer_)) %>%
  mutate(pars = pmap(., list)) %>%
  mutate(ic = pmap(list(S, theta, model), function(S, theta, model) {
    if (model %in% c("hier", "hierHOI")) {
      c(rep(1, times = S), runif(S, min = 0, max = 1.5))
    } else {
      c(rep(1, times = S), runif(S, min = -theta, max = theta))
    }
  })) %>%
  select(pars, ic, model) %>%
  rowid_to_column("datFile") %>%
  mutate(datFile = sprintf("dat_%04d.rds", datFile)) %>%
  split(.$datFile) %>%
  walk(~write_rds(.x, .x$datFile))

# Submit jobs
tibble(datFile = Sys.glob("dat_*.rds")) %>%
  mutate(jobname = str_c("HOI", str_extract(datFile, "\\d+"), ".sh")) %>%
  mutate(projectID = "liu-2018-29", time = "01:30:00", memory = "16000") %>%
  mutate(runstr = str_c("simOnCluster.R ", datFile)) %>%
  mutate(launch = TRUE) %>%
  select(-datFile) %>%
  pwalk(., bash_script)
