library(tidyverse)

clargs <- commandArgs(trailingOnly = TRUE)
if (length(clargs) > 0) outfile <- clargs[1] else outfile <- "./data.rds"

tibble(file = Sys.glob("dat_*.rds")) %>%
  transmute(dat = map(file, read_rds)) %>%
  unnest(dat) %>%
  select(-datFile) %>%
  write_rds(file = outfile, compress = "xz")
