library(tidyverse)


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
    facet_wrap(~ time, ncol = 1, scales = "free_y") +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "none",
          plot.margin = unit(c(0.6, 0.3, 0.3, 0.3), "cm"))
}

plot_snaps <- function(data, model, moments = c(0, Inf), res = 501) {
  if (model %in% c("hier", "hierHOI")) {
    bfun <- function(x) 1 - exp(-x)
    limits <- c(-0.1, NA, NA)
  } else {
    bfun <- function(x) (x > -0.5) * (x < 0.5)
    limits <- c(-0.6, 0.6, NA)
  }
  plot_snapshot(data, bfun, limits, res)
}


read_rds("data_fig1.rds") %>%
  mutate(n = case_when(
    time == 0 & model == "hier" ~ n / 3,
    time == 0 & model == "hierHOI" ~ n / 5,
    time == 0 & model == "intraHOI" ~ n / 3,
    TRUE ~ n
  )) %>%
  mutate(time = as_factor(ifelse(time == 0, "initial state", "final state"))) %>%
  nest(data = !model) %>%
  slice(c(6, 2, 5, 1, 3, 4)) %>%
  mutate(plot = map2(data, model, plot_snaps)) %>%
  pull(plot) %>%
  cowplot::plot_grid(plotlist = ., nrow = 2, align = "hv", vjust = 1.0,
                     labels = c("pairwise","hier","intraHOI","evoHOI","hierHOI","HOI"))
# ggsave("../../figures/examples.pdf", width = 8, height = 8)
