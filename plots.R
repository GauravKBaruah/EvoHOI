library(tidyverse)

plotData <- function(dat, ydata, ylabel, yscaleFun = scale_y_continuous) {
  ggplot(dat, aes(x = indVar, y = {{ydata}}, colour = indVar, fill = indVar)) +
    geom_boxplot(alpha = 0.2, outlier.shape = NA) +
    geom_jitter(shape = 1, position = position_jitter(0.2)) +
    yscaleFun(name = ylabel) +
    labs(x = "Individual variation", y = ylabel) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    facet_grid(w ~ model, labeller = label_parsed) +
    theme_bw() +
    theme(legend.position = "none", panel.grid = element_blank())
}

dat <- read_rds("data_processed.rds")

plotData(dat, richness, "Species richness")
# ggsave("../../figures/richness.pdf", width = 8.4, height = 4.8)
plotData(dat, diversity, "Species diversity")
# ggsave("../../figures/diversity.pdf", width = 8.4, height = 4.8)
plotData(dat, clustering, "Trait clustering")
# ggsave("../../figures/clustering.pdf", width = 8.4, height = 4.8)
plotData(dat, robustness, "Average community robustness", scale_y_log10)
# ggsave("../../figures/robustness.pdf", width = 8.4, height = 4.8)
