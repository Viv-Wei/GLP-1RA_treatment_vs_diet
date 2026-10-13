suppressMessages({library(ggplot2); library(dplyr); library(tidyr); library(patchwork)})

d <- read.csv("tables/tableS_cost_analysis.csv", check.names = FALSE) %>%
  filter(Country != "GLB")
names(d)[names(d) == "GLP-1 cost/death (USD)"] <- "glp1"
names(d)[names(d) == "Diet cost/death (USD)"]  <- "diet"
stopifnot(all(c("Country", "Price scenario", "glp1", "diet") %in% names(d)), n_distinct(d$Country) == 19)

col_diet <- "#2d9b8a"
col_glp1 <- "#e85d5d"
ord <- d %>% filter(`Price scenario` == "Generic $140/yr") %>% arrange(glp1) %>% pull(Country)

d_long <- d %>%
  select(Country, `Price scenario`, `Dietary policy` = diet, `GLP-1RA treatment` = glp1) %>%
  pivot_longer(c(`Dietary policy`, `GLP-1RA treatment`), names_to = "Arm", values_to = "Expenditure") %>%
  mutate(Country = factor(Country, levels = ord),
         Arm = factor(Arm, c("Dietary policy", "GLP-1RA treatment")),
         Exp_k = Expenditure / 1000)

prices <- c("Branded" = "Branded",
            "Generic $140/yr" = "Generic, US$140 a year",
            "Generic $28/yr" = "Generic, US$28 a year")
BRK <- c(0, 5, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000)
lab_k <- function(b) ifelse(b >= 1000, paste0("$", b / 1000, "M"), paste0("$", b, "k"))

plots <- list()
for (i in seq_along(prices)) {
  sub <- d_long %>% filter(`Price scenario` == names(prices)[i])
  top <- max(sub$Exp_k) * 1.05; brk <- c()
  for (b in BRK[BRK > 0 & BRK <= top]) {
    last <- if (length(brk)) tail(brk, 1) else 0
    if (sqrt(b) - sqrt(last) >= 0.12 * sqrt(top)) brk <- c(brk, b)
  }
  plots[[i]] <- ggplot(sub, aes(x = Country, y = sqrt(Exp_k), fill = Arm)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    scale_fill_manual(values = c(col_diet, col_glp1)) +
    scale_y_continuous(breaks = sqrt(brk), labels = lab_k(brk)) +
    coord_polar(start = 0) +
    labs(title = prices[i], y = NULL, x = NULL, fill = NULL) +
    theme_minimal(base_size = 9) +
    theme(axis.text.x = element_text(size = 8, face = "bold"),
          axis.text.y = element_text(size = 7),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
}

p <- plots[[1]] + plots[[2]] + plots[[3]] +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

dir.create("figures", showWarnings = FALSE)
ggsave("figures/suppfig8_expenditure.pdf", p, width = 18, height = 7, device = cairo_pdf)
ggsave("figures/suppfig8_expenditure.png", p, width = 18, height = 7, dpi = 300)
cat("done: figures/suppfig8_expenditure.*\n")
print(d %>% group_by(`Price scenario`) %>%
        summarise(glp1_min = min(glp1), glp1_max = max(glp1), diet_min = min(diet), diet_max = max(diet)))
