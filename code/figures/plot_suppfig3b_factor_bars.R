library(ggplot2); library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
ISO <- ifelse(length(args) >= 1, args[1], "GLB")
COV_LABEL <- ifelse(length(args) >= 2, args[2], "20%")
cov_tag <- switch(COV_LABEL, "20%" = "cov020", "50%" = "cov050", "affordable" = "afford")

d <- read.csv(sprintf("results/f3cv_depth_%s_%s.csv", cov_tag, ISO))
d$seg <- sub("Pathway open, flow insufficient", "Mediated, not offset", d$seg, fixed = TRUE)
d$seg <- sub("Unreachable \\(no pathway\\)", "No quantified pathway", d$seg)


seg_offset <- sprintf("Offset by GLP-1RA treatment (%s coverage)", ifelse(COV_LABEL == "affordable", "maximum affordable", COV_LABEL))
SEGS <- c(seg_offset, "Mediated, not offset", "No quantified pathway")
PAL <- c("#00A087", "#A8DED4", "#E64B35"); names(PAL) <- SEGS

iso_label <- ifelse(ISO == "GLB", "global", ISO)

w <- d %>% mutate(seg = factor(seg, SEGS)) %>%
  group_by(factor) %>% mutate(tot = sum(med)) %>% ungroup() %>%
  mutate(factor = reorder(factor, tot, max))

p <- ggplot(w, aes(med/1000, factor, fill = seg)) +
  geom_col(width = .78, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = PAL, name = NULL,
                    guide = guide_legend(nrow = 3)) +
  labs(x = sprintf("Diet-attributable cardiovascular deaths (thousands), %s, 2023", iso_label),
       y = NULL,
       title = sprintf("GLP-1RA treatment coverage: %s", ifelse(COV_LABEL == "affordable", "maximum affordable", COV_LABEL))) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top", legend.text = element_text(size = 7.2),
        legend.key.size = unit(8.5, "pt"), legend.justification = "left",
        axis.text.y = element_text(size = 8, face = "bold"),
        plot.title = element_text(face = "bold", size = 11))

ggsave(sprintf("figures/figure3_all/foodbar_cv_%s_%s.png", cov_tag, ISO), p,
       width = 5.6, height = 4.4, dpi = 420)
ggsave(sprintf("figures/figure3_all/foodbar_cv_%s_%s.pdf", cov_tag, ISO), p,
       width = 5.6, height = 4.4, device = cairo_pdf)
cat("done:", cov_tag, ISO, "\n")
