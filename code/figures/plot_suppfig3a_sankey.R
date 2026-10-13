library(ggplot2); library(dplyr); library(ggalluvial)

args <- commandArgs(trailingOnly = TRUE)
ISO <- args[1]
COV_LABEL <- ifelse(length(args) >= 2, args[2], "20%")

d <- read.csv(sprintf("results/f3cv_flows_%s_%s.csv", switch(COV_LABEL, "20%" = "cov020", "50%" = "cov050", "affordable" = "afford"), ISO))
d$seg <- sub("Pathway open, flow insufficient", "Mediated, not offset", d$seg, fixed = TRUE)
d$seg <- sub("Unreachable \\(no pathway\\)", "No quantified pathway", d$seg)


seg_offset <- sprintf("Offset by GLP-1RA treatment (%s coverage)", ifelse(COV_LABEL == "affordable", "maximum affordable", COV_LABEL))
seg_beyond <- "Mediated, not offset"
seg_none   <- "No quantified pathway"

d$seg <- gsub("Offset by GLP-1 \\(.*coverage\\)", seg_offset, d$seg)

thresh <- 0.02 * sum(d$deaths)
d <- d %>%
  group_by(factor) %>% mutate(ft = sum(deaths)) %>% ungroup() %>%
  mutate(factor = ifelse(ft < thresh, "Other (minor)", factor),
         mediator = factor(mediator, c("SBP", "FPG", "LDL", "Direct")),
         outcome = dplyr::recode(outcome, "Ischaemic stroke" = "IS", "Type 2 diabetes" = "T2D"),
         outcome = factor(outcome, c("IHD", "IS", "ICH", "SAH", "T2D", "CKD")),
         seg = factor(seg, c(seg_offset, seg_beyond, seg_none))) %>%
  count(factor, mediator, outcome, seg, wt = deaths, name = "deaths") %>%
  mutate(factor = reorder(factor, deaths, sum))

PAL <- setNames(c("#00A087", "#A8DED4", "#E64B35"), c(seg_offset, seg_beyond, seg_none))

iso_label <- ifelse(ISO == "GLB", "global", ISO)

p <- ggplot(d, aes(axis1 = factor, axis2 = mediator, axis3 = outcome, y = deaths/1000)) +
  geom_alluvium(aes(fill = seg), width = 1/9, alpha = .85) +
  geom_stratum(width = 1/9, fill = "grey96", colour = "grey50", linewidth = .22) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 2.4, fontface = "bold") +
  scale_x_discrete(limits = c("Dietary risk factor", "Mediator", "Outcome"),
                   expand = c(.14, .04)) +
  scale_fill_manual(values = PAL, name = NULL, drop = FALSE,
                    guide = guide_legend(nrow = 1)) +
  labs(y = NULL, caption = NULL) +
  theme_minimal(base_size = 11) +
  theme(axis.text.y = element_blank(), panel.grid = element_blank(),
        axis.title.x = element_blank(), legend.position = "top",
        legend.text = element_text(size = 7.4), legend.key.size = unit(9, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        plot.caption = element_text(size = 5.0, colour = "grey40", hjust = 0))

cov_tag <- switch(COV_LABEL, "20%" = "cov020", "50%" = "cov050", "affordable" = "afford")
out <- sprintf("figures/figure3_all/sankey_cv_%s_%s", cov_tag, ISO)
ggsave(paste0(out, ".pdf"), p, width = 5.6, height = 5.4, device = cairo_pdf)
ggsave(paste0(out, ".png"), p, width = 5.6, height = 5.4, dpi = 420)
cat("done:", out, "\n")
