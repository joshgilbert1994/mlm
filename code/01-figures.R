# =====================================================================
#  Figures for Chapter 1: Introduction and Variance Components Models.
#  Everything is computed live from the data, so the figures and the
#  numbers in the text cannot drift apart. Written by Bill Young.
#
#  Run from the root of the book repository. Set DATA to wherever your
#  copy of jebs_sub.csv lives; figures are written to OUT.
# =====================================================================
suppressMessages({library(tidyverse); library(lme4); library(patchwork)})

DATA <- "data/jebs_sub.csv"
OUT  <- "chapters/images/01/"
dir.create(OUT, showWarnings = FALSE)
set.seed(1)  # fixes the jitter in figures 1 and 3

RED <- "red"; BLUE <- "#1F78B4"; GREEN <- "#1B9E3F"; INK <- "#1F2024"
theme_set(theme_bw(base_size = 15))
theme_update(panel.grid.minor = element_blank(),
             plot.title = element_text(face = "bold", size = 17),
             plot.subtitle = element_text(size = 13, colour = "grey25"))

# The data and the variance components model
jebs <- read_csv(DATA, show_col_types = FALSE) |>
  mutate(cluster_id = factor(cluster_id))

GM <- mean(jebs$read_comp)                       # raw grand mean, 47.92
m0 <- lmer(read_comp ~ 1 + (1 | cluster_id), jebs)
vc <- as.data.frame(VarCorr(m0))
TAU2 <- vc$vcov[1]; SIG2 <- vc$vcov[2]; G0 <- fixef(m0)[[1]]

sch <- jebs |>
  group_by(cluster_id) |>
  summarise(n = n(), mean = mean(read_comp), sd = sd(read_comp), .groups = "drop") |>
  mutate(se      = sd / sqrt(n),
         lb      = mean - 1.96 * se,
         ub      = mean + 1.96 * se,
         R       = TAU2 / (TAU2 + SIG2 / n),      # reliability / shrinkage factor
         raw     = mean - G0,                     # raw residual
         eb      = raw * R,                       # Empirical Bayes estimate
         se_eb   = sqrt(TAU2 * (1 - R)),
         eb_lb   = eb - 1.96 * se_eb,
         eb_ub   = eb + 1.96 * se_eb,
         discard = 100 * (1 - R))                 # share thrown away as noise

jebs <- jebs |> group_by(cluster_id) |> mutate(demeaned = read_comp - mean(read_comp)) |> ungroup()

# ---- FIG 1: both layers at once -----------------------------------------
p1 <- ggplot(jebs, aes(cluster_id, read_comp)) +
  geom_jitter(width = .12, height = 0, alpha = .30, size = 1.6, colour = INK) +
  geom_hline(yintercept = GM, linetype = "dashed", linewidth = .8) +
  geom_point(data = sch, aes(cluster_id, mean), colour = RED, size = 5.5) +
  labs(x = "School ID", y = "Reading Score")
ggsave(paste0(OUT, "fig01_two_layers.png"), p1, width = 7.7, height = 4.62, dpi = 200, bg = "white")

# ---- FIG 2: between-school variation only -------------------------------
p2 <- ggplot(sch, aes(cluster_id, mean)) +
  geom_hline(yintercept = GM, linetype = "dashed", linewidth = .8) +
  geom_point(colour = RED, size = 5) +
  labs(x = "School ID", y = "Mean Reading Score")
ggsave(paste0(OUT, "fig02_between.png"), p2, width = 7.7, height = 3.74, dpi = 200, bg = "white")

# ---- FIG 3: within-school variation only --------------------------------
p3 <- ggplot(jebs, aes(cluster_id, demeaned)) +
  geom_jitter(width = .12, height = 0, alpha = .30, size = 1.6, colour = INK) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .8) +
  labs(x = "School ID", y = "Demeaned Reading Score")
ggsave(paste0(OUT, "fig03_within.png"), p3, width = 7.7, height = 3.74, dpi = 200, bg = "white")

# ---- FIG 4: overdispersion -- the same means, now with intervals --------
p4 <- ggplot(sch, aes(cluster_id, mean)) +
  geom_hline(yintercept = GM, linetype = "dashed", linewidth = .8) +
  geom_linerange(aes(ymin = lb, ymax = ub), linewidth = .9, colour = INK) +
  geom_point(colour = RED, size = 5) +
  labs(x = "School ID", y = "Mean Reading Score")
ggsave(paste0(OUT, "fig04_overdispersion.png"), p4, width = 7.7, height = 3.96, dpi = 200, bg = "white")

# ---- FIG 5: the reduced form, drawn ------------------------------------
# Illustrative, not a particular school: one below-average school intercept
# and one student sitting above it.
B0J <- 38; YIJ <- 48
grid <- tibble(x = seq(14, 86, length.out = 900)) |>
  mutate(between = dnorm(x, G0,  sqrt(TAU2)),     # spread of school means
         within  = dnorm(x, B0J, sqrt(SIG2)))     # spread of students in one school
TOP <- max(grid$between)
p5 <- ggplot(grid, aes(x)) +
  geom_area(aes(y = within),  fill = BLUE,  alpha = .10) +
  geom_area(aes(y = between), fill = GREEN, alpha = .10) +
  geom_line(aes(y = within),  colour = BLUE,  linewidth = 1.4) +
  geom_line(aes(y = between), colour = GREEN, linewidth = 1.4) +
  geom_vline(xintercept = G0,  colour = GREEN, linewidth = .8, alpha = .7) +
  geom_vline(xintercept = B0J, colour = BLUE,  linewidth = 1.1, linetype = "dashed") +
  annotate("text", x = G0 + 1.4, y = TOP * 1.07, label = "gamma[0]", parse = TRUE,
           colour = GREEN, size = 8, hjust = 0) +
  annotate("text", x = B0J - 1.2, y = TOP * .62, label = "beta[0][j]", parse = TRUE,
           colour = BLUE, size = 8, hjust = 1) +
  # u_j : grand mean -> this school's intercept
  annotate("segment", x = G0, xend = B0J, y = TOP * .45, yend = TOP * .45,
           colour = GREEN, linewidth = 1.1,
           arrow = arrow(length = unit(.32, "cm"), type = "closed")) +
  annotate("label", x = (G0 + B0J) / 2, y = TOP * .45, label = "u[j]", parse = TRUE,
           colour = INK, size = 6.5, label.size = .4) +
  # e_ij : this school's intercept -> this student
  annotate("segment", x = B0J, xend = YIJ, y = TOP * .17, yend = TOP * .17,
           colour = RED, linewidth = 1.1,
           arrow = arrow(length = unit(.32, "cm"), type = "closed")) +
  annotate("label", x = (B0J + YIJ) / 2, y = TOP * .17, label = "e[ij]", parse = TRUE,
           colour = INK, size = 6.5, label.size = .4) +
  annotate("point", x = YIJ, y = 0, colour = RED, size = 5) +
  annotate("text", x = YIJ, y = -TOP * .07, label = "Y[ij]", parse = TRUE,
           colour = RED, size = 7, vjust = 1) +
  scale_x_continuous(breaks = seq(20, 80, 20)) +
  scale_y_continuous(expand = expansion(mult = c(.12, .10))) +
  labs(x = "Reading Score", y = NULL,
       title = "Structural Decomposition of a Multilevel Model",
       subtitle = "From the grand mean down to one student's score") +
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.line.x = element_line(colour = INK, linewidth = .7))
ggsave(paste0(OUT, "fig05_structural_decomp.png"), p5, width = 7.7, height = 4.84, dpi = 200, bg = "white")

# ---- FIG 6: caterpillar plot --------------------------------------------
ord <- sch |> arrange(eb) |> pull(cluster_id) |> as.character()
cat_d <- sch |> mutate(cluster_id = factor(as.character(cluster_id), levels = ord))
p6 <- ggplot(cat_d, aes(cluster_id, eb)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .8) +
  geom_linerange(aes(ymin = eb_lb, ymax = eb_ub), linewidth = .8, colour = INK) +
  geom_point(size = 3, colour = INK) +
  labs(x = "School ID", y = "Empirical Bayes Estimate")
ggsave(paste0(OUT, "fig06_caterpillar.png"), p6, width = 7.7, height = 3.96, dpi = 200, bg = "white")

# ---- FIG 7: EB against raw residual, sized by n -------------------------
lim <- c(-8.8, 11)
p7 <- ggplot(sch, aes(raw, eb)) +
  geom_abline(slope = 1, intercept = 0, linewidth = .7, colour = INK) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .8) +
  geom_point(aes(size = n), shape = 21, fill = NA, stroke = .9, colour = INK) +
  geom_text(aes(label = cluster_id), hjust = -.25, vjust = -.45, size = 5, colour = INK) +
  scale_size_continuous(range = c(1.2, 8), name = "Sample Size",
                        breaks = c(30, 60, 90)) +
  coord_cartesian(xlim = lim, ylim = lim) +
  labs(x = "Raw Residual", y = "Empirical Bayes Estimate") +
  theme(legend.position = "bottom")
ggsave(paste0(OUT, "fig07_shrinkage.png"), p7, width = 6.16, height = 6.16, dpi = 200, bg = "white")

# ---- FIG 8: what shrinkage does, two panels -----------------------------
# Panel (a) is ordered by the RAW residual; panel (b) follows the same order.
ord_raw <- sch |> arrange(raw) |> pull(cluster_id) |> as.character()
mv  <- sch |> mutate(cluster_id = factor(as.character(cluster_id), levels = ord_raw))
pts <- mv |> select(cluster_id, raw, eb) |>
  pivot_longer(c(raw, eb), names_to = "kind", values_to = "val") |>
  mutate(kind = factor(kind, levels = c("raw", "eb"),
                       labels = c("Raw residual (no shrinkage)", "Empirical Bayes (shrunk)")))
p8a <- ggplot(mv, aes(cluster_id)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = .8) +
  geom_segment(aes(xend = cluster_id, y = raw, yend = eb), colour = "grey35", linewidth = .9) +
  geom_point(data = pts, aes(y = val, fill = kind, colour = kind),
             shape = 21, size = 4.4, stroke = 1.1) +
  geom_text(aes(y = pmax(raw, eb) + .9, label = paste0("n=", n)), size = 4.1, colour = INK) +
  scale_fill_manual(values = c("white", BLUE), name = NULL) +
  scale_colour_manual(values = c(INK, BLUE), name = NULL) +
  labs(x = NULL, y = "Points above / below grand mean",
       subtitle = "(a) Where each school lands. Lines join each raw residual to\nits shrunken estimate, and every estimate sits closer to zero.") +
  theme(legend.position = "bottom")
p8b <- ggplot(mv, aes(cluster_id, discard)) +
  geom_col(fill = "grey75", colour = NA, width = .68) +
  geom_text(aes(label = sprintf("%.1f%%", discard)), vjust = -.5, size = 4.1, colour = INK) +
  scale_y_continuous(expand = expansion(mult = c(0, .12))) +
  labs(x = "School ID", y = expression("Percent discarded as noise  (1 - R"[j]*")"),
       subtitle = "(b) Every school gets shrunk. This is the share of each school's raw\nresidual thrown away as noise: small for the big schools, but never zero.")
p8 <- (p8a / p8b) +
  plot_annotation(title = "What shrinkage does",
                  theme = theme(plot.title = element_text(face = "bold", size = 19)))
ggsave(paste0(OUT, "fig08_shrink_movement.png"), p8, width = 8.14, height = 9.0, dpi = 200, bg = "white")

cat("saved figures\n")
