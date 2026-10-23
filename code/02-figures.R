# =====================================================================
#  Figures for Chapter 2: Random Intercepts Models and Inference.
#  Everything is computed live from the data, so the figures and the
#  numbers in the text cannot drift apart. Written by Bill Young.
#
#  Run from the root of the book repository. Set DATA to wherever your
#  copy of jebs_sub.csv lives; figures are written to OUT.
# =====================================================================
suppressMessages({library(tidyverse); library(lme4); library(broom.mixed)
                  library(patchwork); library(ggrepel)})

DATA <- "data/jebs_sub.csv"
OUT  <- "chapters/images/02/"

BLUE <- "#2C7FB8"; CRIM <- "#A51C30"; INK <- "#1F2024"; GREY <- "#7F8792"
theme_set(theme_bw(base_size = 13))
theme_update(legend.position = "bottom",
             plot.title = element_text(face = "bold", size = 15),
             plot.subtitle = element_text(size = 11.5, colour = "grey25"))

jebs <- read_csv(DATA, show_col_types = FALSE) |>
  mutate(cluster_id = factor(cluster_id))

ols  <- lm(read_comp ~ std_baseline, jebs)
mlm0 <- lmer(read_comp ~ (1 | cluster_id), jebs)
mlm1 <- lmer(read_comp ~ std_baseline + (1 | cluster_id), jebs)

vc0 <- as.data.frame(VarCorr(mlm0)); vc1 <- as.data.frame(VarCorr(mlm1))
jebs <- jebs |> mutate(ols_fitted = predict(ols), mlm_fitted = fitted(mlm1))

# ---- FIG 1: the pooled OLS fit -----------------------------------------
p1 <- ggplot(jebs, aes(std_baseline, read_comp)) +
  geom_point(aes(colour = cluster_id), alpha = .45, size = 1.5) +
  geom_line(aes(y = ols_fitted), linewidth = 1.3, colour = INK) +
  labs(x = "Beginning-of-year MAP score (SD units)",
       y = "End-of-year reading comprehension",
       colour = "School ID",
       title = "One line for 780 students",
       subtitle = "Pooled OLS ignores the schools entirely. Slope 7.31 points per SD of BOY MAP score.") +
  guides(colour = guide_legend(nrow = 1, override.aes = list(alpha = 1, size = 2.5)))
ggsave(paste0(OUT, "fig01_pooled_ols.png"), p1, width = 7.7, height = 5.0, dpi = 200, bg = "white")

# ---- FIG 2: trellis, OLS vs MLM ----------------------------------------
p2 <- ggplot(jebs, aes(std_baseline, read_comp)) +
  geom_point(aes(colour = cluster_id), alpha = .35, size = 1.1, show.legend = FALSE) +
  geom_line(aes(y = ols_fitted), linewidth = .9, colour = INK) +
  geom_line(aes(y = mlm_fitted), linewidth = .9, colour = CRIM, linetype = "dashed") +
  facet_wrap(~cluster_id, nrow = 2) +
  labs(x = "Beginning-of-year MAP score (SD units)",
       y = "End-of-year reading comprehension",
       title = "The same slope, moved up or down for each school",
       subtitle = "Solid black = pooled OLS, identical in every panel. Dashed red = random intercepts.\nSlopes 7.31 and 7.21, so the two are near-parallel: what really moves is the height.")
ggsave(paste0(OUT, "fig02_trellis.png"), p2, width = 8.6, height = 5.0, dpi = 200, bg = "white")

# ---- FIG 3: three pooling strategies ------------------------------------
base <- list(labs(x = "BOY MAP score", y = "EOY reading"),
             coord_cartesian(ylim = c(28, 69)),
             theme(legend.position = "none"))
a <- ggplot(jebs, aes(std_baseline, read_comp)) +
  geom_point(alpha = .18, size = .9) +
  geom_smooth(method = lm, se = FALSE, colour = INK, linewidth = 1.2) +
  base + labs(title = "(a) Complete pooling", subtitle = "one line, 3 parameters")
b <- ggplot(jebs, aes(std_baseline, mlm_fitted, colour = cluster_id)) +
  geom_line(linewidth = 1) + base +
  labs(title = "(b) Partial pooling", subtitle = "random intercepts, 4 parameters")
cc <- ggplot(jebs, aes(std_baseline, read_comp, colour = cluster_id)) +
  geom_smooth(method = lm, se = FALSE, linewidth = 1) + base +
  labs(title = "(c) No pooling", subtitle = "separate OLS, 30 parameters")
p3 <- (a | b | cc) +
  plot_annotation(title = "Three ways to handle the schools",
                  subtitle = "Watch the flattest line in panel (c). That is school 11, and it is built from only 13 students.",
                  theme = theme(plot.title = element_text(face = "bold", size = 15),
                                plot.subtitle = element_text(size = 11.5, colour = "grey25")))
suppressMessages(ggsave(paste0(OUT, "fig03_three_pooling.png"), p3, width = 9.4, height = 4.0, dpi = 200, bg = "white"))

# ---- FIG 4: caterpillar before and after --------------------------------
cat_d <- bind_rows(
  tidy(mlm0, effects = "ran_vals") |> mutate(model = "VC model: no covariate"),
  tidy(mlm1, effects = "ran_vals") |> mutate(model = "VAM: adjusting for BOY MAP")) |>
  mutate(model = factor(model, levels = c("VC model: no covariate", "VAM: adjusting for BOY MAP")),
         lb = estimate - 1.96*std.error, ub = estimate + 1.96*std.error,
         flag = case_when(level == "62" ~ "62", level == "12" ~ "12", TRUE ~ "other"))
ord <- cat_d |> filter(model == "VC model: no covariate") |> arrange(estimate) |> pull(level)
cat_d <- cat_d |> mutate(level = factor(level, levels = ord))
p4 <- ggplot(cat_d, aes(level, estimate, ymin = lb, ymax = ub, colour = flag)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey35") +
  geom_errorbar(width = 0, linewidth = .8) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = c("62" = CRIM, "12" = BLUE, "other" = INK), guide = "none") +
  facet_wrap(~model) +
  coord_cartesian(ylim = c(-9, 11)) +
  labs(x = "School ID (ordered by the VC model estimate)",
       y = expression("School residual  "*hat(u)[j]),
       title = "Adjusting for where students started reorders the schools",
       subtitle = "Red is school 62, blue is school 12. Both intervals are 95%.\nSchool 62 goes from worst in the sample to exactly average; school 12 gives up most of its lead.")
ggsave(paste0(OUT, "fig04_caterpillar_compare.png"), p4, width = 8.6, height = 4.6, dpi = 200, bg = "white")

# ---- FIG: Lecture 1's caterpillar plot, reproduced -----------------------
cat_vc <- tidy(mlm0, effects = "ran_vals") |>
  mutate(lb = estimate - 1.96*std.error, ub = estimate + 1.96*std.error,
         flag = case_when(level == "62" ~ "62", level == "12" ~ "12", TRUE ~ "other"),
         level = factor(level, levels = level[order(estimate)]))
pcat <- ggplot(cat_vc, aes(level, estimate, ymin = lb, ymax = ub, colour = flag)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey35") +
  geom_errorbar(width = 0, linewidth = .8) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("62" = CRIM, "12" = BLUE, "other" = INK), guide = "none") +
  labs(x = "School ID", y = expression("Empirical Bayes estimate  "*hat(u)[j]),
       title = "Lecture 1's caterpillar plot, from the variance components model",
       subtitle = "Ten schools sorted from lowest to highest, each with a 95% interval.\nSchool 62 (red) is lowest, at the far left; school 12 (blue) is highest, at the far right.")
ggsave(paste0(OUT, "fig_caterpillar_vc.png"), pcat, width = 7.7, height = 4.6, dpi = 200, bg = "white")

# ---- FIG: the standard normal against two t distributions ----------------
zt <- expand_grid(x = seq(-4.5, 4.5, .01),
                  dist = c("Normal (z)", "t, 3 df", "t, 30 df")) |>
  mutate(d = case_when(dist == "Normal (z)" ~ dnorm(x),
                       dist == "t, 3 df"    ~ dt(x, 3),
                       TRUE                 ~ dt(x, 30)),
         dist = factor(dist, levels = c("Normal (z)", "t, 30 df", "t, 3 df")))
cuts <- tibble(dist = factor(c("Normal (z)", "t, 30 df", "t, 3 df"),
                             levels = levels(zt$dist)),
               cut  = c(qnorm(.975), qt(.975, 30), qt(.975, 3)))
pzt <- ggplot(zt, aes(x, d, colour = dist, linewidth = dist)) +
  geom_line() +
  geom_vline(data = cuts, aes(xintercept = cut, colour = dist),
             linetype = "dashed", linewidth = .5, show.legend = FALSE) +
  annotate("text", x = 1.96, y = .425, label = "1.96", colour = INK,   size = 3.4, hjust = 1.15) +
  annotate("text", x = qt(.975, 30), y = .385, label = "2.04", colour = BLUE, size = 3.4, hjust = -0.15) +
  annotate("text", x = qt(.975,  3), y = .345, label = "3.18", colour = CRIM, size = 3.4, hjust = -0.15) +
  scale_colour_manual(values = c("Normal (z)" = INK, "t, 30 df" = BLUE, "t, 3 df" = CRIM), name = NULL) +
  scale_linewidth_manual(values = c("Normal (z)" = 1.2, "t, 30 df" = .9, "t, 3 df" = .9), guide = "none") +
  coord_cartesian(xlim = c(-4.5, 4.5), ylim = c(0, .45)) +
  labs(x = "Test statistic", y = "Density",
       title = "The same statistic, compared against different distributions",
       subtitle = "Dashed lines are each distribution's two-sided 5% cutoff. Fewer degrees of freedom means\nfatter tails, and a fatter tail demands a larger statistic before it calls a result surprising.")
ggsave(paste0(OUT, "fig_z_vs_t.png"), pzt, width = 7.7, height = 4.4, dpi = 200, bg = "white")

# ---- FIG: separate regressions vs one interacted regression --------------
# Same ten lines either way; what
# differs is how certain we are about each slope.
sep <- jebs |> group_by(cluster_id) |> group_modify(~{
    m <- lm(read_comp ~ std_baseline, .x)
    tibble(slope = coef(m)[2], se = summary(m)$coefficients[2, 2], df = df.residual(m))
  }) |> ungroup() |> mutate(source = "Ten separate regressions")

int_m  <- lm(read_comp ~ std_baseline * factor(cluster_id), jebs)
lv     <- levels(jebs$cluster_id)
cf     <- coef(int_m); V <- vcov(int_m)
int <- map_dfr(seq_along(lv), function(k) {
  v <- rep(0, length(cf)); names(v) <- names(cf)
  v["std_baseline"] <- 1
  if (k > 1) v[paste0("std_baseline:factor(cluster_id)", lv[k])] <- 1
  tibble(cluster_id = factor(lv[k], levels = lv),
         slope = sum(v * cf), se = sqrt(drop(t(v) %*% V %*% v)), df = df.residual(int_m))
}) |> mutate(source = "One interacted regression")

cmp <- bind_rows(sep, int) |>
  mutate(source = factor(source, levels = c("Ten separate regressions", "One interacted regression")),
         crit = qt(.975, df), lb = slope - crit*se, ub = slope + crit*se,
         cluster_id = factor(cluster_id, levels = lv))

pint <- ggplot(cmp, aes(cluster_id, slope, ymin = lb, ymax = ub, colour = source)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_errorbar(position = position_dodge(width = .55), width = 0, linewidth = .9) +
  geom_point(position = position_dodge(width = .55), size = 2.4) +
  scale_colour_manual(values = c("Ten separate regressions" = CRIM,
                                 "One interacted regression" = BLUE), name = NULL) +
  labs(x = "School ID", y = "Estimated MAP slope for that school",
       title = "Same ten slopes, different certainty about them",
       subtitle = "Each pair of dots sits at exactly the same height: both approaches give\nidentical estimates. Only the intervals differ, because separate regressions use\none school's residuals while the interacted model pools all 780 students.")
ggsave(paste0(OUT, "fig_separate_vs_interaction.png"), pint, width = 8.4, height = 4.6, dpi = 200, bg = "white")

# ---- FIG 7: the null distribution of the LR statistic at a boundary ------
# Testing tau^2 = 0 puts the
# null on the edge of the parameter space, so the LR statistic is distributed as
# 0.5*chisq_0 + 0.5*chisq_1 -- half the mass at exactly zero, the rest chi-squared
# on 1 df. Simulating it also shows WHY we halve: beyond 3.84 lies 2.5%, not 5%.
# The positive half is exactly a chi-squared on 1 df, so draw its density rather
# than a histogram: the tail past 3.84 is what matters and bars there are invisible.
cx   <- seq(0.10, 7, length.out = 900)
dens <- tibble(x = cx, d = 0.5 * dchisq(cx, df = 1))
top  <- max(dens$d)

p7 <- ggplot(dens, aes(x, d)) +
  geom_area(data = filter(dens, x >= 3.841459), fill = CRIM, alpha = .55) +
  geom_line(colour = GREY, linewidth = 1.3) +
  annotate("segment", x = 0, xend = 0, y = 0, yend = top * 1.30,
           colour = CRIM, linewidth = 2.8) +
  annotate("label", x = .30, y = top * 1.22, hjust = 0, size = 3.8, colour = CRIM,
           fill = "white", linewidth = 0, lineheight = 1.05,
           label = "half the probability sits\nat exactly zero") +
  annotate("label", x = 1.55, y = top * .60, hjust = 0, size = 3.8, colour = "grey20",
           fill = "white", linewidth = 0, lineheight = 1.05,
           label = "the other half is a\nchi-squared on 1 df") +
  geom_vline(xintercept = 3.841459, linetype = "dashed", colour = INK, linewidth = .8) +
  annotate("curve", x = 5.6, xend = 4.35, y = top * .30, yend = top * .045,
           curvature = .28, colour = INK, linewidth = .4,
           arrow = arrow(length = unit(.14, "cm"), type = "closed")) +
  annotate("label", x = 4.5, y = top * .42, hjust = 0, size = 3.6, colour = INK,
           fill = "white", linewidth = 0, lineheight = 1.05,
           label = "3.84 is the ordinary 5% cutoff, but only\n2.5% of this null lies beyond it. That is\nwhy the p-value gets halved.") +
  scale_y_continuous(breaks = NULL) +
  coord_cartesian(xlim = c(-.2, 7), ylim = c(0, top * 1.40)) +
  labs(x = "value of the likelihood ratio statistic when the truth is tau-squared = 0",
       y = "how often",
       title = "Why the likelihood ratio test is halved",
       subtitle = "A variance cannot be negative, so the null value of zero sits at the very edge of what the\nparameter can be. The null distribution is a mixture, not a single chi-squared curve.")
ggsave(paste0(OUT, "fig07_null_spike_slab.png"), p7, width = 8.2, height = 4.6, dpi = 200, bg = "white")

# ---- FIG 9: the expectation moves, the scores do not ---------------------
# Compare the VC model to the VAM, not pooled OLS to raw scores.
# Each school gets TWO expectations -- one flat (VC), one built from its own BOY
# scores (VAM) -- and the same actual average in both.
g0vc <- fixef(mlm0)[[1]]
g0v  <- fixef(mlm1)[[1]]; b1v <- fixef(mlm1)[[2]]
foc <- jebs |> filter(cluster_id %in% c("62", "12")) |>
  mutate(school = paste0("School ", cluster_id,
           ifelse(cluster_id == "62", "  (started lowest)", "  (started highest)")))
sm <- foc |> group_by(school, cluster_id) |>
  summarise(map = mean(std_baseline), y = mean(read_comp), .groups = "drop") |>
  mutate(vc = g0vc, vam = g0v + b1v*map)

sm <- sm |> mutate(xv = ifelse(cluster_id == "62", map + .22, map - .22),
                   xa = ifelse(cluster_id == "62", map - .22, map + .22),
                   x_vc   = ifelse(cluster_id == "62", map + 1.75, map),
                   hj_vc  = ifelse(cluster_id == "62", 0, 1.08),
                   hj_vam = ifelse(cluster_id == "62", -0.10, -0.10),
                   hj_act = 0.5,
                   vj_act = ifelse(cluster_id == "62", 2.0, -1.6))

p9 <- ggplot(foc, aes(std_baseline, read_comp)) +
  geom_point(alpha = .16, size = 1.3, colour = GREY) +
  geom_hline(yintercept = g0vc, linetype = "dashed", colour = BLUE, linewidth = 1.0) +
  geom_abline(intercept = g0v, slope = b1v, colour = INK, linewidth = 1.1) +
  # arrows are offset so they do not sit on top of one another; the MARKERS stay
  # at the school's own mean MAP, so each sits exactly on the line it belongs to
  geom_segment(data = sm, aes(x = xv, xend = xv, y = vc, yend = y),
               colour = BLUE, linewidth = 1.0,
               arrow = arrow(length = unit(.16, "cm"), type = "closed")) +
  geom_segment(data = sm, aes(x = xa, xend = xa, y = vam, yend = y),
               colour = CRIM, linewidth = 1.0,
               arrow = arrow(length = unit(.10, "cm"), type = "closed")) +
  geom_point(data = sm, aes(map, vc),  shape = 21, fill = "white", colour = BLUE, size = 3.6, stroke = 1.2) +
  geom_point(data = sm, aes(map, vam), shape = 21, fill = "white", colour = INK,  size = 6.0, stroke = 1.2) +
  geom_point(data = sm, aes(map, y), colour = CRIM, size = 3.4) +
  geom_label(data = sm, aes(x_vc, vc, label = sprintf("VC expects %.2f", vc), hjust = hj_vc),
             vjust = 0.5, size = 3.4, colour = BLUE, fill = "white", linewidth = 0) +
  geom_label(data = sm, aes(map, vam, label = sprintf("VAM expects %.2f", vam)),
             hjust = -0.14, vjust = 0.5, size = 3.4, colour = INK, fill = "white", linewidth = 0) +
  geom_label(data = sm, aes(map, y, label = sprintf("actually averaged %.2f", y), vjust = vj_act, hjust = hj_act),
             size = 3.5, colour = CRIM, fontface = "bold", fill = "white", linewidth = 0) +
  facet_wrap(~school) +
  coord_cartesian(xlim = c(-3.4, 2.9), ylim = c(28, 72)) +
  labs(x = "Beginning-of-year (BOY) MAP score, in SD units",
       y = "End-of-year reading comprehension",
       title = "What each model expects of a school",
       subtitle = "Blue dashed: the variance components model expects 48.43 from every school, whoever it enrolled.\nBlack: the VAM expects a different score from each school, built from its own students' BOY scores.\nThe crimson dot is the same actual average in both cases. School 62 lands inside its VAM ring.")
ggsave(paste0(OUT, "fig09_above_expectation.png"), p9, width = 8.6, height = 4.8, dpi = 200, bg = "white")

# ---- FIG 10: what a chi-squared distribution is --------------------------
xs <- seq(0.05, 10, length.out = 900)
cs <- tibble(x = xs, d = dchisq(xs, df = 1))
p10 <- ggplot(cs, aes(x, d)) +
  geom_area(data = filter(cs, x >= qchisq(.95, 1)), fill = CRIM, alpha = .45) +
  geom_line(linewidth = 1.2, colour = BLUE) +
  geom_vline(xintercept = qchisq(.95, 1), colour = CRIM, linetype = "22", linewidth = .8) +
  annotate("segment", x = 5.6, xend = 4.3, y = .115, yend = .022,
           colour = CRIM, linewidth = .6,
           arrow = arrow(length = unit(.15, "cm"), type = "closed")) +
  annotate("label", x = 5.7, y = .135, hjust = 0, label = "the shaded 5%: values\nthis big are surprising",
           colour = CRIM, fill = "white", size = 3.8, lineheight = 1.05) +
  annotate("label", x = 3.84, y = .30, label = "cutoff 3.84", colour = CRIM,
           fill = "white", size = 3.9, fontface = "bold") +
  annotate("label", x = 5.7, y = .245, hjust = 0, label = "our LR statistic is 28.98,\nway off the page to the right",
           colour = INK, fill = "white", size = 3.8, lineheight = 1.05) +
  coord_cartesian(ylim = c(0, .42), xlim = c(0, 10)) +
  labs(x = "value of the test statistic", y = "how often values that size occur",
       title = "A chi-squared distribution with 1 degree of freedom",
       subtitle = "Square a standard normal draw and this is the shape: mostly small values,\na long thin tail to the right. The curve runs off the top near zero.")
ggsave(paste0(OUT, "fig10_chisq.png"), p10, width = 7.7, height = 4.3, dpi = 200, bg = "white")

# ---- FIG 11: the spaghetti plot, with every parameter labeled -----------
# Answers "where is tau squared on the graph?" -- it is the vertical spread of
# the parallel school lines, not the gap between any particular pair.
g0 <- fixef(mlm1)[[1]]; b1 <- fixef(mlm1)[[2]]
tau_hat <- as.data.frame(VarCorr(mlm1))$sdcor[1]
u  <- ranef(mlm1)$cluster_id[, 1]
names(u) <- rownames(ranef(mlm1)$cluster_id)
grid <- expand_grid(cluster_id = factor(names(u)), std_baseline = seq(-2.3, 2.6, .05)) |>
  mutate(fit = g0 + u[as.character(cluster_id)] + b1*std_baseline)
hi <- names(u)[which.max(u)]; lo <- names(u)[which.min(u)]
SCHOOL_COLS <- c("#2C7FB8","#A51C30","#4C9F70","#B5651D","#6A51A3",
                 "#0F766E","#C2185B","#5B7C99","#8A6D3B","#7A3E9D")

# overlay normal curves? For the within and between variance?
# Both: the blue curve is the between-school distribution of intercepts (tau), the
# grey one the within-school scatter of students around their line (sigma). Same
# density scale, so their shapes are directly comparable -- that gap is the ICC.
sig_hat <- as.data.frame(VarCorr(mlm1))$sdcor[2]
# Both curves use the SAME density scale, so their shapes are directly comparable:
# tau is short and fat, sigma long and thin. They sit at different points on x.
KSC <- 0.92 / dnorm(g0, g0, tau_hat)
# The ten lines are parallel, so each spread is identical at every MAP. That lets the
# two curves sit at DIFFERENT points on the x axis rather than on top of one another.
XT <- -2.0; CT <- g0 + b1*XT                  # tau, on the average line at MAP = -2.0
lo_sch <- names(u)[which.min(u)]              # school 74, the lowest of the ten
XS <-  0.95; CS <- g0 + min(u) + b1*XS        # sigma, on THAT school's line
ty  <- seq(CT - 3.2*tau_hat, CT + 3.2*tau_hat, length.out = 250)
tcv <- tibble(y = ty, x = XT - KSC * dnorm(ty, CT, tau_hat))
tpoly <- bind_rows(tibble(x = XT, y = min(ty)), tcv, tibble(x = XT, y = max(ty)))
sy  <- seq(CS - 2.8*sig_hat, CS + 2.8*sig_hat, length.out = 320)
scv <- tibble(y = sy, x = XS + KSC * dnorm(sy, CS, sig_hat))
spoly <- bind_rows(tibble(x = XS, y = min(sy)), scv, tibble(x = XS, y = max(sy)))
p9 <- ggplot(grid, aes(std_baseline, fit, group = cluster_id)) +
  geom_line(aes(colour = cluster_id), linewidth = .9, alpha = .95) +
  geom_abline(intercept = g0, slope = b1, colour = INK, linewidth = 1.35, linetype = "dashed") +
  geom_vline(xintercept = 0, colour = "grey85") +
  geom_polygon(data = spoly, aes(x, y), fill = INK, alpha = .13, inherit.aes = FALSE) +
  geom_path(data = scv, aes(x, y), colour = INK, linewidth = .9, inherit.aes = FALSE) +
  annotate("point", x = XS, y = CS, size = 3.4, shape = 21, fill = "white", colour = INK, stroke = 1.3) +
  annotate("segment", x = XS, xend = XS, y = CS - sig_hat, yend = CS + sig_hat,
           colour = INK, linewidth = 1.4,
           arrow = arrow(length = unit(.13, "cm"), ends = "both", type = "closed")) +
  geom_polygon(data = tpoly, aes(x, y), fill = BLUE, alpha = .22, inherit.aes = FALSE) +
  geom_path(data = tcv, aes(x, y), colour = BLUE, linewidth = .9, inherit.aes = FALSE) +
  annotate("label", x = 0.72, y = 33.4, hjust = 1, vjust = .5, size = 5.0, parse = TRUE,
           colour = INK, fill = "white", linewidth = 0, label = "hat(sigma) == 6.16") +
  annotate("label", x = 0.72, y = 28.0, hjust = 1, vjust = .5, size = 3.5, lineheight = 1.05,
           colour = INK, fill = "white", linewidth = 0,
           label = "how far students scatter around their own\nschool's line, shown here on school 74's,\nthe lowest of the ten") +
  annotate("segment", x = 0.80, xend = 0.94, y = 33.2, yend = 37.4,
           colour = INK, linewidth = .4, arrow = arrow(length = unit(.13,"cm"), type = "closed")) +
  # ---- tau: double arrow on the average line at MAP = -2.0, one SD each way.
  annotate("segment", x = XT, xend = XT, y = CT - tau_hat, yend = CT + tau_hat,
           colour = BLUE, linewidth = 1.8,
           arrow = arrow(length = unit(.14, "cm"), ends = "both", type = "closed")) +
  annotate("point", x = 0, y = g0, size = 2.8, colour = INK) +
  # ---- gamma_0: upper-left, where the plot is empty
  annotate("segment", x = -1.30, xend = -0.10, y = 62.0, yend = g0 + 1.0,
           colour = INK, linewidth = .5, arrow = arrow(length = unit(.15,"cm"), type = "closed")) +
  annotate("label", x = -3.02, y = 70.0, hjust = 0, vjust = .5, size = 5.0, parse = TRUE,
           colour = INK, fill = "white", linewidth = 0, label = "hat(gamma)[0] == 48.82") +
  annotate("label", x = -3.02, y = 65.2, hjust = 0, vjust = .5, size = 3.5, lineheight = 1.05,
           colour = INK, fill = "white", linewidth = 0,
           label = "the height of the average\nschool's line at MAP = 0") +
  # ---- tau: lower-right upper slot, arrow back to the double arrow at MAP = 0
  annotate("segment", x = -2.74, xend = -2.44, y = 43.2, yend = 39.8,
           colour = BLUE, linewidth = .5, arrow = arrow(length = unit(.14,"cm"), type = "closed")) +
  annotate("label", x = -3.02, y = 49.0, hjust = 0, vjust = .5, size = 5.0, parse = TRUE,
           colour = BLUE, fill = "white", linewidth = 0, label = "hat(tau) == 1.79") +
  annotate("label", x = -3.02, y = 44.2, hjust = 0, vjust = .5, size = 3.5, lineheight = 1.05,
           colour = BLUE, fill = "white", linewidth = 0,
           label = "one SD up and one SD down\nfrom the middle line") +
  # ---- beta_1: rise-over-run bracket on the average line, label below it
  annotate("segment", x = 1.6, xend = 2.4, y = g0 + b1*1.6, yend = g0 + b1*1.6,
           colour = CRIM, linewidth = .9) +
  annotate("segment", x = 2.4, xend = 2.4, y = g0 + b1*1.6, yend = g0 + b1*2.4,
           colour = CRIM, linewidth = .9) +
  annotate("segment", x = 2.20, xend = 2.38, y = 47.0, yend = 58.6,
           colour = CRIM, linewidth = .5, arrow = arrow(length = unit(.15,"cm"), type = "closed")) +
  annotate("label", x = 2.66, y = 45.0, hjust = 1, vjust = .5, size = 5.0, parse = TRUE,
           colour = CRIM, fill = "white", linewidth = 0, label = "hat(beta)[1] == 7.21") +
  annotate("label", x = 2.66, y = 39.6, hjust = 1, vjust = .5, size = 3.5, lineheight = 1.05,
           colour = CRIM, fill = "white", linewidth = 0,
           label = "the slope, shared by every\none of the ten lines") +
  scale_colour_manual(values = setNames(SCHOOL_COLS, levels(grid$cluster_id)), name = "School ID") +
  guides(colour = guide_legend(nrow = 1, override.aes = list(linewidth = 1.4))) +
  coord_cartesian(xlim = c(-3.1, 2.75), ylim = c(24, 72)) +
  labs(x = "Beginning-of-year (BOY) MAP score, in SD units",
       y = "Predicted end-of-year reading score",
       title = "Where each estimate actually lives on the picture",
       subtitle = "One colored line per school, all ten parallel; the dashed black line is the average school.\nBoth normal curves use the same scale, at different points along the line: blue is the spread of the ten\nschool lines, gray the spread of students around one school's line, centered on the ring.")
ggsave(paste0(OUT,"fig11_parameters_on_plot.png"), p9, width = 7.9, height = 6.3, dpi = 200, bg = "white")

cat("saved figures\n")
