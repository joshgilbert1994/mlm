# Figures for Chapter 5: Longitudinal Data and Linear Growth.
# Written by Bill Young. Run from the root of the book repository; set DATA
# to wherever your copy of dep.csv lives. Figures are written to OUT.
suppressMessages({library(tidyverse); library(lme4); library(patchwork); library(ggeffects)})
DATA <- "data/dep.csv"
OUT  <- "chapters/images/05/"; dir.create(OUT, showWarnings = FALSE)
set.seed(2026)

# The lecture's 12 patients, read exactly as the lecture reads them
dep <- read_csv(DATA, show_col_types = FALSE) |>
  filter(id <= 12) |>
  mutate(treat = factor(treat))

INK<-"#1A1D23"; BLU<-"#2C5F8A"; RED<-"#B03A2E"; GRN<-"#2E7D5B"; MUT<-"#5C6672"; RUL<-"#DDE2E6"
ORANGE <- "#D98E04"   # own-slope dots in fig10, chosen to contrast with BLU
th <- theme_bw(base_size = 14) + theme(
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(color = RUL, linewidth = .3),
  plot.title = element_text(face = "bold", size = 14, color = INK),
  plot.subtitle = element_text(size = 11.5, color = MUT),
  axis.title = element_text(size = 11.5, color = INK),
  strip.background = element_rect(fill = "grey94", color = NA),
  strip.text = element_text(face = "bold", size = 12, color = INK),
  legend.position = "bottom", legend.title = element_text(size = 11),
  plot.background = element_rect(fill = "white", color = NA))
theme_set(th)
TRT <- c("0" = BLU, "1" = RED)
trt_lab <- c("0" = "nonspecific therapy (0)", "1" = "disorder-specific therapy (1)")

# The lecture's models
ols <- lm(mdd9 ~ wave, dep)
m0 <- lmer(mdd9 ~ (1 | id), dep)
m1 <- lmer(mdd9 ~ wave + (1 | id), dep)
m2 <- lmer(mdd9 ~ wave + (wave | id), dep)
m5 <- lmer(mdd9 ~ wave * treat + (wave | id), dep)
dep <- dep |> mutate(fitted_m0 = fitted(m0), fitted_m1 = fitted(m1),
                     fitted_m2 = fitted(m2), fitted_m5 = fitted(m5))

v0 <- as.data.frame(VarCorr(m0)); tau2 <- v0$vcov[1]; sig2 <- v0$vcov[2]
v2 <- as.data.frame(VarCorr(m2)); t1 <- v2$sdcor[2]; g1 <- fixef(m2)[2]

save <- function(p, f, w, h) { cat("###", f, "\n")
  ggsave(paste0(OUT, f), p, width = w, height = h, dpi = 200, bg = "white") }

## fig01 --- when each person was measured ------------------------------------
# Not in the lecture. Shows at a glance that the panel is unbalanced.
cnt <- dep |> count(id)
p1 <- ggplot(dep, aes(wave, factor(id, levels = 12:1), color = treat)) +
  geom_point(size = 2.2) +
  geom_text(data = cnt, aes(x = 29.3, y = factor(id, levels = 12:1), label = paste0("n = ", n)),
            inherit.aes = FALSE, hjust = 0, size = 3.7, color = INK) +
  scale_color_manual(values = TRT, labels = trt_lab) +
  scale_x_continuous(breaks = seq(0, 27, 3)) +
  coord_cartesian(xlim = c(0, 32), clip = "off") +
  labs(x = "wave", y = "patient", color = NULL,
       title = "Which sessions each patient was measured at",
       subtitle = "Each dot is one measurement. Gaps are sessions with no score; short rows end at an earlier last session.")
save(p1, "fig01_timeline.png", 8.6, 5.0)

## fig02 --- time means and person means (lecture pp. 12-13) ------------------
tm <- dep |> group_by(wave) |> summarise(m = mean(mdd9), n = n(), .groups = "drop")
pm <- dep |> group_by(id) |> summarise(m = mean(mdd9), n = n(), .groups = "drop")
p2a <- ggplot(tm, aes(wave, m, size = n)) +
  geom_point(shape = 21, stroke = .9, color = INK) +
  scale_size_area(max_size = 5.5, breaks = c(4, 8, 12)) +
  coord_cartesian(ylim = c(0, 7)) +
  labs(x = "wave", y = "mean MDD9 across patients", size = "patients",
       title = "Time means", subtitle = "One dot per session: average over whoever was measured")
p2b <- ggplot(pm, aes(factor(id), m, size = n)) +
  geom_point(shape = 21, stroke = .9, color = INK) +
  scale_size_area(max_size = 5.5, breaks = c(5, 15, 25)) +
  coord_cartesian(ylim = c(0, 7)) +
  labs(x = "patient", y = "mean MDD9 across sessions", size = "sessions",
       title = "Person means", subtitle = "One dot per patient: average over their own sessions")
save(p2a + p2b, "fig02_time_person_means.png", 10.2, 4.4)

## fig03 --- sources of variation (lecture p. 15) -------------------------------
# each patient's mean with an approximate 95% confidence interval: mean +/- 1.96 x SE, SE = own SD / sqrt(n)
pmean <- dep |> group_by(id) |>
  summarise(m = mean(mdd9), n = n(), se = sd(mdd9) / sqrt(n), .groups = "drop") |>
  mutate(lo = m - 1.96 * se, hi = m + 1.96 * se)
print(pmean |> mutate(across(c(m, se, lo, hi), ~ round(.x, 2))), n = 12)
p3 <- ggplot(dep, aes(factor(id), mdd9)) +
  geom_hline(yintercept = mean(dep$mdd9), color = MUT, linetype = "dashed", linewidth = .5) +
  geom_jitter(width = .12, height = .12, alpha = .45, size = 1.7, color = INK) +
  geom_errorbar(data = pmean, aes(factor(id), ymin = lo, ymax = hi), inherit.aes = FALSE,
                color = RED, width = .28, linewidth = .8) +
  geom_point(data = pmean, aes(factor(id), m), color = RED, size = 3.6, shape = 18) +
  scale_y_continuous(breaks = 0:9) + coord_cartesian(ylim = c(-0.3, 9.3)) +
  labs(x = "patient", y = "MDD9 score (0 to 9)",
       title = "Two kinds of variation in the same 232 scores",
       subtitle = "Red: each patient's mean, plus or minus 1.96 standard errors. Dashed line: the overall mean, 3.24.")
save(p3, "fig03_sources.png", 8.6, 4.6)

## fig04 / 05 / 07 / 08 --- the trellis, raw and with each model (lecture p. 16, 20, 29, 33)
trellis <- function(yvar = NULL, ttl, sub) {
  g <- ggplot(dep, aes(wave, mdd9, color = treat)) +
    geom_point(alpha = .55, size = 1.5) +
    facet_wrap(~ id, ncol = 4, labeller = labeller(id = function(x) paste("patient", x))) +
    scale_color_manual(values = TRT, labels = trt_lab) +
    scale_y_continuous(breaks = c(0, 3, 6, 9)) +
    coord_cartesian(ylim = c(-0.4, 9.3)) +
    scale_x_continuous(breaks = c(0, 10, 20)) +
    labs(x = "wave", y = "MDD9 score", color = NULL, title = ttl, subtitle = sub)
  if (!is.null(yvar)) g <- g + geom_line(aes(y = .data[[yvar]]), linewidth = .8)
  g
}
save(trellis(NULL, "Each patient's scores over the study",
             "Twelve patients, one panel each. Look for level, direction and how much each patient varies."),
     "fig04_trellis_raw.png", 9.0, 6.4)
save(trellis("fitted_m0", "Unconditional means model: one flat line per patient",
             "Each line is the patient's estimated mean, pulled part of the way toward the average patient's mean, 3.39."),
     "fig05_trellis_m0.png", 9.0, 6.4)
save(trellis("fitted_m1", "Linear growth model with random intercepts",
             "Every patient's line has the same slope. Only the intercept differs from patient to patient."),
     "fig07_trellis_m1.png", 9.0, 6.4)
save(trellis("fitted_m2", "Linear growth model with random slopes",
             "Each patient now has their own intercept and their own slope."),
     "fig08_trellis_m2.png", 9.0, 6.4)

## fig06 --- reliability of each patient's mean (lecture p. 24, labeled) -----
rel <- dep |> count(id) |> mutate(rel = tau2 / (tau2 + sig2 / n))
nn <- tibble(n = 1:30, rel = tau2 / (tau2 + sig2 / n))
p6 <- ggplot(rel, aes(n, rel)) +
  geom_line(data = nn, color = MUT, linewidth = .6, linetype = "dashed") +
  geom_point(size = 3, color = BLU) +
  geom_text(aes(label = paste0("patient ", id)), data = filter(rel, id %in% c(6, 8)),
            nudge_x = 1.2, hjust = 0, size = 3.8, color = INK) +
  geom_text(aes(label = paste0("patient ", id)), data = filter(rel, id == 4),
            nudge_y = -0.035, size = 3.8, color = INK) +
  scale_y_continuous(limits = c(0.55, 1), breaks = seq(0.6, 1, .1)) +
  scale_x_continuous(breaks = seq(0, 30, 5)) +
  labs(x = expression("number of sessions observed, " * n[j]),
       y = "reliability of the patient's mean",
       title = "More sessions, a more reliable average",
       subtitle = "Dots are our 12 patients; patients with the same number of sessions share a dot.\nThe dashed curve is the same formula at every possible number of sessions.")
save(p6, "fig06_reliability.png", 8.4, 4.8)

## fig09 --- the prediction interval on person slopes --------------------------
lo <- g1 - 1.96 * t1; hi <- g1 + 1.96 * t1
xs <- seq(g1 - 4 * t1, g1 + 4 * t1, length.out = 600)
dn <- tibble(x = xs, y = dnorm(xs, g1, t1))
ytop <- dnorm(g1, g1, t1)
p9 <- ggplot(dn, aes(x, y)) +
  geom_area(data = filter(dn, x >= lo, x <= hi), fill = "#D7E3EE") +
  geom_area(data = filter(dn, x >= 0, x <= hi), fill = "#F2D3CF") +
  geom_line(color = BLU, linewidth = 1.1) +
  geom_vline(xintercept = 0, color = INK, linewidth = .6) +
  geom_segment(x = g1, xend = g1, y = 0, yend = ytop, color = BLU, linetype = "dashed") +
  annotate("text", x = 0.004, y = ytop * 0.88, label = "no change", hjust = 0, size = 4, color = INK) +
  annotate("text", x = g1, y = ytop * 1.06, label = sprintf("average slope %.3f", g1), size = 4, color = BLU) +
  annotate("text", x = lo, y = dnorm(lo, g1, t1) + ytop * .12, label = sprintf("%.2f", lo), size = 4, color = INK) +
  annotate("text", x = hi + 0.004, y = dnorm(hi, g1, t1) + ytop * .05, hjust = 0, label = sprintf("%.2f", hi), size = 4, color = INK) +
  annotate("text", x = 0.105, y = ytop * 0.42, label = "slopes above zero:\nsymptoms rising", size = 3.7, color = RED, lineheight = .9) +
  annotate("segment", x = 0.072, xend = 0.035, y = ytop * 0.36, yend = ytop * 0.12, color = RED, linewidth = .4) +
  labs(x = "a patient's own slope (change in MDD9 per week)", y = NULL,
       title = "The range expected to contain about 95% of patients' slopes",
       subtitle = "Shaded: the interval from -0.23 to 0.06. Red: the part of that interval above zero.") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank())
save(p9, "fig09_slope_interval.png", 8.4, 4.4)

## fig10 --- each patient's own slope and their random slopes estimate -----------
own <- dep |> group_by(id) |>
  summarise(n = n(), own = coef(lm(mdd9 ~ wave))[2], .groups = "drop") |>
  mutate(model = coef(m2)$id[as.character(id), 2],
         lab = paste0("patient ", id, "  (n = ", n, ")"),
         lab = fct_reorder(lab, n))
XMAX <- 0.12   # patient 6's own slope (+2.50) is far off this scale; its segment runs off the right edge
own_plot <- own |> mutate(own_shown = pmin(own, 1))
p10 <- ggplot(own_plot) +
  geom_vline(xintercept = fixef(m2)[2], color = RED, linetype = "dashed") +
  geom_vline(xintercept = 0, color = MUT, linewidth = .4) +
  geom_segment(aes(x = own_shown, xend = model, y = lab, yend = lab),
               color = "#A7B0BA", linewidth = .9) +
  geom_point(aes(own_shown, lab), data = filter(own_plot, own <= XMAX), color = ORANGE, size = 3.0) +
  geom_point(aes(model, lab), color = BLU, size = 3.0) +
  annotate("label", x = XMAX, y = "patient 6  (n = 3)", label = "own slope +2.50, off the scale",
           hjust = 1, vjust = -0.35, size = 3.6, color = INK, fill = "white", border.colour = NA) +
  annotate("label", x = fixef(m2)[2], y = 12.75, label = "average -0.081",
           color = RED, size = 3.8, vjust = 0, fill = "white", border.colour = NA) +
  coord_cartesian(xlim = c(-0.52, XMAX), ylim = c(0.6, 13.0), clip = "on") +
  labs(x = "slope (change in MDD9 per week)", y = NULL,
       title = "Each patient's own slope and their random slopes estimate",
       subtitle = "Orange dot: a separate regression on that patient alone. Blue dot: the random slopes model.")
save(p10, "fig10_shrinkage.png", 8.8, 5.6)

## fig11 --- four ways to draw the same 12 patients -------------------------
grid <- expand_grid(id = 1:12, wave = 0:27) |>
  left_join(dep |> group_by(id) |> summarise(last = max(wave), .groups = "drop"), by = "id") |>
  filter(wave <= last) |>
  mutate(id = as.integer(id))
cf0 <- coef(m0)$id; cf1 <- coef(m1)$id; cf2 <- coef(m2)$id
ownc <- dep |> group_by(id) |> summarise(a = coef(lm(mdd9 ~ wave))[1], b = coef(lm(mdd9 ~ wave))[2], .groups = "drop")
lines4 <- bind_rows(
  grid |> mutate(model = "A. Unconditional means", y = cf0[as.character(id), 1]),
  grid |> mutate(model = "B. Random intercepts", y = cf1[as.character(id), 1] + cf1[as.character(id), 2] * wave),
  grid |> mutate(model = "C. Random slopes", y = cf2[as.character(id), 1] + cf2[as.character(id), 2] * wave),
  grid |> left_join(ownc, by = "id") |> mutate(model = "D. Separate regressions", y = a + b * wave)) |>
  left_join(dep |> distinct(id, treat), by = "id")
p11 <- ggplot(lines4, aes(wave, y, group = id, color = treat)) +
  geom_line(linewidth = .7, alpha = .85) +
  facet_wrap(~ model, nrow = 1) +
  scale_color_manual(values = TRT, labels = trt_lab) +
  geom_hline(yintercept = 0, color = MUT, linewidth = .3, linetype = "dotted") +
  coord_cartesian(ylim = c(-0.9, 9)) +
  scale_x_continuous(breaks = c(0, 10, 20)) +
  labs(x = "wave", y = "fitted MDD9", color = NULL,
       title = "The same 12 patients under four different models",
       subtitle = "Each line is one patient's fitted trajectory over the sessions they were observed.")
save(p11, "fig11_four_models.png", 11.0, 4.6)

## fig12 --- appendix: treatment by time (lecture pp. 50-51) -----------------
p12a <- ggplot(dep, aes(wave, fitted_m5, color = treat, group = id)) +
  geom_line(linewidth = .7) +
  scale_color_manual(values = TRT, labels = trt_lab) +
  guides(color = "none") +
  geom_hline(yintercept = 0, color = MUT, linewidth = .3, linetype = "dotted") +
  coord_cartesian(ylim = c(-1.6, 9)) +
  labs(x = "wave", y = "fitted MDD9", color = NULL,
       title = "Each patient's fitted line", subtitle = "Treatment-by-time model, 12 patients")
ge <- as.data.frame(ggpredict(m5, terms = c("wave [0:27]", "treat")))
p12b <- ggplot(ge, aes(x, predicted, color = group, fill = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .18, color = NA) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = TRT, labels = trt_lab) +
  scale_fill_manual(values = TRT, labels = trt_lab) +
  geom_hline(yintercept = 0, color = MUT, linewidth = .3, linetype = "dotted") +
  coord_cartesian(ylim = c(-1.6, 9)) +
  labs(x = "wave", y = "predicted MDD9", color = NULL, fill = NULL,
       title = "The average patient in each group", subtitle = "Shaded bands are 95% confidence intervals")
save((p12a + p12b) + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
     "fig12_treatment.png", 10.4, 4.8)

cat("\nAll figures written to", OUT, "\n")
