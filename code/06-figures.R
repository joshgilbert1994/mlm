# Figures for Chapter 6: Non-Linear Growth, Autocorrelation, and Dynamic Models.
# Written by Bill Young. Run from the root of the book repository; set DATA
# to wherever your copy of dep.csv lives. Figures are written to OUT.
suppressMessages({library(tidyverse); library(lme4); library(patchwork); library(ggeffects)})
DATA <- "data/dep.csv"
OUT  <- "chapters/images/06/"; dir.create(OUT, showWarnings = FALSE)
set.seed(2026)

# The full data, read exactly as the lecture reads it
dep <- read_csv(DATA, show_col_types = FALSE) |>
  mutate(treat = factor(treat), wave_f = factor(wave), wave2 = wave^2,
         wc = (wave - 13.5) / 10, wc2 = wc^2) |>
  arrange(id, wave) |>
  group_by(id) |>
  mutate(mdd9_lag = lag(mdd9),
         c_mdd9_lag = mdd9_lag - mean(mdd9_lag, na.rm = TRUE)) |>
  ungroup()

INK<-"#1A1D23"; BLU<-"#2C5F8A"; RED<-"#B03A2E"; GRN<-"#2E7D5B"; MUT<-"#5C6672"; RUL<-"#DDE2E6"
ORANGE <- "#D98E04"
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
bq <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))

# Each canvas is scaled by (display width in the sheet) / (width it was designed for), so the
# printed text stays the same size after the display widths were reduced to fit 42 pages.
FAC <- c(fig01 = 80/90, fig02 = 72/80, fig03 = 86/94, fig04 = 80/90, fig05_ = 76/84, fig05b = 80/92,
         fig06 = 72/78, fig07 = 82/92, fig08 = 74/88, fig09b = 70/80, fig09_ = 72/78, fig10 = 74/80,
         fig12 = 78/88, fig13 = 76/84, fig14 = 80/90, fig15 = 80/90)
save <- function(p, f, w, h) {
  k <- names(FAC)[startsWith(f, names(FAC))]; s <- if (length(k)) FAC[[k[1]]] else 1
  cat("###", f, round(s, 3), "\n")
  ggsave(paste0(OUT, f), p, width = w * s, height = h * s, dpi = 200, bg = "white") }

# ---- the lecture's models (warnings are printed so nothing is hidden) ----------------
m_trt  <- lmer(mdd9 ~ 1 + treat * wave + (1 + wave | id), dep)           # lecture m3 / dep2 m1
m_lin  <- lmer(mdd9 ~ 1 + wave + (1 + wave | id), dep)                   # lecture m2
m_quad <- lmer(mdd9 ~ 1 + wave + wave2 + (1 + wave + wave2 | id), dep)   # lecture m4, fitted exactly as in the lecture (it prints a convergence warning)
m_mid  <- lmer(mdd9 ~ 1 + wave + wave2 + (1 + wave | id), dep)          # quadratic average only (the lecture's intermediate option)
m_ri   <- lmer(mdd9 ~ 1 + treat * wave + (1 | id), dep)
m_ar   <- lmer(mdd9 ~ 1 + treat * wave + (1 + wave | id) + ar1(0 + wave_f | id), dep)   # default settings, as in the lecture
dcomp  <- dep |> filter(!is.na(mdd9_lag))
m_dyn  <- lmer(mdd9 ~ 1 + treat * wave + c_mdd9_lag + (1 + wave | id), dcomp)
m_fe   <- lmer(mdd9 ~ 1 + wave_f + (1 + wave | id), dep)  # appendix m5
m_fet  <- lmer(mdd9 ~ 1 + wave_f * treat + (1 + wave | id), dep)         # appendix m6
m_cr1  <- lmer(mdd9 ~ 1 + wave + (1 + wave | id) + (1 | wave), dep)

dep <- dep |> mutate(f_trt = fitted(m_trt), f_lin = fitted(m_lin), f_quad = fitted(m_quad),
                     f_fe = fitted(m_fe))

# the lecture's quadratic fit, on the week scale
fq <- fixef(m_quad)
g0 <- fq[1]; g1 <- fq[2]; g2 <- fq[3]
cat("quadratic, week scale:", round(c(g0, g1, g2), 6), " turning point:", round(-g1 / (2 * g2), 2), "\n")
slope_at <- function(w) g1 + 2 * g2 * w
mean_at  <- function(w) g0 + g1 * w + g2 * w^2

## fig01 --- the treatment comparison (lecture PDF 7 and 9) --------------------------
p1a <- ggplot(dep, aes(wave, f_trt, color = treat, group = id)) +
  geom_line(alpha = .3, linewidth = .45) +
  scale_color_manual(values = TRT, labels = trt_lab) +
  guides(color = "none") +
  coord_cartesian(ylim = c(-1, 9.5)) +
  labs(x = "wave (weeks)", y = "predicted MDD9", title = "Each patient's predicted line",
       subtitle = "All 256 patients, treatment-by-time model")
ge <- as.data.frame(ggpredict(m_trt, terms = c("wave [0:27]", "treat")))
p1b <- ggplot(ge, aes(x, predicted, color = group, fill = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .18, color = NA) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = TRT, labels = trt_lab) +
  scale_fill_manual(values = TRT, labels = trt_lab) +
  coord_cartesian(ylim = c(1.5, 5.5)) +
  labs(x = "wave (weeks)", y = "predicted MDD9", color = NULL, fill = NULL,
       title = "The average patient in each group", subtitle = "Shaded bands are 95% confidence intervals")
save((p1a + p1b) + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
     "fig01_treatment.png", 10.4, 4.8)

## fig02 --- marginal trend: straight line and smooth curve (lecture PDF 12) -----------
wm <- dep |> group_by(wave) |> summarise(m = mean(mdd9), n = n(), .groups = "drop")
p2 <- ggplot(dep, aes(wave, mdd9)) +
  geom_smooth(method = lm, formula = y ~ x, color = NA, fill = "grey75") +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), color = NA, fill = "#F2D3CF") +
  geom_smooth(method = lm, formula = y ~ x, se = FALSE, color = INK) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE, color = RED) +
  geom_point(data = wm, aes(wave, m, size = n), shape = 21, color = INK, fill = "white", stroke = .8) +
  scale_size_area(max_size = 4, breaks = c(150, 200, 240)) +
  coord_cartesian(ylim = c(2.3, 5.2)) +
  labs(x = "wave (weeks)", y = "MDD9 score", size = "patients measured",
       title = "Average symptoms by week: a straight line and a smooth curve",
       subtitle = "Black: straight line. Red: smooth curve with no assumed shape. Circles: the mean at each week.")
save(p2, "fig02_marginal_trend.png", 8.6, 4.4)

## fig03 --- what the quadratic terms mean (not in the lecture) ------------------------
cur <- tibble(w = seq(0, 38, .1)) |> mutate(y = mean_at(w), inside = w <= 27)
tan <- tibble(w0 = c(0, 13.5, 27)) |> mutate(y0 = mean_at(w0), b = slope_at(w0),
                                             lab = ifelse(w0 == 0, sprintf("week %s: slope %.3f", w0, b), sprintf("week %s:\nslope %.3f", w0, b)))
tl <- tan |> rowwise() |> reframe(w0 = w0, w = pmax(c(w0 - 4.5, w0 + 4.5), 0), y = y0 + b * (pmax(c(w0 - 4.5, w0 + 4.5), 0) - w0))
tp <- -g1 / (2 * g2)
p3a <- ggplot() +
  geom_line(data = filter(cur, inside), aes(w, y), color = BLU, linewidth = 1.3) +
  geom_line(data = filter(cur, !inside), aes(w, y), color = BLU, linewidth = .9, linetype = "dashed") +
  geom_line(data = tl, aes(w, y, group = w0), color = RED, linewidth = .8) +
  geom_point(data = tan, aes(w0, y0), color = RED, size = 2.6) +
  geom_text(data = tan |> mutate(lx = c(1, 14.8, 16.2), ly = c(4.95, 3.86, 2.6), hj = c(0, 0, 0)),
            aes(lx, ly, label = lab, hjust = hj), size = 4.4, color = RED, lineheight = .9) +
  geom_vline(xintercept = 27, color = MUT, linetype = "dotted") +
  annotate("text", x = 27.4, y = 4.55, label = "last week\nof data", hjust = 0, size = 4.2, color = MUT, lineheight = .9) +
  annotate("point", x = tp, y = mean_at(tp), color = INK, size = 2.4) +
  annotate("text", x = tp - 1.5, y = mean_at(tp) - .2, label = sprintf("lowest point,\nweek %.0f", tp), size = 4.2, color = INK, hjust = .5, lineheight = .9) +
  coord_cartesian(xlim = c(-1, 39), ylim = c(2.4, 5.1)) +
  labs(x = "wave (weeks)", y = "predicted MDD9 (average patient)",
       title = "The average curve and its slope at three weeks",
       subtitle = "Each red segment has the slope of the curve at that week.")
sl <- tibble(w = seq(0, 38, .1)) |> mutate(s = slope_at(w), inside = w <= 27)
p3b <- ggplot(sl, aes(w, s)) +
  geom_hline(yintercept = 0, color = MUT, linewidth = .4) +
  geom_line(data = filter(sl, inside), color = RED, linewidth = 1.3) +
  geom_line(data = filter(sl, !inside), color = RED, linewidth = .9, linetype = "dashed") +
  geom_point(data = tan, aes(w0, b), color = RED, size = 2.6) +
  annotate("text", x = 1, y = -0.012, hjust = 0, size = 4.4, color = INK,
           label = sprintf("the slope increases by\n2 x %.5f = %.4f each week", g2, 2 * g2), lineheight = .95) +
  coord_cartesian(xlim = c(-1, 39), ylim = c(-0.12, 0.01)) +
  labs(x = "wave (weeks)", y = "slope (change in MDD9 per week)",
       title = "The slope itself changes with time",
       subtitle = "Slope at week w =\nlinear coefficient + 2 x quadratic coefficient x w")
save(p3a + p3b, "fig03_quadratic_anatomy.png", 11.0, 4.8)

## fig04 --- the quadratic model: spaghetti and individual curvature --------------------
p4a <- ggplot(dep, aes(wave, f_quad, group = id)) +
  geom_line(alpha = .25, linewidth = .45, color = INK) +
  coord_cartesian(ylim = c(-1, 9.5)) +
  labs(x = "wave (weeks)", y = "predicted MDD9", title = "Each patient's predicted curve",
       subtitle = "Quadratic growth model, all 256 patients")
curv <- coef(m_quad)$id |> rownames_to_column("id") |> mutate(q = wave2)
share_pos <- mean(curv$q > 0)
cat("share of patients with positive quadratic coefficient:", round(share_pos, 3), "\n")
p4b <- ggplot(curv, aes(q)) +
  geom_histogram(bins = 30, fill = "grey80", color = "white") +
  geom_vline(xintercept = 0, color = INK, linewidth = .6) +
  geom_vline(xintercept = g2, color = RED, linetype = "dashed", linewidth = .8) +
  annotate("text", x = max(curv$q), y = Inf, vjust = 1.4, hjust = 1, size = 3.7, color = RED,
           label = sprintf("dashed red line:\naverage %.5f", g2), lineheight = .9) +
  annotate("text", x = min(curv$q), y = Inf, vjust = 1.4, hjust = 0, size = 3.5, color = INK,
           label = sprintf("below 0: %.0f%%\nabove 0: %.0f%%", 100 * (1 - share_pos), 100 * share_pos), lineheight = .9) +
  labs(x = "patient's quadratic coefficient", y = "patients",
       title = "Each patient's quadratic coefficient",
       subtitle = "Above 0: the slope increases each week.\nBelow 0: the slope decreases each week.")
save(p4a + p4b, "fig04_quadratic_spaghetti.png", 10.6, 4.6)

## fig05 --- trellis: linear vs quadratic, patients 1-12 (lecture PDF 17) ---------------
d12 <- dep |> filter(id <= 12)
p5 <- ggplot(d12, aes(wave)) +
  geom_point(aes(y = mdd9), alpha = .35, size = 1.4, color = INK) +
  geom_line(aes(y = f_lin), color = RED, linetype = "dashed", linewidth = .8) +
  geom_line(aes(y = f_quad), color = INK, linewidth = .9) +
  facet_wrap(~ id, ncol = 4, labeller = labeller(id = function(x) paste("patient", x))) +
  scale_y_continuous(breaks = c(0, 3, 6, 9)) + coord_cartesian(ylim = c(-0.4, 9.3)) +
  scale_x_continuous(breaks = c(0, 10, 20)) +
  labs(x = "wave (weeks)", y = "MDD9 score",
       title = "Straight line or curve, patient by patient",
       subtitle = "Dashed red: linear growth model with random slopes. Black: quadratic growth model.")
save(p5, "fig05_trellis_quadratic.png", 9.0, 5.6)

## fig05b --- the intermediate option: average curvature only vs each patient's own curvature ---
set.seed(7); ids40 <- sample(unique(dep$id), 40)
dep <- dep |> mutate(f_mid = fitted(m_mid))
fmid <- fixef(m_mid); avgmid <- tibble(wave = seq(0, 27, .25)) |> mutate(y = fmid[1] + fmid[2] * wave + fmid[3] * wave^2)
avgq <- tibble(wave = seq(0, 27, .25)) |> mutate(y = mean_at(wave))
mk5 <- function(v, av, ttl, sub) ggplot(filter(dep, id %in% ids40), aes(wave, .data[[v]], group = id)) +
  geom_line(alpha = .45, linewidth = .5, color = INK) +
  geom_line(data = av, aes(wave, y), inherit.aes = FALSE, color = RED, linewidth = 1.3) +
  coord_cartesian(ylim = c(-1, 9.5)) +
  labs(x = "wave (weeks)", y = "predicted MDD9", title = ttl, subtitle = sub)
p5b <- mk5("f_mid", avgmid, "Intermediate model", "Every patient has the same curvature") +
  mk5("f_quad", avgq, "Quadratic growth model", "Each patient has their own curvature") +
  plot_annotation(subtitle = "The same 40 randomly chosen patients in both panels. Red: the average patient's curve.",
                  theme = theme(plot.subtitle = element_text(size = 11.5, color = MUT)))
save(p5b, "fig05b_intermediate.png", 10.6, 4.9)

## fig06 --- two departures for one measurement (not in the lecture) ---------------------
cl <- coef(m_lin)$id; fl <- fixef(m_lin)
p1d <- dep |> filter(id == 1) |> mutate(own = cl["1", 1] + cl["1", 2] * wave,
                                        avg = fl[1] + fl[2] * wave, e = mdd9 - own)
pick <- p1d |> filter(wave >= 8, wave <= 20) |> slice_max(e, n = 1, with_ties = FALSE)
cat("patient 1 picked week", pick$wave, ":", round(c(pick$mdd9, pick$avg, pick$own, pick$own - pick$avg, pick$e), 3), "\n")
xb <- pick$wave
p6 <- ggplot(p1d, aes(wave)) +
  geom_line(aes(y = avg), color = MUT, linetype = "dashed", linewidth = .9) +
  geom_line(aes(y = own), color = BLU, linewidth = 1.1) +
  geom_point(aes(y = mdd9), color = INK, size = 2.2, alpha = .6) +
  geom_point(data = pick, aes(y = mdd9), color = RED, size = 3.6) +
  annotate("segment", x = xb + .5, xend = xb + .5, y = pick$avg, yend = pick$own, color = BLU, linewidth = 1,
           arrow = arrow(length = unit(.12, "in"), ends = "both")) +
  annotate("segment", x = xb + .5, xend = xb + .5, y = pick$own, yend = pick$mdd9, color = RED, linewidth = 1,
           arrow = arrow(length = unit(.1, "in"), ends = "both")) +
  annotate("label", x = xb + 1, y = (pick$avg + pick$own) / 2, hjust = 0, size = 3.8, color = BLU, fill = "white", border.colour = NA,
           label = sprintf("patient's departure from the average line: %.2f", pick$own - pick$avg)) +
  annotate("label", x = xb + 1, y = pick$mdd9 + .5, hjust = 0, size = 3.8, color = RED, fill = "white", border.colour = NA,
           label = sprintf("residual e: %.2f (distance from patient 1's predicted line)", pick$e)) +
  annotate("text", x = 0.3, y = fl[1] + .4, label = "average patient's line", hjust = 0, size = 3.7, color = MUT) +
  annotate("text", x = 0.3, y = 5.55, label = "patient 1's predicted line", hjust = 0, size = 3.7, color = BLU) +
  scale_y_continuous(breaks = 0:9) + coord_cartesian(ylim = c(2.5, 9), xlim = c(0, 28)) +
  labs(x = "wave (weeks)", y = "MDD9 score",
       title = "One measurement: the patient's departure and the residual",
       subtitle = "Patient 1 under the linear growth model with random slopes. Red dot: one measurement.")
save(p6, "fig06_two_departures.png", 8.8, 4.3)

## fig07 --- what autocorrelation looks like (lecture PDF 22, simulated) -----------------
set.seed(43)
sim_ac <- function(phi) tibble(t = 1:100) |>
  mutate(e = as.numeric(arima.sim(model = if (phi == 0) list() else list(ar = phi), n = 100)),
         y = .2 * t + e, phi = phi)
sims <- bind_rows(sim_ac(-.9), sim_ac(0), sim_ac(.9)) |>
  mutate(lab = factor(sprintf("autocorrelation = %s", phi), levels = sprintf("autocorrelation = %s", c(-.9, 0, .9))))
p7 <- ggplot(sims, aes(t, y)) +
  geom_line(alpha = .45, color = INK) + geom_point(size = 1.1, color = INK) +
  geom_abline(intercept = 0, slope = .2, color = BLU, linewidth = .9) +   # the true common trend
  facet_wrap(~ lab, nrow = 1) +
  labs(x = "time", y = "y", title = "Simulated series with the same trend and three autocorrelations",
       subtitle = "Blue: the trend line. Left: points alternate sides of the line. Middle: no pattern. Right: long stretches on one side.")
save(p7, "fig07_autocorrelation_sim.png", 11.0, 4.2)

## fig08 --- the level-1 residuals and their lag (lecture PDF 24-25) ---------------------
dep <- dep |> mutate(res = mdd9 - f_trt) |> group_by(id) |> mutate(res_lag = lag(res)) |> ungroup()
rr <- cor(dep$res, dep$res_lag, use = "complete.obs"); cat("lag residual r:", round(rr, 4), "\n")
p8a <- ggplot(dep, aes(res)) + geom_histogram(bins = 30, fill = "grey80", color = "white") +
  labs(x = expression("residual (distance from the patient's predicted line), " * e[ij]), y = "measurements",
       title = "The residuals", subtitle = "Treatment-by-time model, 5,517 measurements")
p8b <- ggplot(filter(dep, !is.na(res_lag)), aes(res_lag, res)) +
  geom_point(alpha = .12, size = .9, color = INK) +
  geom_smooth(method = lm, formula = y ~ x, se = FALSE, color = RED, linewidth = 1) +
  annotate("label", x = -4.6, y = 5.6, hjust = 0, size = 4, label = sprintf("r = %.2f", rr), border.colour = NA) +
  coord_cartesian(xlim = c(-5, 6.3), ylim = c(-5, 6.3)) +
  labs(x = expression("residual at the previous measurement, " * e[i-1*","*j]),
       y = expression("residual at this measurement, " * e[ij]),
       title = "This residual and the previous one",
       subtitle = "Each dot: two consecutive residuals of one patient")
save(p8a + p8b, "fig08_residual_lag.png", 10.6, 4.8)

## fig09 --- correlation by distance: data and three models (not in the lecture) --------
wide <- dep |> select(id, wave, mdd9) |> pivot_wider(names_from = wave, values_from = mdd9) |>
  select(all_of(as.character(0:27)))   # columns in week order (pivot_wider orders them by first appearance)
C <- cor(wide, use = "pairwise.complete.obs")
obs <- map_dfr(1:15, ~ tibble(k = .x, r = mean(diag(C[-(1:.x), -((28 - .x + 1):28), drop = FALSE]))))
vc_ri <- as.data.frame(VarCorr(m_ri)); t0_ri <- vc_ri$vcov[1]; s2_ri <- vc_ri$vcov[2]
vc_rs <- as.data.frame(VarCorr(m_trt))
T0 <- vc_rs$vcov[1]; T1 <- vc_rs$vcov[2]; T01 <- vc_rs$vcov[3]; S2 <- vc_rs$vcov[4]
vca <- VarCorr(m_ar); A0 <- vca$id[1, 1]; A1 <- vca$id[2, 2]; A01 <- vca$id[1, 2]
ARv <- vca$id.1[1, 1]; phi <- attr(vca$id.1, "correlation")[1, 2]; RS2 <- sigma(m_ar)^2
cat("AR pieces: phi", round(phi, 4), " AR var", round(ARv, 4), " resid var", round(RS2, 4), "\n")
imp <- function(k, rs, a0, a1, a01, s2, arv = 0, ph = 0) {
  mean(sapply(0:(27 - k), function(t) {
    s <- t + k
    cv <- a0 + (t + s) * a01 + t * s * a1 + arv * ph^k
    v1 <- a0 + 2 * t * a01 + t^2 * a1 + arv + s2
    v2 <- a0 + 2 * s * a01 + s^2 * a1 + arv + s2
    cv / sqrt(v1 * v2) }))
}
cbd <- bind_rows(
  obs |> mutate(src = "observed data"),
  tibble(k = 1:15, r = t0_ri / (t0_ri + s2_ri), src = "random intercepts"),
  tibble(k = 1:15) |> mutate(r = map_dbl(k, imp, a0 = T0, a1 = T1, a01 = T01, s2 = S2), src = "random slopes"),
  tibble(k = 1:15) |> mutate(r = map_dbl(k, imp, a0 = A0, a1 = A1, a01 = A01, s2 = RS2, arv = ARv, ph = phi),
                             src = "random slopes + AR(1)"))
print(as.data.frame(cbd |> filter(k %in% c(1, 2, 3, 5, 10, 15)) |> mutate(r = round(r, 3)) |> pivot_wider(names_from = src, values_from = r)))
cbd$src <- factor(cbd$src, levels = c("observed data", "random intercepts", "random slopes", "random slopes + AR(1)"))
p9 <- ggplot(cbd, aes(k, r, color = src, linetype = src)) +
  geom_line(linewidth = 1.1) + geom_point(data = filter(cbd, src == "observed data"), size = 2.4) +
  scale_color_manual(values = c(INK, MUT, BLU, RED)) +
  scale_linetype_manual(values = c("solid", "dashed", "solid", "solid")) +
  scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13, 15)) +
  coord_cartesian(ylim = c(.5, .85)) +
  labs(x = "weeks apart", y = "correlation between the two scores", color = NULL, linetype = NULL,
       title = "How alike are two scores from the same patient?",
       subtitle = "Black: computed from the scores, with no model. Gray, blue and red: what each model implies.\nEach value is averaged over all pairs of weeks that far apart.") +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2))
save(p9, "fig09_corr_by_distance.png", 8.8, 4.6)

## fig10 --- four structures for the same average line (simulated, not real patients) ----
set.seed(611)
types <- c("Independent residuals (OLS)", "Random intercepts", "Random slopes", "Positive autocorrelation (AR(1))")
wk <- 0:27; pop <- 6 - 0.12 * wk
mk <- function(type, person) {
  if (type == types[1]) y <- pop + rnorm(28, 0, 1.0)
  if (type == types[2]) y <- pop + c(2.2, -2.0)[person] + rnorm(28, 0, .45)
  if (type == types[3]) y <- pop + c(.4, -.3)[person] + c(.13, -.06)[person] * wk + rnorm(28, 0, .45)
  if (type == types[4]) { repeat { e <- as.numeric(arima.sim(list(ar = .9), n = 28, sd = .6))
      if (sum(diff(sign(e)) != 0) %in% 3:6 && all(abs(e) < 2.6)) break }
    y <- pop + e }
  own <- pop
  if (type == types[2]) own <- pop + c(2.2, -2.0)[person]
  if (type == types[3]) own <- pop + c(.4, -.3)[person] + c(.13, -.06)[person] * wk
  tibble(wave = wk, y = y, own = own, person = paste("person", c("A", "B")[person]), type = type)
}

s4 <- map_dfr(types, \(tt) map_dfr(1:2, \(pp) mk(tt, pp))) |> mutate(type = factor(type, levels = types))
p10 <- ggplot(s4, aes(wave, y, color = person)) +
  geom_line(data = tibble(wave = wk, y = pop), aes(wave, y), inherit.aes = FALSE, color = MUT, linewidth = 1, linetype = "dashed") +
  geom_line(data = filter(s4, type %in% types[2:3]), aes(wave, own, color = person), linewidth = 1.1) +
  geom_line(alpha = .3, linewidth = .45) + geom_point(size = 2) +
  facet_wrap(~ type, ncol = 2) +
  scale_color_manual(values = c(BLU, ORANGE)) +
  coord_cartesian(ylim = c(0, 9)) + scale_y_continuous(breaks = c(0, 3, 6, 9)) +
  labs(x = "week", y = "outcome", color = NULL,
       title = "Four ways a person's scores can relate to the average line",
       subtitle = "Simulated data, two people per panel. Dashed gray: the average line.\nSolid colored lines: each person's predicted line where it differs from the average line.")
save(p10, "fig10_four_structures.png", 8.8, 5.8)

## fig11 --- a straight line that is too rigid (not in the lecture) ---------------------
dep <- dep |> mutate(r_lin = mdd9 - f_lin, r_quad = mdd9 - f_quad)
lagr <- dep |> group_by(id) |> filter(n() >= 20) |>
  summarise(n = n(), rl = cor(r_lin[-1], r_lin[-n()]), rq = cor(r_quad[-1], r_quad[-n()]), .groups = "drop") |>
  mutate(drop = rl - rq) |> arrange(desc(drop))
print(head(lagr, 6))
pooled <- dep |> group_by(id) |> mutate(rl1 = lag(r_lin), rq1 = lag(r_quad)) |> ungroup()
cat("pooled lag r, linear:", round(cor(pooled$r_lin, pooled$rl1, use = "complete.obs"), 3),
    " quadratic:", round(cor(pooled$r_quad, pooled$rq1, use = "complete.obs"), 3), "\n")
PID <- 191
pp <- dep |> filter(id == PID)
runs <- function(r) sum(diff(sign(r)) != 0)
mkp <- function(fv, rv, ttl) {
  x <- pp[[rv]]; rlag <- cor(x[-1], x[-length(x)])
  ggplot(pp, aes(wave)) +
    geom_segment(aes(xend = wave, y = .data[[fv]], yend = mdd9, color = .data[[rv]] > 0), linewidth = .7) +
    geom_line(aes(y = .data[[fv]]), color = INK, linewidth = 1) +
    geom_point(aes(y = mdd9, color = .data[[rv]] > 0), size = 2.4) +
    scale_color_manual(values = c("FALSE" = BLU, "TRUE" = RED), labels = c("below the predicted line or curve", "above the predicted line or curve")) +
    scale_y_continuous(breaks = 0:9) + coord_cartesian(ylim = c(0, 9)) +
    labs(x = "wave (weeks)", y = "MDD9 score", color = NULL, title = ttl,
         subtitle = sprintf("Correlation with the previous residual: %.2f\nChanges of side: %d", rlag, runs(x)))
}
p11 <- (mkp("f_lin", "r_lin", sprintf("Patient %d, linear growth model", PID)) +
        mkp("f_quad", "r_quad", sprintf("Patient %d, quadratic growth model", PID))) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
# fig11 is no longer used in the sheet (condensed Oct 9); the numbers above are still printed
# save(p11, "fig11_rigid_trend.png", 10.6, 4.8)

## fig12 --- why the lag is demeaned (not in the lecture) --------------------------------
re0 <- ranef(m_trt)$id |> rownames_to_column("id") |> transmute(id = as.numeric(id), u0 = `(Intercept)`)
lg <- dcomp |> left_join(re0, by = "id")
r_raw <- cor(lg$mdd9_lag, lg$u0); r_dm <- cor(lg$c_mdd9_lag, lg$u0)
cat("lag vs random intercept r, raw:", round(r_raw, 3), " demeaned:", round(r_dm, 3), "\n")
mkl <- function(xv, xl, ttl, r) ggplot(lg, aes(.data[[xv]], u0)) +
  geom_jitter(width = .18, height = 0, alpha = .07, size = .8, color = INK) +
  geom_smooth(method = lm, formula = y ~ x, se = FALSE, color = RED) +
  annotate("label", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, label = sprintf("r = %.2f", r), size = 4, border.colour = NA) +
  labs(x = xl, y = expression("patient's random intercept departure, " * hat(u)[0*j]), title = ttl)
p12 <- mkl("mdd9_lag", "previous score (raw)", "Raw previous score", r_raw) +
  mkl("c_mdd9_lag", "previous score minus the patient's average previous score", "Demeaned previous score", r_dm) +
  plot_annotation(subtitle = "Each dot is one measurement. The random intercepts come from the treatment-by-time model.",
                  theme = theme(plot.subtitle = element_text(size = 11.5, color = MUT)))
save(p12, "fig12_lag_vs_intercept.png", 10.6, 4.8)

## fig13 --- the dynamic model, patients 1-12 (lecture PDF 38) ---------------------------
dcomp <- dcomp |> mutate(f_dyn = fitted(m_dyn))
p13 <- ggplot(filter(dcomp, id <= 12), aes(wave, mdd9, color = treat)) +
  geom_point(alpha = .5, size = 1.4) + geom_line(aes(y = f_dyn), linewidth = .8) +
  facet_wrap(~ id, ncol = 4, scales = "free_y", labeller = labeller(id = function(x) paste("patient", x))) +
  scale_color_manual(values = TRT, labels = trt_lab) +
  scale_x_continuous(breaks = c(0, 10, 20)) +
  labs(x = "wave (weeks)", y = "MDD9 score", color = NULL,
       title = "The dynamic model, patient by patient",
       subtitle = "Each prediction: the patient's predicted line plus 0.30 times their demeaned previous score.")
save(p13, "fig13_dynamic_trellis.png", 9.0, 5.7)

## fig14 --- appendix: time fixed effects (lecture PDF 51-52) ----------------------------
ff <- fixef(m_fe); wmeans <- tibble(wave = 0:27, fe = ff[1] + c(0, ff[-1]))
avg3 <- wmeans |> mutate(lin = fixef(m_lin)[1] + fixef(m_lin)[2] * wave, quad = mean_at(wave)) |>
  pivot_longer(c(fe, lin, quad), names_to = "model", values_to = "y") |>
  mutate(model = factor(model, levels = c("fe", "quad", "lin"),
                        labels = c("time fixed effects (one mean per week)", "quadratic growth model", "linear growth model (straight line)")))
p14a <- ggplot(avg3, aes(wave, y, color = model, linetype = model)) +
  geom_line(linewidth = 1.1) + geom_point(data = filter(avg3, model == levels(avg3$model)[1]), size = 1.8) +
  scale_color_manual(values = c(INK, BLU, RED)) + scale_linetype_manual(values = c("solid", "solid", "dashed")) +
  labs(x = "wave (weeks)", y = "predicted MDD9 (average patient)", color = NULL, linetype = NULL,
       title = "Three versions of the average trajectory") +
  guides(color = guide_legend(ncol = 1), linetype = guide_legend(ncol = 1))
ids8 <- c(1, 3, 5, 9, 11, 12)
p14b <- ggplot(filter(dep, id %in% ids8), aes(wave, f_fe, group = id, color = factor(id))) +
  geom_line(linewidth = .9) +
  scale_color_brewer(palette = "Dark2") + guides(color = "none") +
  labs(x = "wave (weeks)", y = "predicted MDD9", title = "Six patients under time fixed effects",
       subtitle = "Each curve: the week means plus the patient's straight line")
save(p14a + p14b, "fig14_time_fixed_effects.png", 10.8, 5.0)

## fig15 --- appendix: treatment by week and the event study (lecture PDF 54-55) ---------
gt <- as.data.frame(ggpredict(m_fet, terms = c("wave_f", "treat"))) |> mutate(x = as.numeric(as.character(x)))
p15a <- ggplot(gt, aes(x + ifelse(group == "1", .2, -.2), predicted, color = group)) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0, alpha = .7) + geom_point(size = 1.8) +
  scale_color_manual(values = TRT, labels = trt_lab) +
  labs(x = "wave (weeks)", y = "predicted MDD9", color = NULL, title = "Each group's mean at each week",
       subtitle = "With 95% confidence intervals")
ev <- broom.mixed::tidy(m_fet, effects = "fixed", conf.int = TRUE) |> filter(str_detect(term, ":treat")) |>
  mutate(wave = parse_number(term))
print(ev |> filter(conf.high < 0 | conf.low > 0))
p15b <- ggplot(ev, aes(wave, estimate, ymin = conf.low, ymax = conf.high)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = MUT) +
  geom_errorbar(width = 0, color = INK) + geom_point(size = 1.8, color = INK) +
  labs(x = "wave (weeks)", y = "interaction (change in gap since wave 0)", title = "The event-study plot",
       subtitle = "Each interaction, with its 95% confidence interval")
save((p15a + p15b) + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
     "fig15_event_study.png", 10.8, 5.0)

## fig16 --- appendix: week as a crossed random effect (not in the lecture's plots) -------
cr <- tibble(wave = 0:27, line = fixef(m_cr1)[1] + fixef(m_cr1)[2] * wave,
             blup = ranef(m_cr1)$wave[as.character(0:27), 1]) |> mutate(pulled = line + blup) |>
  left_join(wmeans, by = "wave")
cat("SD of week departures:", round(as.data.frame(VarCorr(m_cr1)) |> filter(grp == "wave") |> pull(sdcor), 4), "\n")
p16 <- ggplot(cr, aes(wave)) +
  geom_line(aes(y = line), color = RED, linewidth = 1, linetype = "dashed") +
  geom_segment(aes(xend = wave, y = fe, yend = pulled), color = "#A7B0BA", linewidth = .7,
               arrow = arrow(length = unit(.05, "in"))) +
  geom_point(aes(y = fe), shape = 21, fill = "white", color = INK, size = 2.6, stroke = .8) +
  geom_point(aes(y = pulled), color = BLU, size = 2.4) +
  labs(x = "wave (weeks)", y = "predicted MDD9 (average patient)",
       title = "Week means pulled toward the straight line",
       subtitle = "Open circles: one mean per week (time fixed effects). Blue: week as a random effect. Dashed red: the straight line.")
# fig16 is no longer used in the sheet (condensed Oct 9)
# save(p16, "fig16_crossed_weeks.png", 8.8, 4.3)


## fig09b --- the two pieces of an AR(1) residual (simulated; not in the lecture) ---------
phi_hat <- attr(VarCorr(m_ar)$id.1, "correlation")[1, 2]
sd_ar <- attr(VarCorr(m_ar)$id.1, "stddev")[1]; sd_ind <- sigma(m_ar)
set.seed(91)
a <- numeric(28); a[1] <- rnorm(1, 0, sd_ar)
for (t in 2:28) a[t] <- phi_hat * a[t - 1] + rnorm(1, 0, sd_ar * sqrt(1 - phi_hat^2))
eps <- rnorm(28, 0, sd_ind)
tp9 <- tibble(week = rep(0:27, 3), value = c(a, eps, a + eps),
              piece = factor(rep(c("Autocorrelated piece, a", "Independent piece, ε", "Residual = a + ε"), each = 28),
                             levels = c("Autocorrelated piece, a", "Independent piece, ε", "Residual = a + ε")))
p9b <- ggplot(tp9, aes(week, value)) +
  geom_hline(yintercept = 0, color = MUT, linewidth = .5) +
  geom_line(color = "#A7B0BA", linewidth = .5) +
  geom_point(aes(color = value > 0), size = 2.2) +
  scale_color_manual(values = c("FALSE" = BLU, "TRUE" = RED), guide = "none") +
  facet_wrap(~ piece, ncol = 1) +
  labs(x = "week", y = "value (symptoms)",
       title = "The two pieces of a residual under the AR(1) model",
       subtitle = sprintf("Simulated values for one patient, using the model's estimates: φ = %.2f,\nSD of the autocorrelated piece %.2f, SD of the independent piece %.2f.", phi_hat, sd_ar, sd_ind))
save(p9b, "fig09b_two_pieces.png", 8.6, 5.2)

## fig17 --- the dynamic effect, one patient and all patients (not in the lecture) ---------
b2 <- fixef(m_dyn)[["c_mdd9_lag"]]
dcomp <- dcomp |> mutate(f_line = f_dyn - b2 * c_mdd9_lag,      # the patient's predicted line at that week
                         gap = mdd9 - f_line)                   # measurement minus predicted line
cat("slope of (score - predicted line) on the demeaned previous score:", round(coef(lm(gap ~ c_mdd9_lag, dcomp))[2], 4),
    " beta2:", round(b2, 4), "\n")
p10d <- filter(dcomp, id == 10)
p10m <- filter(dep, id == 10)
cat("patient 10, usual level:", round(mean(p10d$mdd9_lag), 3), "\n"); print(filter(p10d, wave %in% c(7, 9)) |>
  select(wave, mdd9_lag, c_mdd9_lag, f_line, f_dyn) |> mutate(across(everything(), ~ round(.x, 3))))
w7 <- filter(p10d, wave == 7); w9 <- filter(p10d, wave == 9)
p17a <- ggplot(p10d, aes(wave)) +
  geom_line(aes(y = f_line), linetype = "22", color = INK, linewidth = .7) +
  geom_segment(aes(xend = wave, y = f_line, yend = f_dyn, color = c_mdd9_lag > 0), linewidth = 1.1) +
  geom_point(data = p10m, aes(y = mdd9), color = "#A7B0BA", size = 2.4) +
  geom_line(aes(y = f_dyn), color = INK, linewidth = .45) +
  geom_point(aes(y = f_dyn), color = INK, size = 1.6) +
  scale_color_manual(values = c("FALSE" = BLU, "TRUE" = RED), guide = "none") +
  annotate("text", x = 8.6, y = 6.9, hjust = 0, vjust = 1, size = 3.6, lineheight = .95, color = INK,
           label = "week 7: previous score 6 (week 5), about 5\nabove the usual level (0.95), so the prediction\nis 0.30 × 5 = 1.5 above the predicted line") +
  annotate("segment", x = 8.5, xend = 7.35, y = 6.2, yend = 2.95, color = MUT, linewidth = .4,
           arrow = arrow(length = unit(.12, "cm"))) +
  annotate("text", x = 10.0, y = -0.45, hjust = 0, vjust = 1, size = 3.6, lineheight = .95, color = INK,
           label = "week 9: previous score 0, about 1 below\nthe usual level, so the prediction is\n0.30 × 1 = 0.3 below the predicted line") +
  annotate("segment", x = 9.9, xend = 9.2, y = -0.5, yend = 0.72, color = MUT, linewidth = .4,
           arrow = arrow(length = unit(.12, "cm"))) +
  scale_x_continuous(breaks = seq(0, 25, 5)) + scale_y_continuous(breaks = seq(0, 6, 2)) +
  coord_cartesian(ylim = c(-1.75, 6.9)) +
  labs(x = "wave (weeks)", y = "MDD9 score", title = "One patient (patient 10)",
       subtitle = "Gray dots: measurements. Dashed: predicted line. Solid black:\npredictions. Bars: 0.30 × (previous score − usual level).")
bins <- dcomp |> mutate(bin = round(c_mdd9_lag)) |> group_by(bin) |>
  summarise(n = n(), gap = mean(gap), .groups = "drop") |> filter(n >= 20)
print(bins |> mutate(line = round(b2 * bin, 2), gap = round(gap, 2)))
cat("measurements outside the plotted range of fig17 (right):",
    sum(dcomp$c_mdd9_lag < -5.5 | dcomp$c_mdd9_lag > 6.4 | abs(dcomp$gap) > 4.3), "of", nrow(dcomp), "\n")
p17b <- ggplot(dcomp, aes(c_mdd9_lag, gap)) +
  geom_hline(yintercept = 0, color = MUT, linewidth = .4) + geom_vline(xintercept = 0, color = MUT, linewidth = .4) +
  geom_jitter(width = .12, height = .12, alpha = .05, size = .7, color = INK) +
  geom_abline(intercept = 0, slope = b2, color = RED, linewidth = 1) +
  geom_point(data = bins, aes(bin, gap), shape = 21, fill = "white", color = INK, size = 2.8, stroke = 1) +
  annotate("label", x = 3.6, y = -3.3, label = sprintf("slope = %.2f", b2), color = RED, size = 4, border.colour = NA) +
  scale_x_continuous(breaks = -5:6) + scale_y_continuous(breaks = -4:4) +
  coord_cartesian(xlim = c(-5.5, 6.4), ylim = c(-4.3, 4.3)) +
  labs(x = "previous score minus the patient's usual level", y = "score minus the patient's predicted line",
       title = "All 5,261 measurements",
       subtitle = "Small dots: measurements. Circles: averages.\nRed line: slope 0.30, intercept 0.")
save(p17a + p17b + plot_layout(widths = c(1.08, 1)), "fig17_dynamic_effect.png", 10.6, 4.9)

## fig18 --- one dynamic effect, one per therapy group, one per patient (not in the lecture) --
m_dyt <- lmer(mdd9 ~ 1 + treat * (wave + c_mdd9_lag) + (1 + wave | id), dcomp)        # Section 5.7
m_rsl <- lmer(mdd9 ~ 1 + treat * wave + mdd9_lag + (1 + wave + mdd9_lag | id), dcomp) # Section 5.8, as in the lecture (raw previous score)
g20 <- fixef(m_dyt)[["c_mdd9_lag"]]; g21 <- fixef(m_dyt)[["treat1:c_mdd9_lag"]]
g20r <- fixef(m_rsl)[["mdd9_lag"]]; t2 <- attr(VarCorr(m_rsl)$id, "stddev")[["mdd9_lag"]]
lo <- g20r - 1.96 * t2; hi <- g20r + 1.96 * t2
pat <- coef(m_rsl)$id$mdd9_lag
cat("5.7:", round(c(g20, g21, g20 + g21), 4), " 5.8:", round(c(g20r, t2, lo, hi), 4),
    " patients' estimated effects: range", round(range(pat), 3), " SD", round(sd(pat), 4),
    " below/above the band:", sum(pat < lo), sum(pat > hi), "\n")
lab18 <- function(p, ttl, sub) p + geom_hline(yintercept = 0, color = MUT, linewidth = .4) +
  geom_vline(xintercept = 0, color = MUT, linewidth = .4) +
  scale_x_continuous(breaks = -3:3) + scale_y_continuous(breaks = -2:2) +
  coord_cartesian(xlim = c(-3, 3), ylim = c(-2.4, 2.45)) +
  labs(x = "previous score minus usual level", y = "difference in predicted score", title = ttl, subtitle = sub)
p18a <- lab18(ggplot() + geom_abline(intercept = 0, slope = b2, color = INK, linewidth = 1.1) +
                annotate("text", x = -2.9, y = 1.95, hjust = 0, label = sprintf("slope: %.2f", b2), size = 3.7),
              "Section 5.5", "One dynamic effect\nfor every patient")
p18b <- lab18(ggplot() + geom_abline(intercept = 0, slope = g20, color = BLU, linewidth = 1.1) +
                geom_abline(intercept = 0, slope = g20 + g21, color = RED, linewidth = 1.1, linetype = "42") +
                annotate("text", x = -2.9, y = 1.95, hjust = 0, size = 3.7, color = RED,
                         label = sprintf("disorder-specific: %.3f", g20 + g21)) +
                annotate("text", x = -2.9, y = 1.5, hjust = 0, size = 3.7, color = BLU,
                         label = sprintf("nonspecific: %.3f", g20)),
              "Section 5.7", "One dynamic effect\nper therapy group")
p18c <- lab18(ggplot() +
                annotate("polygon", x = c(-3.4, 3.4, 3.4, -3.4), y = c(-3.4 * hi, 3.4 * hi, 3.4 * lo, -3.4 * lo), fill = ORANGE, alpha = .13) +
                geom_abline(data = tibble(s = pat), aes(intercept = 0, slope = s), color = "#A7B0BA", alpha = .35, linewidth = .3) +
                geom_abline(intercept = 0, slope = c(lo, hi), color = ORANGE, linewidth = .6, linetype = "22") +
                geom_abline(intercept = 0, slope = g20r, color = INK, linewidth = 1.1) +
                annotate("text", x = 3.2, y = -2.15, hjust = 1, size = 3.7, color = "#9A6503",
                         label = sprintf("95%% range: %.2f to %.2f", lo, hi)) +
                annotate("text", x = -2.9, y = 1.95, hjust = 0, size = 3.7, label = sprintf("average: %.2f", g20r)),
              "Section 5.8", "One dynamic effect per patient\n(gray lines; raw previous score)")
p18 <- (p18a + p18b + p18c) + plot_annotation(
  title = "Each line's slope is a dynamic effect",
  theme = theme(plot.title = element_text(face = "bold", size = 14, color = INK)))
save(p18, "fig18_three_dynamic_models.png", 11.0, 4.4)

cat("done\n")
