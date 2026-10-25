# Figures for Chapter 7: Binary Outcomes and Multilevel Logistic Regression.
# Written by Bill Young. Run from the root of the book repository; set DATA
# to wherever your copy of thai.csv lives. Figures are written to OUT.
suppressMessages({library(tidyverse); library(lme4); library(broom.mixed); library(patchwork); library(ggeffects)})
DATA <- "data/thai.csv"
OUT  <- "chapters/images/07/"; dir.create(OUT, showWarnings = FALSE)
set.seed(2026)

# The data, read exactly as the lecture reads it
thai <- read_csv(DATA, show_col_types = FALSE) |>
  group_by(schoolid) |>
  mutate(prop_male = mean(male)) |>
  ungroup() |>
  mutate(male = factor(male))

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
SEX <- c("0" = BLU, "1" = ORANGE)
sex_lab <- c("0" = "girls (male = 0)", "1" = "boys (male = 1)")
pct <- function(x, d = 1) paste0(formatC(100 * x, format = "f", digits = d), "%")
FAC <- c(fig01 = 0.9268, fig02 = 0.9024, fig03 = 0.8837, fig04 = 0.8974, fig05 = 0.8666, fig06_ = 0.8605, fig06b = 0.9250, fig07 = 0.8666, fig08 = 0.9250, fig09 = 0.8666, fig10 = 0.8810, fig11_ = 0.8947, fig11b = 0.9250, fig12 = 0.8666, fig13 = 0.8666, fig14 = 0.8421)
save <- function(p, f, w, h) {
  k <- names(FAC)[startsWith(f, names(FAC))]; sc <- if (length(k)) FAC[[k[length(k)]]] else 1
  cat("###", f, round(sc, 3), "\n")
  ggsave(paste0(OUT, f), p, width = w * sc, height = h * sc, dpi = 200, bg = "white") }

# ---- the lecture's seven models ---------------------------------------------------------
m1 <- lm(retained ~ male + sch_ses, thai)
m2 <- lm(retained ~ male + sch_ses + pped, thai)
m3 <- lm(retained ~ male*pped + sch_ses, thai)
m4 <- glm(retained ~ male + sch_ses + pped, thai, family = binomial)
m5 <- glmer(retained ~ male + sch_ses + pped + (1|schoolid), thai, family = binomial)
m6 <- glmer(retained ~ male + sch_ses + pped + (male|schoolid), thai, family = binomial)
m7 <- glmer(retained ~ male + prop_male + sch_ses + pped + (1|schoolid), thai, family = binomial)

g0  <- fixef(m5)[["(Intercept)"]]; bm <- fixef(m5)[["male1"]]
tau <- as.data.frame(VarCorr(m5))$sdcor[1]
cat("m5: g0", g0, " male", bm, " tau", tau, "\n")

## fig01 --- a 0/1 outcome and its mean (lecture PDF 6-7) -------------------------------
cnt <- thai |> count(retained) |> mutate(lab = scales::comma(n))
p1a <- ggplot(cnt, aes(factor(retained), n)) +
  geom_col(fill = c(MUT, RED), width = .6) +
  geom_text(aes(label = lab), vjust = -.4, size = 4.2) +
  scale_x_discrete(labels = c("0" = "0: not retained", "1" = "1: retained")) +
  scale_y_continuous(limits = c(0, 7200), labels = scales::comma) +
  labs(x = NULL, y = "number of students", title = "The outcome, retained",
       subtitle = "Mean = 1,067 / 7,516 = 0.142: 14.2% were retained")
grp <- thai |> group_by(male, pped) |> summarise(p = mean(retained), n = n(), .groups = "drop") |>
  mutate(pp = factor(pped, labels = c("no preschool (pped = 0)", "preschool (pped = 1)")))
p1b <- ggplot(grp, aes(pp, p, fill = male)) +
  geom_col(position = position_dodge(.7), width = .62) +
  geom_text(aes(label = pct(p)), position = position_dodge(.7), vjust = -.4, size = 4) +
  scale_fill_manual(values = SEX, labels = sex_lab, name = NULL) +
  scale_y_continuous(limits = c(0, .25), labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "percent retained", title = "Percent retained in each group",
       subtitle = "Each bar is the mean of retained in that group")
save(p1a + p1b + plot_layout(widths = c(.8, 1.2)), "fig01_outcome.png", 10.6, 4.6)

## fig02 --- the LPM picture (lecture PDF 10), with observed proportions --------------------
lpm_lines <- expand_grid(male = factor(c(0, 1)), sch_ses = seq(-0.77, 1.49, length.out = 50))
lpm_lines$fit <- predict(m1, lpm_lines)
qs <- thai |> mutate(q = ntile(sch_ses, 5)) |> group_by(q, male) |>
  summarise(sch_ses = mean(sch_ses), obs = mean(retained), .groups = "drop")
set.seed(7)
p2 <- ggplot(thai, aes(sch_ses, retained, color = male)) +
  geom_jitter(height = .01, width = .05, alpha = .22, size = .9) +
  geom_line(data = lpm_lines, aes(y = fit), linewidth = 1.2) +
  geom_point(data = qs, aes(y = obs), shape = 21, fill = "white", size = 3.4, stroke = 1.3) +
  scale_color_manual(values = SEX, labels = sex_lab, name = NULL) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  labs(x = "school SES (sch_ses)", y = "retained (dots); probability (lines, circles)",
       title = "The linear probability model, m1",
       subtitle = "Lines: m1's predicted probability. Circles: observed proportion retained, by fifth of school SES")
save(p2, "fig02_lpm.png", 9.2, 5.0)

## fig03 --- what the LPM residuals look like (beyond the lecture) -------------------------
r2 <- tibble(fit = fitted(m2), res = resid(m2), y = factor(thai$retained))
p3a <- ggplot(r2, aes(fit, res, color = y)) +
  geom_jitter(width = 0, height = .01, alpha = .2, size = .8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = MUT) +
  scale_color_manual(values = c("0" = MUT, "1" = RED), labels = c("0" = "not retained (0)", "1" = "retained (1)"), name = NULL) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 3))) +
  labs(x = "predicted probability from m2", y = "residual (observed minus predicted)",
       title = "Residuals of the LPM, m2",
       subtitle = "Retained: 1 minus the prediction. Not retained: minus the prediction")
vv <- tibble(p = seq(0, 1, .005), v = p * (1 - p))
p3b <- ggplot(vv, aes(p, v)) +
  annotate("rect", xmin = min(r2$fit), xmax = max(r2$fit), ymin = -Inf, ymax = Inf, fill = BLU, alpha = .12) +
  geom_line(linewidth = 1.1, color = INK) +
  annotate("text", x = .134, y = .235, label = "m2's predictions:\n0.045 to 0.223", color = BLU, size = 3.6) +
  labs(x = "probability of a 1, p", y = "variance, p(1 - p)",
       title = "The variance depends on p",
       subtitle = "Largest at p = 0.5, smaller toward 0 and 1")
save(p3a + p3b + plot_layout(widths = c(1.25, 1)), "fig03_lpm_residuals.png", 10.6, 4.7)

## fig04 --- the logit and the S-curve (lecture PDF 18-19) --------------------------------
lt <- tibble(logit = seq(-4, 4, by = .05), prob = plogis(logit))
anc <- tibble(logit = -3:3, prob = plogis(logit), odds = exp(logit)) |>
  mutate(lab = sprintf("p = %.3f\nodds = %.3f", prob, odds),
         hj = ifelse(logit < 0, 1.1, -0.1), vj = ifelse(logit < 0, -0.15, 1.1))
p4 <- ggplot(lt, aes(logit, prob)) +
  geom_hline(yintercept = c(0, 1), color = MUT, linetype = "dotted") +
  geom_line(linewidth = 1.2, color = INK) +
  geom_point(data = anc, color = RED, size = 2.6) +
  geom_text(data = anc, aes(label = lab, hjust = hj, vjust = vj), size = 3.3, color = RED, lineheight = .9) +
  scale_x_continuous(breaks = -4:4) +
  scale_y_continuous(breaks = seq(0, 1, .25), limits = c(0, 1)) +
  labs(x = "log-odds (logit)", y = "probability",
       title = "Probability against log-odds",
       subtitle = "Every log-odds value corresponds to a probability between 0 and 1")
save(p4, "fig04_logit_curve.png", 8.6, 5.0)

## fig05 --- equal steps on one scale are not equal on the other (beyond the lecture) ----
lt2 <- tibble(logit = seq(-4.5, 4.5, by = .05), prob = plogis(logit))
steps <- function(p0, l1) tibble(p0 = p0, l0 = qlogis(p0), l1 = l1, p1 = plogis(l1))
seg <- steps(c(.50, .95), qlogis(c(.53, .98)))
p5a <- ggplot(lt2, aes(logit, prob)) +
  geom_line(linewidth = 1, color = MUT) +
  geom_segment(data = seg, aes(x = l0, xend = l1, y = p0, yend = p0), color = BLU, linewidth = 1.2,
               arrow = arrow(length = unit(.1, "cm"))) +
  geom_segment(data = seg, aes(x = l1, xend = l1, y = p0, yend = p1), color = ORANGE, linewidth = 1.2,
               arrow = arrow(length = unit(.1, "cm"))) +
  geom_point(data = seg, aes(l0, p0), color = BLU, size = 2) +
  annotate("text", x = 0.45, y = .40, label = "50% to 53%:\n0.12 in log-odds", hjust = 0, size = 3.6, color = INK, lineheight = .9) +
  annotate("text", x = 4.4, y = .84, label = "95% to 98%:\n0.95 in log-odds", hjust = 1, size = 3.6, color = INK, lineheight = .9) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "log-odds", y = "probability", title = "The same 3-point change in probability",
       subtitle = "is 0.12 in log-odds at 50% and 0.95 at 95%")
st <- steps(c(.05, .50, .95), qlogis(c(.05, .50, .95)) + 0.43) |>
  mutate(lab = sprintf("+%.1f points", 100 * (p1 - p0)))
p5b <- ggplot(lt2, aes(logit, prob)) +
  geom_line(linewidth = 1, color = MUT) +
  geom_segment(data = st, aes(x = l0, xend = l1, y = p0, yend = p0), color = BLU, linewidth = 1.2,
               arrow = arrow(length = unit(.1, "cm"))) +
  geom_segment(data = st, aes(x = l1, xend = l1, y = p0, yend = p1), color = ORANGE, linewidth = 1.2,
               arrow = arrow(length = unit(.1, "cm"))) +
  geom_point(data = st, aes(l0, p0), color = BLU, size = 2) +
  geom_text(data = st, aes(x = ifelse(p0 > .9, l0 - .3, l1 + .25), y = ifelse(p0 > .9, p0 + .02, (p0 + p1) / 2), label = lab),
            hjust = ifelse(st$p0 > .9, 1, 0), size = 3.6, color = ORANGE) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "log-odds", y = "probability", title = "The same 0.43 change in log-odds",
       subtitle = "is a different number of points at 5%, 50% and 95%")
save(p5a + p5b, "fig05_equal_steps.png", 10.8, 4.7)

## fig06 --- the logit model on the probability scale, with the LPM (lecture PDF 22) ------
ge <- as.data.frame(ggpredict(m4, terms = c("sch_ses [-0.77:1.49 by=0.02]", "male", "pped")))
lp <- expand_grid(sch_ses = seq(-0.77, 1.49, .02), male = factor(c(0, 1)), pped = c(0, 1))
lp$fit <- predict(m2, lp)
lp <- lp |> mutate(group = male, facet = factor(pped))
ge <- ge |> mutate(facet = factor(facet))
pp_lab <- c("0" = "no preschool (pped = 0)", "1" = "preschool (pped = 1)")
p6 <- ggplot(ge, aes(x, predicted, color = group, fill = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = .15, color = NA) +
  geom_line(linewidth = 1.2) +
  geom_line(data = lp, aes(sch_ses, fit, color = group), linetype = "dashed", linewidth = .8, inherit.aes = FALSE) +
  facet_wrap(~ facet, labeller = as_labeller(pp_lab)) +
  scale_color_manual(values = SEX, labels = sex_lab, name = NULL) +
  scale_fill_manual(values = SEX, labels = sex_lab, name = NULL) +
  scale_y_continuous(breaks = seq(.05, .25, .05), labels = scales::percent_format(accuracy = 1)) +
  labs(x = "school SES (sch_ses)", y = "predicted probability of retention",
       title = "The logit model m4 on the probability scale",
       subtitle = "Solid: m4, with 95% confidence bands. Dashed: the LPM m2")
save(p6, "fig06_logit_probability.png", 9.6, 5.0)

## fig06b --- research question 3 on two scales (beyond the lecture) ----------------------
nd <- expand_grid(male = factor(c(0, 1)), pped = c(0, 1), sch_ses = 0)
nd$lo4 <- predict(m4, nd); nd$p4 <- plogis(nd$lo4); nd$p3 <- predict(m3, nd)
nd <- nd |> mutate(x = factor(pped, labels = c("no preschool", "preschool")))
labL <- nd |> group_by(male) |> summarise(y = mean(lo4), d = diff(lo4))
labR <- nd |> group_by(male) |> summarise(y = mean(p4), d = diff(p4))
p6ba <- ggplot(nd, aes(x, lo4, color = male, group = male)) +
  geom_line(linewidth = 1.2) + geom_point(size = 2.6) +
  geom_text(data = labL, aes(x = 1.5, y = y + .07, label = sprintf("%.2f", d)), inherit.aes = FALSE,
            color = c(BLU, ORANGE), size = 3.8, vjust = 0) +
  scale_color_manual(values = SEX, labels = sex_lab, name = NULL) +
  labs(x = NULL, y = "predicted log-odds of retention", title = "Log-odds scale (m4)",
       subtitle = "Parallel lines: the same preschool difference")
p6bb <- ggplot(nd, aes(x, p4, color = male, group = male)) +
  geom_line(aes(y = p3), linetype = "dashed", linewidth = .8) +
  geom_line(linewidth = 1.2) + geom_point(size = 2.6) +
  geom_text(data = labR, aes(x = c(1.45, 1.55), y = y + c(-.012, .012), label = sprintf("%.1f points", 100 * d)), inherit.aes = FALSE,
            color = c(BLU, ORANGE), size = 3.8, hjust = c(1, 0)) +
  scale_color_manual(values = SEX, labels = sex_lab, name = NULL, guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "predicted probability of retention", title = "Probability scale",
       subtitle = "Solid: m4. Dashed: the LPM m3, with an interaction")
save((p6ba + p6bb) + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
     "fig06b_rq3_two_scales.png", 10.4, 4.8)

## fig07 --- tau on the log-odds and probability scales (beyond the lecture) ----------------
uu <- tibble(z = seq(-3, 3, .01), lo = g0 + z * tau, dens = dnorm(z))
marks <- tibble(z = -2:2, lo = g0 + z * tau, p = plogis(lo)) |>
  mutate(lab_lo = sprintf("%.2f", lo), lab_p = pct(p),
         zlab = c("-2 SD", "-1 SD", "average\nschool", "+1 SD", "+2 SD"))
p7a <- ggplot(uu, aes(lo, dens)) +
  geom_area(fill = BLU, alpha = .18) + geom_line(color = BLU, linewidth = 1) +
  geom_segment(data = marks, aes(x = lo, xend = lo, y = 0, yend = dnorm(z)), color = INK, linewidth = .5) +
  geom_label(data = marks, aes(x = lo + c(-.25, -.2, 0, .2, .25), y = dnorm(z) + .04, label = lab_lo,
                               hjust = c(1, 1, .5, 0, 0)), size = 3.6, fill = "white", linewidth = 0, label.padding = unit(.1, "lines")) +
  geom_text(data = marks, aes(x = lo, y = -.04, label = zlab), size = 3.6, color = MUT, lineheight = .85) +
  scale_y_continuous(limits = c(-.07, .46)) +
  labs(x = "school's log-odds of retention for the reference girl", y = "density (higher = more schools)",
       title = "Schools on the log-odds scale", subtitle = "Normal: mean -2.24, SD 1.27") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
pd <- tibble(p = seq(.0005, .85, .0005)) |> mutate(lo = qlogis(p), dens = dnorm((lo - g0) / tau) / tau / (p * (1 - p)))
p7b <- ggplot(pd, aes(p, dens)) +
  geom_area(fill = ORANGE, alpha = .18) + geom_line(color = ORANGE, linewidth = 1) +
  geom_vline(data = marks, aes(xintercept = p), color = INK, linewidth = .4, linetype = "dashed") +
  geom_label(data = marks, aes(x = p, y = c(8.1, 7.5, 5.2, 2.6, 1.6), label = lab_p), size = 3.6, hjust = -0.05,
             fill = "white", linewidth = 0, label.padding = unit(.1, "lines")) +
  coord_cartesian(xlim = c(0, .8), ylim = c(0, 8.5)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "school's probability of retention for the reference girl", y = "density (higher = more schools)",
       title = "The same schools on the probability scale", subtitle = "Dashed lines: the five schools marked on the left") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
save(p7a + p7b, "fig07_tau_two_scales.png", 10.8, 4.8)

## fig08 --- caterpillar plot (lecture PDF 27), 90% intervals as in the lecture ---------------
# Each school's own SE. (augment(ranef(m5)) with broom.mixed 0.2.9.7 sorts the estimates but
# leaves std.error/lb/ub in the original school order, so its intervals belong to other schools.)
cat_df <- tidy(m5, effects = "ran_vals", conf.int = TRUE, conf.level = .90) |>
  arrange(estimate) |> mutate(level = factor(level, levels = level))
cat("fig08: above", sum(cat_df$conf.low > 0), " below", sum(cat_df$conf.high < 0), " incl0", sum(cat_df$conf.low <= 0 & cat_df$conf.high >= 0), "\n")
p8 <- ggplot(cat_df, aes(x = level, y = estimate, ymin = conf.low, ymax = conf.high)) +
  geom_errorbar(width = 0, alpha = .3) +
  geom_point(size = .9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "school (356, sorted by estimate)", y = expression(paste("EB estimate of ", u[j], " (log-odds)")),
       title = "Caterpillar plot of the school intercepts, m5",
       subtitle = "Each point is a school's estimated departure from the average school, with a 90% interval") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid.major.x = element_blank())
save(p8, "fig08_caterpillar.png", 8.8, 4.6)

## fig09 --- school intercepts against observed retention rates (beyond the lecture) ---------
re <- ranef(m5)$schoolid |> rownames_to_column("schoolid") |> rename(u = `(Intercept)`) |>
  mutate(schoolid = as.numeric(schoolid))
sch <- thai |> group_by(schoolid) |> summarise(n = n(), r = sum(retained), p = mean(retained)) |>
  left_join(re, by = "schoolid") |> mutate(none = r == 0)
p9a <- ggplot(sch, aes(p, u)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = MUT) +
  geom_point(aes(size = n, color = none), alpha = .55) +
  scale_color_manual(values = c("FALSE" = BLU, "TRUE" = RED), guide = "none") +
  scale_size_continuous(range = c(.8, 4.5), name = "students\nin school") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "percent of the school's students retained", y = expression(paste("EB estimate of ", u[j], " (log-odds)")),
       title = "All 356 schools", subtitle = "Red: schools with no student retained")
zz9 <- sch |> filter(none)
ends <- bind_rows(zz9 |> slice_min(n, n = 1, with_ties = FALSE), zz9 |> slice_max(n, n = 1, with_ties = FALSE))
p9b <- ggplot(zz9, aes(n, u)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = MUT) +
  geom_point(color = RED, alpha = .6, size = 2) +
  geom_label(data = ends, aes(label = sprintf("%d students: %.2f", n, u)), size = 3.5, color = INK, fill = alpha("white", .7), linewidth = 0, label.padding = unit(.1, "lines"),
            hjust = c(0, 1), nudge_x = c(.8, -.8)) +
  labs(x = "number of students in the school", y = NULL,
       title = "90 schools with no student retained", subtitle = "More students: estimate farther below 0")
save(p9a + p9b + plot_layout(widths = c(1.3, 1)), "fig09_shrinkage.png", 10.8, 4.9)

## fig10 --- random intercepts vs random slopes for male (beyond the lecture) ---------------
set.seed(11)
pick <- sample(unique(thai$schoolid), 40)
r5 <- ranef(m5)$schoolid |> rownames_to_column("schoolid") |> filter(schoolid %in% pick) |>
  transmute(schoolid, u0 = `(Intercept)`, u1 = 0, model = "m5: random intercepts only")
r6 <- ranef(m6)$schoolid |> rownames_to_column("schoolid") |> filter(schoolid %in% pick) |>
  transmute(schoolid, u0 = `(Intercept)`, u1 = male1, model = "m6: + random slope for male")
f5 <- fixef(m5); f6 <- fixef(m6)
seg10 <- bind_rows(
  r5 |> mutate(girl = f5[1] + u0, boy = f5[1] + f5[2] + u0),
  r6 |> mutate(girl = f6[1] + u0, boy = f6[1] + f6[2] + u0 + u1)) |>
  pivot_longer(c(girl, boy), names_to = "sex", values_to = "lo") |>
  mutate(x = ifelse(sex == "girl", 0, 1))
avg10 <- tibble(model = c("m5: random intercepts only", "m6: + random slope for male"),
                girl = c(f5[1], f6[1]), boy = c(f5[1] + f5[2], f6[1] + f6[2])) |>
  pivot_longer(c(girl, boy), names_to = "sex", values_to = "lo") |> mutate(x = ifelse(sex == "girl", 0, 1))
cols40 <- setNames(scales::hue_pal(l = 55)(40), as.character(sort(as.numeric(pick))))
p10 <- ggplot(seg10, aes(x, lo, group = schoolid, color = schoolid)) +
  geom_line(alpha = .85) +
  geom_line(data = avg10, aes(x, lo, group = model), color = INK, linewidth = 1.6, inherit.aes = FALSE) +
  facet_wrap(~ model) +
  scale_color_manual(values = cols40, guide = "none") +
  scale_x_continuous(breaks = c(0, 1), labels = c("girl", "boy"), expand = expansion(.12)) +
  labs(x = NULL, y = "predicted log-odds of retention",
       title = "Forty schools' predicted log-odds for a girl and a boy",
       subtitle = "No preschool, school SES = 0. Black: the average school. Each line's slope is that school's gender gap")
save(p10, "fig10_random_slopes.png", 9.8, 5.0)

## fig11 --- school proportion male (lecture PDF 32) -------------------------------------
sp <- thai |> distinct(schoolid, prop_male) |> mutate(edge = prop_male %in% c(0, 1))
p11 <- ggplot(sp, aes(prop_male, fill = edge)) +
  geom_histogram(binwidth = .025, boundary = 0, color = "white") +
  scale_fill_manual(values = c("FALSE" = MUT, "TRUE" = RED), guide = "none") +
  annotate("text", x = .03, y = 12, label = "1 school with\nonly girls", color = RED, size = 3.6, hjust = 0) +
  annotate("text", x = .97, y = 12, label = "6 schools with\nonly boys", color = RED, size = 3.6, hjust = 1) +
  labs(x = "proportion of the school's students who are boys (prop_male)", y = "number of schools",
       title = "School proportion male across the 356 schools",
       subtitle = "Mean 0.50, SD 0.15 across schools")
save(p11, "fig11_prop_male.png", 8.6, 4.4)

## fig11b --- the contextual effect as the gap between two slopes in m7 (after Topic 3's R&B figure) --
f7 <- fixef(m7); bW <- f7[["male1"]]; bC <- f7[["prop_male"]]; b0 <- f7[["(Intercept)"]]; bB <- bW + bC
pA <- .25; pB <- .75                                            # two schools: 25% and 75% boys
yy <- function(p, x) b0 + bC * p + bW * x                       # no preschool, SES 0, u_j = 0
wl <- tibble(sch = rep(c("A", "B"), each = 2), p = rep(c(pA, pB), each = 2), x = rep(0:1, 2)) |> mutate(y = yy(p, x))
av <- tibble(x = c(pA, pB), y = b0 + bB * c(pA, pB))            # each school's average
gp <- tibble(x = 0:1, y = yy(pA, 0:1), yend = yy(pB, 0:1))      # same gender, the two schools
glab <- sprintf("beta[C] %%*%% 0.5 == %.2f", bC * (pB - pA))
cat("fig11b: A", round(yy(pA, 0:1), 4), " B", round(yy(pB, 0:1), 4), " averages", round(av$y, 4), " gap", round(bC * (pB - pA), 4), "\n")
p11b <- ggplot() +
  annotate("segment", x = -.08, xend = 1.08, y = b0 - .08 * bB, yend = b0 + 1.08 * bB, color = RED, linewidth = 1.1) +
  geom_line(data = wl, aes(x, y, group = sch), color = GRN, linewidth = 1.2) +
  geom_segment(data = gp, aes(x = x, xend = x, y = y, yend = yend), color = BLU, linewidth = 1.4) +
  geom_point(data = wl, aes(x, y), color = GRN, size = 2.6) +
  geom_point(data = av, aes(x, y), color = RED, size = 3.4) +
  annotate("text", x = 1.05, y = yy(pA, 1), label = "School A: 25% boys", hjust = 0, color = GRN, size = 3.8) +
  annotate("text", x = 1.05, y = yy(pB, 1), label = "School B: 75% boys", hjust = 0, color = GRN, size = 3.8) +
  annotate("text", x = 1.05, y = mean(yy(c(pA, pB), 1)) + .03, label = "two boys:", hjust = 0, color = BLU, size = 3.8) +
  annotate("text", x = 1.05, y = mean(yy(c(pA, pB), 1)) - .03, label = glab, parse = TRUE, hjust = 0, color = BLU, size = 3.8) +
  annotate("text", x = -.05, y = mean(yy(c(pA, pB), 0)) + .03, label = "two girls:", hjust = 1, color = BLU, size = 3.8) +
  annotate("text", x = -.05, y = mean(yy(c(pA, pB), 0)) - .03, label = glab, parse = TRUE, hjust = 1, color = BLU, size = 3.8) +
  annotate("text", x = -.11, y = b0 - .11 * bB, label = "slope~beta[B] == 0.06", parse = TRUE, hjust = 1, color = RED, size = 3.8) +
  annotate("text", x = .40, y = yy(pA, .40) + .10, label = "slope~beta[W] == 0.54", parse = TRUE, hjust = 1, color = GRN, size = 3.8) +
  annotate("text", x = pA, y = av$y[1] - .06, label = "A's average", vjust = 1, color = RED, size = 3.4) +
  annotate("text", x = pB, y = av$y[2] + .06, label = "B's average", vjust = 0, color = RED, size = 3.4) +
  scale_x_continuous(breaks = c(0, .25, .5, .75, 1), labels = c("0\ngirl", "0.25", "0.50", "0.75", "1\nboy")) +
  coord_cartesian(xlim = c(-.42, 1.42)) +
  labs(x = "male (0 = girl, 1 = boy); a school's average of male is its proportion of boys", y = "predicted log-odds of retention",
       title = "The contextual effect is the gap between the two slopes",
       subtitle = "Within a school: slope 0.54. Between school averages: slope 0.06. Gap: 0.06 \u2212 0.54 = \u22120.48.")
save(p11b, "fig11b_mundlak_lines.png", 9.0, 5.2)

## fig12 --- the lecture's simulation, with a seed (lecture PDF 34) -------------------------
set.seed(34)
n <- 10000; j <- 20
l2 <- tibble(cluster = 1:j, u_j = rnorm(j, 0, 2))
dat <- tibble(cluster = sample(1:j, n, replace = TRUE), x_ij = runif(n, -4, 4)) |>
  left_join(l2, by = "cluster") |>
  mutate(y_ij = x_ij + u_j, p_ij = plogis(y_ij))
grid12 <- tibble(x_ij = seq(-4, 4, .05))
avgc <- grid12 |> mutate(p = sapply(x_ij, function(x) mean(plogis(x + l2$u_j))))
cs_curve <- grid12 |> mutate(p = plogis(x_ij), y = x_ij)
p12a <- ggplot(dat, aes(x_ij, y_ij)) +
  geom_line(aes(group = cluster), alpha = .45, color = MUT) +
  geom_smooth(method = lm, se = FALSE, formula = y ~ x, color = BLU, linewidth = 1.3) +
  geom_line(data = cs_curve, aes(x_ij, y), color = RED, linetype = "dashed", linewidth = 1) +
  labs(x = "X", y = "Y", title = "Continuous outcome",
       subtitle = expression(paste("Blue: one line through all the data. Red dashed: ", u[j], " = 0")))
p12b <- ggplot(dat, aes(x_ij, p_ij)) +
  geom_line(aes(group = cluster), alpha = .45, color = MUT) +
  geom_smooth(se = FALSE, method = "glm", formula = y ~ x, method.args = list(family = "quasibinomial"),
              color = BLU, linewidth = 1.3) +
  geom_line(data = avgc, aes(x_ij, p), color = INK, linetype = "dotted", linewidth = 1.1) +
  geom_line(data = cs_curve, aes(x_ij, p), color = RED, linetype = "dashed", linewidth = 1) +
  labs(x = "X", y = "probability that Y = 1", title = "Dichotomous outcome",
       subtitle = "Black dotted: the average of the 20 curves")
save(p12a + p12b, "fig12_simulation.png", 10.8, 4.8)

## fig13 --- cluster-specific vs population-averaged with our numbers (beyond the lecture) --
zz <- tibble(z = seq(-2.5, 2.5, .01)) |>
  mutate(girl_lo = g0 + z * tau, boy_lo = g0 + bm + z * tau,
         girl = plogis(girl_lo), boy = plogis(boy_lo))
marg <- function(eta) integrate(function(u) plogis(eta + u) * dnorm(u, 0, tau), -Inf, Inf)$value
pa <- c(girl = marg(g0), boy = marg(g0 + bm))
cat("population-averaged probabilities (reference values):", round(pa, 4), " log-odds:", round(qlogis(pa), 3), "\n")
lo_long <- zz |> select(z, girl = girl_lo, boy = boy_lo) |> pivot_longer(-z, names_to = "sex", values_to = "v")
pr_long <- zz |> select(z, girl, boy) |> pivot_longer(-z, names_to = "sex", values_to = "v")
SX2 <- c(girl = BLU, boy = ORANGE); sx2_lab <- c(girl = "girls", boy = "boys")
xlab13 <- expression(paste("school's ", u[j], " in SDs (0 = average school)"))
p13a <- ggplot(lo_long, aes(z, v, color = sex)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = qlogis(pa), color = SX2[c("girl", "boy")], linetype = "dashed", linewidth = .8) +
  annotate("segment", x = 1.5, xend = 1.5, y = g0 + 1.5 * tau, yend = g0 + bm + 1.5 * tau, color = INK,
           arrow = arrow(length = unit(.1, "cm"), ends = "both")) +
  annotate("text", x = -0.2, y = 1.0, label = "solid: gap 0.53 in every school", size = 3.6) +
  annotate("text", x = -2.5, y = qlogis(pa[["boy"]]) + .75, label = "dashed: averaged over\nschools, gap 0.435", lineheight = .9,
           hjust = 0, size = 3.4) +
  scale_color_manual(values = SX2, labels = sx2_lab, breaks = c("girl", "boy"), name = NULL) +
  labs(x = xlab13, y = "predicted log-odds", title = "Log-odds scale",
       subtitle = "Solid: each school. Dashed: averaged over schools")
pts <- tibble(z = 0, v = c(plogis(g0), plogis(g0 + bm)), sex = c("girl", "boy"))
p13b <- ggplot(pr_long, aes(z, v, color = sex)) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = pa, color = SX2[c("girl", "boy")], linetype = "dashed", linewidth = .8) +
  geom_point(data = pts, size = 3) +
  annotate("text", x = .12, y = plogis(g0) - .035, label = sprintf("average school: %.3f", plogis(g0)), hjust = 0, size = 3.3, color = BLU) +
  annotate("text", x = -.12, y = .178, label = sprintf("average school: %.3f", plogis(g0 + bm)), hjust = 1, size = 3.3, color = ORANGE) +
  annotate("text", x = -2.5, y = pa[["boy"]] + .03, label = sprintf("boys, averaged: %.3f", pa[["boy"]]),
           hjust = 0, size = 3.3, color = ORANGE) +
  annotate("text", x = -2.5, y = pa[["girl"]] - .03, label = sprintf("girls, averaged: %.3f", pa[["girl"]]),
           hjust = 0, size = 3.3, color = BLU) +
  scale_color_manual(values = SX2, labels = sx2_lab, breaks = c("girl", "boy"), name = NULL, guide = "none") +
  scale_y_continuous(limits = c(0, .85)) +
  labs(x = xlab13, y = "predicted probability", title = "Probability scale",
       subtitle = "The gap differs from school to school")
save((p13a + p13b) + plot_layout(guides = "collect") & theme(legend.position = "bottom"),
     "fig13_cs_vs_pa.png", 10.8, 5.0)

## fig14 --- the latent-variable view (appendix; beyond the lecture) -----------------------
lat <- expand_grid(y = seq(-10, 7, .02), sch = c("average school", "school 1 SD above average")) |>
  mutate(center = ifelse(startsWith(sch, "average"), g0, g0 + tau), d = dlogis(y, center, 1))
areas <- lat |> distinct(sch, center) |> mutate(p = plogis(center), lab = sprintf("area above 0 = %.3f", p),
                                                clab = sprintf("center = %.2f", center))
p14 <- ggplot(lat, aes(y, d)) +
  geom_area(data = filter(lat, y >= 0), fill = RED, alpha = .35) +
  geom_line(linewidth = 1, color = INK) +
  geom_vline(xintercept = 0, color = RED, linewidth = .8) +
  geom_vline(aes(xintercept = center), linetype = "dashed", color = MUT) +
  geom_text(data = areas, aes(x = 0.4, y = .2, label = lab), hjust = 0, size = 3.8, color = RED) +
  geom_text(data = areas, aes(x = center - .15, y = .02, label = clab), hjust = 1, size = 3.6, color = MUT) +
  facet_wrap(~ sch, ncol = 1) +
  labs(x = "latent propensity to be retained, y*", y = "density (how common each value is)",
       title = "The latent-variable view of the multilevel logit model",
       subtitle = "Reference girl. Retained when y* is above 0 (red line). Dashed line: the school's log-odds")
save(p14, "fig14_latent.png", 8.8, 5.6)

cat("done\n")
