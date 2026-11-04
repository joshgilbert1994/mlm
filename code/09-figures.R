# Figures for Chapter 9: Higher-Level Models and Meta-Analysis.
# Written by Bill Young. Run from the root of the book repository; set DATA
# to wherever your copy of relyea.csv lives. Figures are written to OUT.
# Packages: tidyverse, lme4, broom.mixed, patchwork. (meta, metafor, moderndive and latex2exp are not
# installed: the meta-analysis is computed by hand below, geom_parallel_slopes() is replaced by two lines
# with the common slope from lm(bgk_score ~ treat + t_pretest), and every symbol is plotmath.)
suppressMessages({library(tidyverse); library(lme4); library(broom.mixed); library(patchwork)})
DATA <- "data/relyea.csv"
OUT  <- "chapters/images/09/"; dir.create(OUT, showWarnings = FALSE)

# The data, read exactly as the lecture reads it
relyea <- read_csv(DATA, show_col_types = FALSE) |>
  mutate(treat = factor(treat)) |>
  select(id, cluster_id, block_id, treat, bgk_score, t_pretest)

INK<-"#1A1D23"; BLU<-"#2C5F8A"; RED<-"#B03A2E"; GRN<-"#2E7D5B"; MUT<-"#5C6672"; RUL<-"#DDE2E6"
ORANGE <- "#D98E04"
ARM <- c("adaptive (treat = 1)" = BLU, "core (treat = 0)" = ORANGE)
arm_lab <- function(t) factor(ifelse(t == "1", "adaptive (treat = 1)", "core (treat = 0)"), levels = names(ARM))
th <- theme_bw(base_size = 14) + theme(
  panel.grid.minor = element_blank(),
  panel.grid.major = element_line(color = RUL, linewidth = .3),
  plot.title = element_text(face = "bold", size = 14, color = INK),
  plot.subtitle = element_text(size = 11.5, color = MUT),
  axis.title = element_text(size = 11.5, color = INK),
  axis.text = element_text(size = 11.2, color = INK),
  strip.background = element_rect(fill = "grey94", color = NA),
  strip.text = element_text(face = "bold", size = 11.2, color = INK),
  legend.position = "bottom", legend.title = element_text(size = 11), legend.text = element_text(size = 11.2),
  plot.background = element_rect(fill = "white", color = NA))
theme_set(th)
f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)
# minus sign for printed numbers (plain text): use the hyphen-minus; plotmath draws its own minus

# Canvas sizes. Each figure is drawn on a canvas about 1.3 times its printed width, so the 11.2-pt
# tick labels print at about 8.5 pt. PCT is the width to use in the qmd (percent of the 6.5 in text
# block); FAC rescales a canvas (width and height together) if a figure is later shrunk or enlarged in
# the qmd: change FAC by the same ratio as PCT so the printed text size stays the same.
PCT <- c(fig01 = 78, fig02 = 80, fig03 = 80, fig04 = 72, fig05 = 78, fig06 = 78, fig07 = 80,
         fig08 = 64, fig09 = 66, fig10 = 66, fig11 = 80, fig12 = 70, fig13 = 80)
FAC <- c(fig01 = 1, fig02 = 1, fig03 = 1, fig04 = 1, fig05 = 1, fig06 = 1, fig07 = 1,
         fig08 = 64/72, fig09 = 1, fig10 = 1, fig11 = 1, fig12 = 1, fig13 = 1)
save <- function(p, f, w, h) {
  k <- substr(f, 1, 5); sc <- FAC[[k]]
  pt <- 11.2 * (PCT[[k]] / 100 * 6.5) / (w * sc)
  cat(sprintf("### %s canvas %.2f x %.2f in; printed at %d%% (%.2f in): tick labels %.1f pt\n",
              f, w * sc, h * sc, PCT[[k]], PCT[[k]] / 100 * 6.5, pt))
  ggsave(paste0(OUT, f), p, width = w * sc, height = h * sc, dpi = 300, bg = "white") }

# ---- sample sizes and the lecture's models ----------------------------------------------
cat("N students", nrow(relyea), " classrooms", n_distinct(relyea$cluster_id), " schools", n_distinct(relyea$block_id), "\n")
ols <- lm(bgk_score ~ treat + t_pretest, relyea)
m1 <- lmer(bgk_score ~ 1 + (1 | cluster_id) + (1 | block_id), data = relyea)
m2 <- suppressMessages(lmer(bgk_score ~ 1 + treat + t_pretest + (1 | cluster_id) + (1 | block_id), data = relyea))
m3 <- suppressMessages(lmer(bgk_score ~ 1 + treat + t_pretest + (1 | cluster_id) + (1 + treat | block_id), data = relyea))
print(round(coef(ols), 4))
vc1 <- as.data.frame(VarCorr(m1))
psi2 <- vc1$vcov[vc1$grp == "block_id"]; tau2 <- vc1$vcov[vc1$grp == "cluster_id"]; sig2 <- vc1$vcov[vc1$grp == "Residual"]
tot <- psi2 + tau2 + sig2
cat(sprintf("m1: psi2 %.4f tau2 %.4f sig2 %.4f total %.4f; ICC same classroom %.4f; same school %.4f\n",
            psi2, tau2, sig2, tot, (psi2 + tau2) / tot, psi2 / tot))
ate <- fixef(m3)[["treat1"]]
cat(sprintf("m3: average treatment effect %.4f\n", ate))

# Per-school regressions, as the lecture's agg pipeline (lm(bgk_score ~ treat) within each school)
schools <- sort(unique(relyea$block_id))
agg <- map_dfr(schools, function(s) {
  dat <- filter(relyea, block_id == s)
  fit <- lm(bgk_score ~ treat, dat)
  tt <- tidy(fit, conf.int = TRUE) |> filter(term == "treat1")
  tibble(block_id = s, nobs = nrow(dat), n_T = sum(dat$treat == "1"), n_C = sum(dat$treat == "0"),
         mean_T = mean(dat$bgk_score[dat$treat == "1"]), mean_C = mean(dat$bgk_score[dat$treat == "0"]),
         estimate = tt$estimate, std.error = tt$std.error, lm_lo = tt$conf.low, lm_hi = tt$conf.high,
         sigma = glance(fit)$sigma)
}) |>
  mutate(d = estimate / sigma,
         se_d_std = sqrt((n_T + n_C) / (n_T * n_C) + d^2 / (2 * (n_T + n_C))))   # standard large-sample SE of d
cat(sprintf("mean(d) %.8f var(d) %.7f; standard SE(d) range %.3f to %.3f; raw differences %.3f to %.3f\n",
            mean(agg$d), var(agg$d), min(agg$se_d_std), max(agg$se_d_std), min(agg$estimate), max(agg$estimate)))

# ---- the PRINTED meta-analysis: metagen(d, se_d_std) on Cohen's d with its standard SE, by hand ---
TE <- agg$d; v <- agg$se_d_std^2
rll <- function(t2) { w <- 1 / (v + t2); mu <- sum(w * TE) / sum(w)          # REML log-likelihood of tau^2
  -0.5 * (sum(log(v + t2)) + log(sum(w)) + sum(w * (TE - mu)^2)) }
opt <- optimize(rll, c(0, 10), maximum = TRUE, tol = 1e-12)
t2 <- if (rll(0) >= opt$objective) 0 else opt$maximum
wf <- 1 / v; muf <- sum(wf * TE) / sum(wf); sef <- 1 / sqrt(sum(wf))
wr <- 1 / (v + t2); mur <- sum(wr * TE) / sum(wr); ser <- 1 / sqrt(sum(wr))
Q <- sum(wf * (TE - muf)^2); pQ <- pchisq(Q, 25, lower.tail = FALSE); I2 <- (Q - 25) / Q
agg <- agg |> mutate(w_common = 100 * wf / sum(wf), w_random = 100 * wr / sum(wr),
                     ci_lo = d - 1.96 * se_d_std, ci_hi = d + 1.96 * se_d_std)
cat(sprintf("REML tau^2 %.4f (tau %.4f); common %.4f [%.4f, %.4f]; random %.4f [%.4f, %.4f]; Q %.2f p %.4f; I2 %.1f%%\n",
            t2, sqrt(t2), muf, muf - 1.96 * sef, muf + 1.96 * sef, mur, mur - 1.96 * ser, mur + 1.96 * ser, Q, pQ, 100 * I2))
cat(sprintf("weights: common %.2f%% to %.2f%%; random %.2f%% to %.2f%%\n",
            min(agg$w_common), max(agg$w_common), min(agg$w_random), max(agg$w_random)))
stopifnot(abs(t2 - 0.0290) < 5e-5, abs(mur - 0.0567) < 5e-5, abs(muf - 0.0604) < 5e-5,
          abs(mur - 1.96 * ser + 0.0584) < 5e-5, abs(mur + 1.96 * ser - 0.1717) < 5e-5)

## fig01 --- the design: three of the 26 schools ------------------------------------------------
show_sch <- c(22, 5, 3)
cls <- relyea |> filter(block_id %in% show_sch) |> count(block_id, cluster_id, treat) |>
  mutate(block_id = factor(block_id, levels = show_sch)) |> arrange(block_id, desc(treat), cluster_id)
TW <- 1; TG <- .22; PADX <- .28; SG <- .75; NC <- 4; DY <- .21
cls <- cls |> group_by(block_id) |> mutate(pos = row_number()) |> ungroup()
sch_w <- cls |> count(block_id, name = "ncl") |> mutate(w = ncl * TW + (ncl - 1) * TG + 2 * PADX,
                                                        x0 = cumsum(lag(w + SG, default = 0)), x1 = x0 + w)
cls <- cls |> left_join(sch_w, by = "block_id") |>
  mutate(xmin = x0 + PADX + (pos - 1) * (TW + TG), xmax = xmin + TW,
         arm = arm_lab(treat))
TH <- ceiling(max(cls$n) / NC) * DY + .2
stud <- cls |> rowwise() |> reframe(block_id, cluster_id, arm, k = seq_len(n), xmin) |>
  mutate(x = xmin + .17 + ((k - 1) %% NC) * ((TW - .34) / (NC - 1)),
         y = 0.18 + ((k - 1) %/% NC) * DY)
sch_n <- relyea |> count(block_id) |> filter(block_id %in% show_sch) |>
  mutate(block_id = factor(block_id, levels = show_sch)) |> left_join(sch_w, by = "block_id")
p1 <- ggplot() +
  geom_rect(data = sch_n, aes(xmin = x0, xmax = x1, ymin = -.2, ymax = TH + .2), fill = "grey93", color = "grey70") +
  geom_text(data = sch_n, aes(x = (x0 + x1) / 2, y = TH + .62, label = paste0("school ", block_id)),
            size = 4, color = INK, fontface = "bold") +
  geom_text(data = sch_n, aes(x = (x0 + x1) / 2, y = TH + .36, label = paste0(n, " students")),
            size = 3.8, color = INK) +
  geom_rect(data = cls, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = TH, fill = arm), alpha = .25, color = NA) +
  geom_rect(data = cls, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = TH, color = arm), fill = NA, linewidth = .7) +
  geom_point(data = stud, aes(x, y, color = arm), size = .9) +
  scale_fill_manual(values = ARM, name = NULL) + scale_color_manual(values = ARM, name = NULL) +
  scale_x_continuous(expand = expansion(add = .1)) + scale_y_continuous(expand = expansion(add = c(.05, .15))) +
  labs(title = "Three of the 26 schools in the trial",
       subtitle = "Gray box: a school (level 3); colored box: a classroom (level 2); dot: a student (level 1)") +
  theme_void(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 14, color = INK),
        plot.subtitle = element_text(size = 11.2, color = MUT),
        legend.position = "bottom", legend.text = element_text(size = 11.2),
        plot.margin = margin(6, 6, 6, 6))
save(p1, "fig01_design.png", 6.6, 3.3)
print(cls |> select(block_id, cluster_id, treat, n))

## fig02 --- BGK scores: distribution and against the classroom pretest mean (lecture PDF 9-10) -------
b <- coef(ols)
lines2 <- relyea |> group_by(treat) |> summarise(x0 = min(t_pretest), x1 = max(t_pretest)) |>
  mutate(a = b[["(Intercept)"]] + ifelse(treat == "1", b[["treat1"]], 0), y0 = a + b[["t_pretest"]] * x0,
         y1 = a + b[["t_pretest"]] * x1, arm = arm_lab(treat))
p2a <- ggplot(relyea, aes(bgk_score)) +
  geom_histogram(bins = 30, fill = "grey55", color = "white", linewidth = .2) +
  labs(x = "BGK score (SD units)", y = "number of students", title = "BGK scores",
       subtitle = "1,914 students")
p2b <- ggplot(relyea |> mutate(arm = arm_lab(treat)), aes(t_pretest, bgk_score, color = arm)) +
  geom_point(alpha = .25, size = .8) +
  geom_segment(data = lines2, aes(x = x0, xend = x1, y = y0, yend = y1, color = arm), linewidth = 1.1) +
  scale_color_manual(values = ARM, name = NULL) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2.5, linewidth = 1.1))) +
  labs(x = "classroom mean MAP pretest (SD units)", y = "BGK score (SD units)",
       title = "BGK and classroom pretest", subtitle = "Lines: single-level regression (OLS)")
save(p2a + p2b + plot_layout(widths = c(1, 1.25), guides = "collect") & theme(legend.position = "bottom"),
     "fig02_bgk.png", 6.8, 3.45)

## fig03 --- each school is a small randomized trial (lecture PDF 11 and 31) --------------------
p3 <- ggplot(relyea |> mutate(arm = factor(arm_lab(treat), levels = rev(names(ARM))),
                              sch = factor(paste("school", block_id), levels = paste("school", schools))),
             aes(arm, bgk_score, fill = arm)) +
  geom_boxplot(outlier.size = .5, linewidth = .35, width = .7, color = INK, alpha = .85) +
  facet_wrap(~ sch, ncol = 7) +
  scale_fill_manual(values = ARM, name = NULL, breaks = rev(names(ARM)),
                    labels = c("core (treat = 0)" = "core\n(treat = 0)", "adaptive (treat = 1)" = "adaptive\n(treat = 1)")) +
  scale_y_continuous(breaks = c(-2, 0, 2)) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = NULL, y = "BGK score (SD units)") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.spacing = unit(3, "pt"), strip.text = element_text(size = 10.5, margin = margin(2, 0, 2, 0)),
        legend.position = "inside", legend.position.inside = c(.995, .12), legend.justification = c(1, .5),
        legend.key.spacing.y = unit(6, "pt"),
        legend.background = element_blank())
save(p3, "fig03_schools.png", 6.8, 5.1)

## fig04 --- the data-generating process, the lecture's picture (PDF 16) with plotmath -------------
delta_0 <- 100; r_k <- 10; u_jk <- -5; e_ijk <- 3                       # the lecture's values
school_mean <- delta_0 + r_k; class_mean <- school_mean + u_jk; student_score <- class_mean + e_ijk
sd_school <- 15; sd_class <- 10; sd_student <- 7; plot_range <- c(50, 150)
max_y <- 1 / (min(sd_school, sd_class, sd_student) * sqrt(2 * pi))
xx <- seq(plot_range[1], plot_range[2], length.out = 600)
cur <- bind_rows(
  tibble(x = xx, y = dnorm(xx, delta_0, sd_school), lev = "school"),
  tibble(x = xx, y = dnorm(xx, school_mean, sd_class), lev = "classroom"),
  tibble(x = xx, y = dnorm(xx, class_mean, sd_student), lev = "student")) |>
  mutate(lev = factor(lev, levels = c("school", "classroom", "student")))
LCOL <- c(school = GRN, classroom = ORANGE, student = BLU)
LTY <- c(school = "solid", classroom = "dashed", student = "solid")
arr <- arrow(length = unit(0.16, "cm"), type = "closed")
p4 <- ggplot(cur, aes(x, y)) +
  geom_area(aes(fill = lev, group = lev), alpha = .07, position = "identity") +
  geom_line(aes(color = lev, linetype = lev), linewidth = 1) +
  geom_vline(xintercept = delta_0, color = GRN, alpha = .6) +
  geom_vline(xintercept = school_mean, color = ORANGE, linetype = "dashed", alpha = .7) +
  geom_vline(xintercept = class_mean, color = BLU, linetype = "dotted", alpha = .8) +
  annotate("segment", x = delta_0, xend = school_mean, y = max_y * .70, yend = max_y * .70, arrow = arr, color = GRN, linewidth = .8) +
  annotate("label", x = (delta_0 + school_mean) / 2, y = max_y * .74, vjust = 0, label = 'r[k] == "+10"', parse = TRUE,
           size = 4, fill = "white", border.colour = NA, label.padding = unit(1, "pt")) +
  annotate("segment", x = school_mean, xend = class_mean, y = max_y * .42, yend = max_y * .42, arrow = arr, color = ORANGE, linewidth = .8) +
  annotate("label", x = school_mean + 1, y = max_y * .42, hjust = 0, label = 'u[jk] == "\u22125"', parse = TRUE,
           size = 4, fill = "white", border.colour = NA, label.padding = unit(1, "pt")) +
  annotate("segment", x = class_mean, xend = student_score, y = max_y * .16, yend = max_y * .16, arrow = arr, color = RED, linewidth = .8) +
  annotate("label", x = (class_mean + student_score) / 2, y = max_y * .185, vjust = 0, label = 'e[ijk] == "+3"', parse = TRUE,
           size = 4, fill = "white", border.colour = NA, label.padding = unit(1, "pt")) +
  annotate("point", x = student_score, y = 0, color = RED, size = 3.2) +
  annotate("text", x = student_score - .5, y = -max_y * .075, label = "Y[ijk] == 108", parse = TRUE, hjust = 0, color = RED,
           size = 4) +
  annotate("label", x = delta_0 - .8, y = max_y * 1.12, label = "delta[0] == 100", parse = TRUE, hjust = 1, color = GRN,
           size = 4, fill = "white", border.colour = NA, label.padding = unit(1, "pt")) +
  annotate("label", x = school_mean + .8, y = max_y * 1.12, label = "gamma[0*k] == 110", parse = TRUE, hjust = 0, color = ORANGE,
           size = 4, fill = "white", border.colour = NA, label.padding = unit(1, "pt")) +
  annotate("label", x = class_mean, y = max_y * 1.26, label = "beta[0*jk] == 105", parse = TRUE, color = BLU,
           size = 4, fill = "white", border.colour = NA, label.padding = unit(1, "pt")) +
  annotate("text", x = 51, y = max_y * .97, hjust = 0, color = GRN, size = 3.9, parse = TRUE,
           label = '"school means:"~~N(delta[0], psi^2)') +
  annotate("text", x = 51, y = max_y * .85, hjust = 0, color = ORANGE, size = 3.9, parse = TRUE,
           label = '"classroom means:"~~N(gamma[0*k], tau^2)') +
  annotate("text", x = 51, y = max_y * .73, hjust = 0, color = BLU, size = 3.9, parse = TRUE,
           label = '"student scores:"~~N(beta[0*jk], sigma^2)') +
  scale_color_manual(values = LCOL, guide = "none") + scale_fill_manual(values = LCOL, guide = "none") +
  scale_linetype_manual(values = LTY, guide = "none") +
  scale_x_continuous(breaks = seq(50, 150, 25)) +
  coord_cartesian(xlim = plot_range, ylim = c(0, max_y * 1.33), expand = FALSE, clip = "off") +
  scale_y_continuous(breaks = NULL) +
  labs(title = "One student's score: grand mean plus three deviations",
       subtitle = expression({Y[ijk] == delta[0] + r[k] + u[jk] + e[ijk]}*":   "*{108 == 100 + 10 - 5 + 3}),
       x = "score (illustrative units)", y = NULL) +
  theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line.x = element_line(color = INK),
        plot.margin = margin(6, 14, 6, 14))
save(p4, "fig04_dgp.png", 6.1, 4.1)

## fig05 --- 98 classroom means, grouped by school ----------------------------------------------
cm <- relyea |> group_by(block_id, cluster_id) |> summarise(m = mean(bgk_score), n = n(), .groups = "drop")
sm <- relyea |> group_by(block_id) |> summarise(sm = mean(bgk_score), .groups = "drop") |> arrange(sm) |>
  mutate(sch = factor(block_id, levels = block_id), xi = row_number())
cm <- cm |> left_join(sm, by = "block_id")
within_sd <- cm |> group_by(block_id) |> summarise(s = sd(m)) |> pull(s)
cat(sprintf("fig05: school means %.2f to %.2f (SD %.2f); classroom means %.2f to %.2f; median within-school SD of classroom means %.2f\n",
            min(sm$sm), max(sm$sm), sd(sm$sm), min(cm$m), max(cm$m), median(within_sd)))
p5 <- ggplot(cm, aes(xi, m)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = MUT) +
  geom_segment(data = sm, aes(x = xi - .38, xend = xi + .38, y = sm, yend = sm), color = INK, linewidth = 1) +
  geom_point(color = BLU, size = 1.9, alpha = .8) +
  scale_x_continuous(breaks = sm$xi, labels = sm$block_id, expand = expansion(add = .6)) +
  labs(x = "school (ordered by its mean BGK score)", y = "mean BGK score (SD units)",
       title = "Classroom means grouped by school",
       subtitle = "Dot: a classroom's mean; black bar: the school's mean") +
  theme(panel.grid.major.x = element_blank(), axis.text.x = element_text(size = 10.5))
save(p5, "fig05_classrooms.png", 6.6, 3.6)

## fig06 --- the three kinds of student pairs (schematic) --------------------------------------
icc_cl <- (psi2 + tau2) / tot; icc_sc <- psi2 / tot
rooms <- tibble(sch = c("A", "A", "B", "B"), room = c(1, 2, 1, 2),
                xmin = c(.35, 2.45, 5.35, 7.45)) |> mutate(xmax = xmin + 1.8)
dots <- rooms |> rowwise() |> reframe(sch, room, xmin, k = 1:6) |>
  mutate(x = xmin + .3 + ((k - 1) %% 3) * .6, y = .85 + ((k - 1) %/% 3) * .45)
pick <- function(s, r, k) dots |> filter(sch == s, room == r, k == !!k) |> select(x, y)
P1a <- pick("A", 1, 1); P1b <- pick("A", 1, 3)        # same classroom
P2a <- pick("A", 1, 6); P2b <- pick("A", 2, 4)        # same school, different classroom
P3a <- pick("A", 2, 6); P3b <- pick("B", 1, 5)        # different schools
seg6 <- function(P, Q, pr) tibble(pr = pr, x = P$x, y = P$y, xend = Q$x, yend = Q$y)
PCOL <- c("1" = RED, "2" = BLU, "3" = GRN)
hl <- bind_rows(P1a, P1b, P2a, P2b, P3a, P3b) |> mutate(pr = rep(c("1", "2", "3"), each = 2))
lab_xy <- tibble(pr = c("1", "2", "3"), x = c(1.25, 2.3, 4.8), y = c(.55, 1.75, 2.0))
mark <- function(d) list(geom_point(data = d, aes(x, y, color = pr), shape = 21, size = 7, fill = "white", stroke = 1.1),
                         geom_text(data = d, aes(x, y, label = pr, color = pr), size = 3.9, fontface = "bold"))
p6a <- ggplot() +
  annotate("rect", xmin = .15, xmax = 4.45, ymin = .05, ymax = 1.95, fill = "grey93", color = "grey70") +
  annotate("rect", xmin = 5.15, xmax = 9.45, ymin = .05, ymax = 1.95, fill = "grey93", color = "grey70") +
  annotate("text", x = c(2.3, 7.3), y = 2.18, label = c("school A", "school B"), size = 4, fontface = "bold", color = INK) +
  geom_rect(data = rooms, aes(xmin = xmin, xmax = xmax, ymin = .25, ymax = 1.55), fill = "white", color = MUT) +
  geom_text(data = rooms, aes(x = (xmin + xmax) / 2, y = 1.75, label = paste("classroom", room)), size = 3.7, color = MUT) +
  geom_point(data = dots, aes(x, y), color = "grey55", size = 2.4) +
  geom_curve(data = seg6(P1a, P1b, "1"), aes(x = x, y = y, xend = xend, yend = yend), curvature = .5, color = RED, linewidth = .9) +
  geom_curve(data = seg6(P2a, P2b, "2"), aes(x = x, y = y, xend = xend, yend = yend), curvature = -.45, color = BLU, linewidth = .9) +
  geom_curve(data = seg6(P3a, P3b, "3"), aes(x = x, y = y, xend = xend, yend = yend), curvature = -.35, color = GRN, linewidth = .9) +
  geom_point(data = hl, aes(x, y, color = pr), size = 3.2) +
  mark(lab_xy) +
  scale_color_manual(values = PCOL, guide = "none") +
  coord_cartesian(xlim = c(0, 9.6), ylim = c(0, 2.32), expand = FALSE, clip = "off") +
  theme_void(base_size = 14)
keytxt <- tibble(y = c(3, 2, 1), x = 0, pr = c("1", "2", "3"),
  lab = c(sprintf('"same classroom (so same school):  "~~frac(hat(psi)^2 + hat(tau)^2, hat(psi)^2 + hat(tau)^2 + hat(sigma)^2) == %s', f2(icc_cl, 3)),
          sprintf('"same school, different classroom:  "~~frac(hat(psi)^2, hat(psi)^2 + hat(tau)^2 + hat(sigma)^2) == %s', f2(icc_sc, 3)),
          '"different schools:  "~~0'))
p6b <- ggplot(keytxt) +
  mark(keytxt) +
  geom_text(aes(x = .25, y = y, label = lab), parse = TRUE, hjust = 0, size = 3.9, color = INK) +
  scale_color_manual(values = PCOL, guide = "none") +
  coord_cartesian(xlim = c(-.15, 6), ylim = c(.5, 3.5), expand = FALSE, clip = "off") +
  theme_void(base_size = 14)
p6 <- (p6a / p6b + plot_layout(heights = c(1.15, 1))) +
  plot_annotation(title = "Correlation between two students' BGK scores under m1",
                  theme = theme(plot.title = element_text(face = "bold", size = 14, color = INK),
                                plot.margin = margin(6, 6, 6, 6)))
save(p6, "fig06_pairs.png", 6.6, 4.4)

## fig07 --- school-specific effects: raw differences and m3's estimates (lecture PDF 26, corrected) ----
eb <- tidy(m3, effects = "ran_vals", conf.int = TRUE, conf.level = .95) |>
  filter(group == "block_id", term == "treat1") |>
  mutate(block_id = as.numeric(level), eb = estimate + ate, eb_lo = conf.low + ate, eb_hi = conf.high + ate)
ebi <- tidy(m3, effects = "ran_vals") |> filter(group == "block_id") |> select(level, term, estimate) |>
  pivot_wider(names_from = term, values_from = estimate)
cat(sprintf("fig07: slope BLUP / intercept BLUP ratio %.4f to %.4f\n",
            min(ebi$treat1 / ebi$`(Intercept)`), max(ebi$treat1 / ebi$`(Intercept)`)))
cat(sprintf("fig07: EB effects %.4f to %.4f; raw differences %.3f to %.3f; cor(raw, EB) %.2f\n",
            min(eb$eb), max(eb$eb), min(agg$estimate), max(agg$estimate),
            cor(agg$estimate, eb$eb[match(agg$block_id, eb$block_id)])))
ord7 <- eb |> arrange(desc(eb)) |> pull(block_id)                      # lowest EB effect at the top
f7 <- agg |> mutate(z_lo = estimate - 1.96 * std.error, z_hi = estimate + 1.96 * std.error) |> select(block_id, nobs, estimate, z_lo, z_hi) |> left_join(eb |> select(block_id, eb, eb_lo, eb_hi), by = "block_id") |>
  mutate(sch = factor(block_id, levels = ord7))
print(f7 |> arrange(eb) |> mutate(across(where(is.numeric), ~ round(., 3))), n = 26)
refl <- list(geom_vline(xintercept = 0, linetype = "dotted", color = INK),
             geom_vline(xintercept = ate, linetype = "dashed", color = RED),
             scale_y_discrete(expand = expansion(add = c(.7, 2.6))))
glab <- function(off) annotate("label", x = ate + off, y = 27.5, hjust = 0, label = sprintf("hat(delta)[1] == %s", f2(ate, 3)),
                               parse = TRUE, color = RED, size = 3.9, fill = "white", border.colour = NA, label.padding = unit(1, "pt"))
SZ7 <- list(scale_size_area(max_size = 3.2, breaks = c(30, 75, 150), limits = c(0, 150), name = "students in the school"),
            guides(size = guide_legend(override.aes = list(color = "grey45"))))
p7a <- ggplot(f7, aes(y = sch)) + refl +
  geom_errorbar(aes(xmin = z_lo, xmax = z_hi), width = 0, orientation = "y", color = MUT) +
  geom_point(aes(x = estimate, size = nobs), color = INK) + glab(.04) + SZ7 +
  labs(x = "adaptive mean minus core mean\n(BGK SD units)", y = "school", title = "Raw difference",
       subtitle = "estimate ± 1.96 SE")
p7b <- ggplot(f7, aes(y = sch)) + refl +
  geom_errorbar(aes(xmin = eb_lo, xmax = eb_hi), width = 0, orientation = "y", color = BLU) +
  geom_point(aes(x = eb, size = nobs), color = BLU) + glab(.004) + SZ7 +
  scale_x_continuous(breaks = c(0, .1, .2)) +
  labs(x = expression(atop(hat(delta)[1] + hat(r)[1*k], "(BGK SD units)")), y = NULL,
       title = "m3's estimate", subtitle = "95% interval") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
save(p7a + p7b + plot_layout(widths = c(1.35, 1), guides = "collect") & theme(legend.position = "bottom"), "fig07_school_effects.png", 6.8, 6.2)

## fig08 --- Cohen's d for one school (school 24) ----------------------------------------------
s24 <- agg |> filter(block_id == 24)
cat(sprintf("fig08: school 24 n_T %d n_C %d; mean_T %.4f mean_C %.4f diff %.4f pooled SD %.4f d %.4f\n",
            s24$n_T, s24$n_C, s24$mean_T, s24$mean_C, s24$estimate, s24$sigma, s24$d))
mT <- s24$mean_T; mC <- s24$mean_C; sp <- s24$sigma
xx8 <- seq(-3, 3.2, length.out = 500)
cur8 <- bind_rows(tibble(x = xx8, y = dnorm(xx8, mT, sp), arm = "adaptive (treat = 1)"),
                  tibble(x = xx8, y = dnorm(xx8, mC, sp), arm = "core (treat = 0)")) |>
  mutate(arm = factor(arm, levels = names(ARM)))
pk <- dnorm(0, 0, sp)
set.seed(924)
st24 <- relyea |> filter(block_id == 24) |> mutate(arm = arm_lab(treat),
  y = ifelse(treat == "1", -.06, -.13) * pk / .41 + runif(n(), -.018, .018))
p8 <- ggplot(cur8, aes(x, y)) +
  geom_area(aes(fill = arm), alpha = .15, position = "identity") +
  geom_line(aes(color = arm), linewidth = 1) +
  geom_segment(data = tibble(x = c(mT, mC), arm = factor(names(ARM), levels = names(ARM))),
               aes(x = x, xend = x, y = 0, yend = pk * 1.0, color = arm), linewidth = .8, linetype = "dashed") +
  geom_point(data = st24, aes(bgk_score, y, color = arm), size = 1.3, alpha = .85) +
  annotate("segment", x = mC, xend = mT, y = pk * 1.10, yend = pk * 1.10, color = INK, linewidth = .7,
           arrow = arrow(length = unit(.14, "cm"), ends = "both", type = "closed")) +
  annotate("text", x = mC - .08, y = pk * 1.10, hjust = 1, size = 4, color = INK,
           label = sprintf("difference in means = %s", f2(s24$estimate, 3))) +
  annotate("segment", x = mC, xend = mC + sp, y = pk * 1.25, yend = pk * 1.25, color = INK, linewidth = .7,
           arrow = arrow(length = unit(.14, "cm"), ends = "both", type = "closed")) +
  annotate("text", x = mC - .08, y = pk * 1.25, hjust = 1, size = 4, color = INK,
           label = sprintf("1 pooled SD = %s", f2(sp, 3))) +
  annotate("text", x = 1.33, y = pk * .86, hjust = 0, size = 4, color = INK, parse = TRUE,
           label = sprintf('{d == frac(%s - %s, %s)} == %s', f2(mT, 3), f2(mC, 3), f2(sp, 3), f2(s24$d, 2))) +
  scale_color_manual(values = ARM, name = NULL) + scale_fill_manual(values = ARM, name = NULL) +
  scale_x_continuous(breaks = -3:3,
                     sec.axis = sec_axis(~ (. - mC) / sp, breaks = -3:3,
                                         name = "pooled SDs above or below the core mean")) +
  scale_y_continuous(breaks = NULL, expand = expansion(mult = c(.02, .04))) +
  coord_cartesian(xlim = c(-3, 3.2), clip = "off") +
  labs(x = "BGK score (SD units)", y = NULL, title = "School 24: Cohen's d") +
  theme(panel.grid.major.y = element_blank())
save(p8, "fig08_cohen_d.png", 6.1, 4.5)

## fig09 --- the 26 values of d (lecture PDF 34, drawn with d) ---------------------------------------
dd <- agg |> mutate(bin = round(d / .1) * .1) |> arrange(bin, d) |> group_by(bin) |> mutate(stk = row_number()) |> ungroup()
p9 <- ggplot(dd, aes(bin, stk)) +
  geom_vline(xintercept = mean(agg$d), linetype = "dashed", color = RED) +
  geom_point(size = 5.2, color = BLU) +
  annotate("text", x = mean(agg$d) + .05, y = 6, hjust = 0, size = 3.9, color = RED,
           label = sprintf("mean of the 26 values = %s", f2(mean(agg$d), 2))) +
  scale_y_continuous(breaks = 1:6, limits = c(.5, 6.4)) +
  scale_x_continuous(breaks = seq(-1, .8, .2), limits = c(-1, .8)) +
  labs(x = "Cohen's d (school's own pooled SD units)", y = "number of schools",
       title = "Cohen's d in the 26 schools", subtitle = "Each dot is a school (values rounded to 0.1)") +
  theme(panel.grid.major.y = element_blank())
save(p9, "fig09_d.png", 5.6, 3.3)

## fig10 --- SE(d) against school size (lecture PDF 36, standard formula) ----------------------------
nn <- seq(24, 150, 1)
lab10 <- tibble(nobs = c(31, 60, 149), se_d_std = c(.397, .292, .151), hj = c(0, 0, 1),
                lab = c("schools 8, 19, 22 and 26", "school 1", "school 14"))
p10 <- ggplot(agg, aes(nobs, se_d_std)) +
  geom_line(data = tibble(nobs = nn, se_d_std = 2 / sqrt(nn)), linetype = "dashed", color = MUT) +
  geom_point(color = BLU, size = 2.2) +
  geom_text(data = lab10, aes(label = lab, hjust = hj), size = 3.8, color = INK) +
  annotate("text", x = 112, y = 2 / sqrt(112) + .028, label = "2 / sqrt(n)", parse = TRUE, color = MUT, size = 3.9) +
  scale_y_continuous(limits = c(.14, .42)) + scale_x_continuous(limits = c(20, 155)) +
  labs(x = "school size (number of students)", y = "standard error of d",
       title = "Standard error of d by school size",
       subtitle = "Dashed: equal-sized groups and d = 0")
save(p10, "fig10_se_n.png", 5.6, 3.7)

## fig11 --- common-effect vs random-effects data-generating process (conceptual simulation) ---------
# Illustrative values, not the trial's. Seed 53 is the first seed (searching 1, 2, ...) whose draws are
# typical: every |e_j / SE_j| < 1.6, every |u_j| < 0.36, SD of the u_j between 0.17 and 0.23, 3 to 5 u_j > 0.
set.seed(53)
K11 <- 8; g0 <- .20; tau11 <- .20
se11 <- c(.06, .08, .10, .12, .14, .16, .18, .20)
e11 <- rnorm(K11, 0, se11); u11 <- rnorm(K11, 0, tau11)
sim <- tibble(j = 1:K11, se = se11, e = e11, u = u11) |>
  mutate(true_fe = g0, obs_fe = g0 + e, true_re = g0 + u, obs_re = g0 + u + e)
print(sim |> mutate(across(where(is.numeric), ~ round(., 3))))
YT <- K11 + 1.0                         # baseline of the top curve
UO <- .22                               # u_j segments are drawn just above each school's row
SHP <- c("true effect of the school" = 1, "observed effect" = 16)
r8 <- filter(sim, j == K11)
p11a <- ggplot(sim, aes(y = j)) +
  geom_vline(xintercept = g0, color = RED, linewidth = .6) +
  geom_segment(aes(x = true_fe, xend = obs_fe, yend = j), color = MUT, linewidth = .8) +
  geom_point(aes(x = true_fe, shape = "true effect of the school"), size = 2.6, color = RED, stroke = 1) +
  geom_point(aes(x = obs_fe, shape = "observed effect"), size = 2.6, color = INK) +
  annotate("text", x = (r8$true_fe + r8$obs_fe) / 2, y = K11 - .38, label = "e[j]", parse = TRUE, size = 3.9, color = MUT) +
  annotate("segment", x = g0, xend = g0, y = YT, yend = YT + 1.6, color = RED, linewidth = 1.6) +
  annotate("label", x = g0, y = YT + 2.25, size = 3.9, color = RED, label = '"one true effect: "*gamma[0]', parse = TRUE,
           fill = "white", border.colour = NA, label.padding = unit(1, "pt")) +
  labs(title = "Common-effect model", subtitle = expression(d[j] == gamma[0] + e[j]))
dens <- tibble(x = seq(-.45, .85, length.out = 300)) |> mutate(y = YT + 1.6 * dnorm(x, g0, tau11) / dnorm(0, 0, tau11))
p11b <- ggplot(sim, aes(y = j)) +
  geom_ribbon(data = dens, aes(x = x, ymin = YT, ymax = y), inherit.aes = FALSE, fill = RED, alpha = .12) +
  geom_line(data = dens, aes(x = x, y = y), inherit.aes = FALSE, color = RED, linewidth = .6) +
  geom_vline(xintercept = g0, color = RED, linewidth = .6) +
  geom_segment(aes(x = g0, xend = true_re, y = j + UO, yend = j + UO), color = RED, linewidth = .8) +
  geom_segment(aes(x = true_re, xend = obs_re, yend = j), color = MUT, linewidth = .8) +
  geom_point(aes(x = true_re, shape = "true effect of the school"), size = 2.6, color = RED, stroke = 1) +
  geom_point(aes(x = obs_re, shape = "observed effect"), size = 2.6, color = INK) +
  annotate("text", x = (g0 + r8$true_re) / 2, y = K11 + UO + .3, label = "u[j]", parse = TRUE, size = 3.9, color = RED) +
  annotate("text", x = (r8$true_re + r8$obs_re) / 2, y = K11 - .38, label = "e[j]", parse = TRUE, size = 3.9, color = MUT) +
  annotate("label", x = g0, y = YT + 2.25, size = 3.9, color = RED, parse = TRUE,
           label = '"true effects: "*N(gamma[0], tau^2)', fill = "white", border.colour = NA, label.padding = unit(1, "pt")) +
  labs(title = "Random-effects model", subtitle = expression(d[j] == gamma[0] + u[j] + e[j]))
lims11 <- list(scale_x_continuous(limits = c(-.45, .85), breaks = seq(-.4, .8, .2)),
               scale_y_continuous(breaks = 1:K11, limits = c(.5, YT + 2.6), expand = c(0, 0)),
               scale_shape_manual(values = SHP, name = NULL, breaks = names(SHP)),
               labs(x = "effect size (SD units)", y = "school"),
               theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()))
p11 <- (p11a + lims11) + (p11b + lims11 + labs(y = NULL)) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
save(p11, "fig11_fe_re.png", 6.8, 4.7)

## fig12 --- common-effect and random-effects weights against school size ---------------------------
w12 <- agg |> select(block_id, nobs, w_common, w_random) |>
  pivot_longer(c(w_common, w_random), names_to = "type", values_to = "w") |>
  mutate(type = recode(type, w_common = "common effect: 1 / V", w_random = "random effects: 1 / (V + tau^2)"))
WCOL <- c("common effect: 1 / V" = BLU, "random effects: 1 / (V + tau^2)" = ORANGE)
WLAB <- c(expression("common effect:  "*1/V[j]), expression("random effects:  "*1/(V[j] + hat(tau)^2)))
lab12 <- agg |> filter(block_id %in% c(14, 19)) |>
  mutate(lab = sprintf("school %d: %s%% and %s%%", block_id, f2(w_common), f2(w_random)))
p12 <- ggplot(agg, aes(nobs)) +
  geom_hline(yintercept = 100 / 26, linetype = "dotted", color = INK) +
  annotate("text", x = 151, y = 100 / 26 - .32, hjust = 1, size = 3.7, color = INK,
           label = sprintf("equal weights: 100%% / 26 = %s%%", f2(100 / 26))) +
  geom_segment(aes(xend = nobs, y = w_common, yend = w_random), color = "grey70") +
  geom_point(data = w12, aes(y = w, color = type, shape = type), size = 2.3) +
  geom_text(data = filter(lab12, block_id == 14), aes(x = nobs - 3, y = w_common, label = sub(": ", ":\n", lab)), hjust = 1, size = 3.7, lineheight = .9) +
  geom_text(data = filter(lab12, block_id == 19), aes(x = nobs + 4, y = w_common, label = lab), hjust = 0, size = 3.7) +
  scale_color_manual(values = WCOL, labels = WLAB, name = NULL) +
  scale_shape_manual(values = c(16, 17), labels = WLAB, name = NULL) +
  scale_y_continuous(limits = c(0.5, 8.3), breaks = seq(0, 8, 2), labels = function(x) paste0(x, "%")) +
  scale_x_continuous(limits = c(20, 152)) +
  guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  labs(x = "school size (number of students)", y = "weight in the pooled estimate",
       title = "Each school's weight in the meta-analysis",
       subtitle = "Gray segment joins a school's two weights")
save(p12, "fig12_weights.png", 5.9, 4.2)

## fig13 --- forest plot of the printed meta-analysis (Cohen's d, standard SEs) -------------------
f13 <- agg |> arrange(desc(d)) |> mutate(row = row_number() + 3)       # lowest estimate at the top
pool <- tibble(row = c(2, 1), lab = c("common effect", "random effects"),
               est = c(muf, mur), lo = c(muf - 1.96 * sef, mur - 1.96 * ser), hi = c(muf + 1.96 * sef, mur + 1.96 * ser),
               w = c(NA, 100))
dia <- pool |> rowwise() |> reframe(lab, row, x = c(lo, est, hi, est), y = row + c(0, .32, 0, -.32))
ylabs <- c(setNames(paste("school", f13$block_id), f13$row), setNames(pool$lab, pool$row))
YMAX <- max(f13$row) + .7
p13a <- ggplot() +
  geom_vline(xintercept = 0, color = INK, linewidth = .4) +
  geom_vline(xintercept = mur, linetype = "dashed", color = RED) +
  geom_hline(yintercept = 2.8, color = "grey70", linewidth = .3) +
  geom_errorbar(data = f13, aes(y = row, xmin = ci_lo, xmax = ci_hi), width = 0, orientation = "y", color = MUT) +
  geom_point(data = f13, aes(y = row, x = d, size = w_random), shape = 15, color = BLU) +
  geom_polygon(data = dia, aes(x, y, group = lab), fill = RED, color = RED) +
  scale_size_area(max_size = 3.6, guide = "none") +
  scale_y_continuous(breaks = as.numeric(names(ylabs)), labels = ylabs, limits = c(.4, YMAX), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(-1.5, 1.5, .5)) +
  labs(x = "Cohen's d, adaptive minus core\n(pooled within-school SD units)", y = NULL) +
  theme(panel.grid.major.y = element_blank(), axis.text.y = element_text(size = 10.5))
txt13 <- bind_rows(f13 |> transmute(row, ci = sprintf("%s [%s, %s]", f2(d, 3), f2(ci_lo, 3), f2(ci_hi, 3)), w = paste0(f2(w_random, 1), "%")),
                   pool |> transmute(row, ci = sprintf("%s [%s, %s]", f2(est, 3), f2(lo, 3), f2(hi, 3)),
                                     w = ifelse(is.na(w), "", "100%")))
txt13 <- txt13 |> mutate(ci = gsub("-", "−", ci))
p13b <- ggplot(txt13, aes(y = row)) +
  geom_text(aes(x = 0, label = ci), hjust = 0, size = 3.7, color = INK) +
  geom_text(aes(x = 1.62, label = w), hjust = 1, size = 3.7, color = INK) +
  annotate("text", x = c(0, 1.62), y = YMAX - .1, label = c("d [95% CI]", "weight"),
           hjust = c(0, 1), vjust = 0, size = 3.7, fontface = "bold", color = INK) +
  scale_y_continuous(limits = c(.4, YMAX), expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 1.65), clip = "off") +
  theme_void(base_size = 14) + theme(plot.margin = margin(14, 4, 4, 2))
p13 <- p13a + p13b + plot_layout(widths = c(1.3, 1)) +
  plot_annotation(title = "Forest plot of the 26 schools' Cohen's d",
                  subtitle = expression("Square size: random-effects weight;  "*{hat(tau)^2 == "0.0290"}*",  "*{I^2 == "41.8%"}),
                  theme = theme(plot.title = element_text(face = "bold", size = 14, color = INK),
                                plot.subtitle = element_text(size = 11.2, color = MUT)))
save(p13, "fig13_forest.png", 6.8, 6.4)

cat("done\n")
