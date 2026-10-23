# Figures for Chapter 3: Within, Between, and Contextual Effects.
# Written by Bill Young. Run from the root of the book repository; set DATA
# to wherever your copy of jebs_full.csv lives. Figures are written to OUT.
# Uses only tidyverse + lme4 + scales, so it runs on a plain R install.
suppressMessages({library(tidyverse); library(lme4); library(scales)})
DATA <- "data/jebs_full.csv"
OUT  <- "chapters/images/03/"
dir.create(OUT, showWarnings = FALSE)

jebs <- read_csv(DATA, show_col_types = FALSE) |>
  mutate(cluster_id = factor(cluster_id))
agg <- jebs |> distinct(cluster_id, .keep_all = TRUE)

INK<-"#1A1D23"; BLU<-"#2C5F8A"; RED<-"#B03A2E"; GRN<-"#2E7D5B"; MUT<-"#5C6672"; RUL<-"#DDE2E6"
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

w    <- lm(demeaned_read ~ demeaned_map, jebs)
b    <- lm(sch_read ~ sch_map, agg)
ri   <- lmer(read_comp ~ std_baseline + (1|cluster_id), jebs)
mund <- lmer(read_comp ~ std_baseline + sch_map + (1|cluster_id), jebs)
rw   <- lmer(read_comp ~ demeaned_map + sch_map + (1|cluster_id), jebs)
m1   <- lmer(read_comp ~ treat*std_baseline + (1|cluster_id), jebs)
mi   <- lmer(read_comp ~ std_baseline*sch_map + (1|cluster_id), jebs)

# fixed-effects prediction with a confidence band, no extra packages
predband <- function(model, newdata) {
  X  <- model.matrix(delete.response(terms(model)), newdata)
  fit <- as.vector(X %*% fixef(model))
  se  <- sqrt(rowSums((X %*% as.matrix(vcov(model))) * X))
  newdata |> mutate(predicted = fit, lo = fit - 1.96*se, hi = fit + 1.96*se)
}

## fig01 -- the two regressions, side by side
d1 <- bind_rows(
  agg  |> transmute(x = sch_map,      y = sch_read,      panel = "Between: 110 school averages, slope 7.96"),
  jebs |> transmute(x = demeaned_map, y = demeaned_read, panel = "Within: 7,797 student deviations, slope 6.96"))
p1 <- ggplot(d1, aes(x, y)) +
  geom_point(data = ~filter(.x, grepl("^Between", panel)), color = RED, size = 1.6, alpha = .8) +
  geom_point(data = ~filter(.x, grepl("^Within",  panel)), color = INK, size = .5, alpha = .10) +
  geom_smooth(method = "lm", se = FALSE, color = BLU, linewidth = .9, formula = y ~ x) +
  facet_wrap(~panel, scales = "free") +
  labs(x = "BOY MAP score: school mean (left) or deviation from it (right)",
       y = "EOY reading score:\nschool mean (left), deviation (right)",
       title = "The same relationship, computed two ways")
ggsave(paste0(OUT,"fig01_between_within.png"), p1, width = 8.83, height = 3.77, dpi = 200, bg = "white")

## fig02 -- the four slopes on one set of axes (lecture now fits a total regression too)
tot  <- lm(read_comp ~ std_baseline, jebs)
gm_x <- mean(jebs$std_baseline); gm_y <- mean(jebs$read_comp)
sl <- tibble(lab = factor(c("Within (6.96)","Between (7.96)","Random intercepts (6.99)","Total (7.46)"),
                          levels = c("Within (6.96)","Random intercepts (6.99)","Total (7.46)","Between (7.96)")),
             slope = c(coef(w)[2], coef(b)[2], fixef(ri)[2], coef(tot)[2])) |>
      mutate(int = gm_y - slope*gm_x)
p2 <- ggplot(jebs, aes(std_baseline, read_comp)) +
  geom_point(color = MUT, alpha = .06, size = .45) +
  geom_abline(data = sl, aes(intercept = int, slope = slope, color = lab, linetype = lab), linewidth = 1) +
  scale_color_manual(values = c(GRN, INK, BLU, RED), name = NULL) +
  scale_linetype_manual(values = c("solid","longdash","dotted","solid"), name = NULL) +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  coord_cartesian(xlim = c(-3,3), ylim = c(25,75)) +
  labs(x = "BOY MAP score (SD units)", y = "EOY reading score",
       title = "Four slopes for the same data",
       subtitle = "All four are drawn through the overall mean so the slopes can be compared.\nRandom intercepts sits almost on top of within; total sits between within and between.")
ggsave(paste0(OUT,"fig02_three_slopes.png"), p2, width = 6.81, height = 4.35, dpi = 200, bg = "white")

## fig03 -- precision weights
prec <- tibble(slope = c(coef(w)[2], coef(b)[2]),
               se    = c(summary(w)$coef[2,2], summary(b)$coef[2,2]),
               lab   = c("Within","Between")) |>
        mutate(wt = (1/se^2)/sum(1/se^2), lab = factor(lab, levels = c("Within","Between")))
p3 <- ggplot(prec, aes(lab, wt, fill = lab)) +
  geom_col(width = .55) +
  geom_text(aes(label = percent(wt, accuracy = 0.1)), vjust = -0.5, size = 4.2, color = INK) +
  scale_fill_manual(values = c(GRN, RED), guide = "none") +
  scale_y_continuous(limits = c(0,1.05), labels = percent) +
  labs(x = NULL, y = "Share of the weight",
       title = "Why the blend is almost the within slope",
       subtitle = "Weight is proportional to 1 over the squared SE.\nWithin SE 0.079; between SE 0.522.")
ggsave(paste0(OUT,"fig03_precision_weights.png"), p3, width = 6.07, height = 3.58, dpi = 200, bg = "white")

## fig04 -- the Mundlak picture
jebs$mund_fit <- fitted(mund)
p4 <- ggplot(jebs, aes(std_baseline, mund_fit, group = cluster_id)) +
  geom_line(color = MUT, alpha = .30, linewidth = .3) +
  geom_point(data = agg, aes(sch_map, sch_read), inherit.aes = FALSE, color = RED, size = 1.5, alpha = .85) +
  geom_abline(intercept = coef(b)[1], slope = coef(b)[2], color = RED, linewidth = .9, linetype = "longdash") +
  geom_abline(intercept = fixef(rw)[1], slope = fixef(rw)[3], color = BLU, linewidth = .9, linetype = "dotted") +
  coord_cartesian(xlim = c(-3,3), ylim = c(25,75)) +
  labs(x = "BOY MAP score (SD units)", y = "EOY reading score",
       title = "Three slopes at once",
       subtitle = "Gray: the 110 school lines, slope 6.96. Red dots: observed school means.\nDashed red: mean-on-means, 7.96. Dotted blue: MLM between, 8.72.")
ggsave(paste0(OUT,"fig04_mundlak.png"), p4, width = 6.99, height = 4.23, dpi = 200, bg = "white")

## fig05 -- treatment, two views
d5 <- bind_rows(mutate(jebs, panel = "Every student, jittered"),
                mutate(jebs, panel = "The same data as a boxplot"))
p5 <- ggplot(d5, aes(factor(treat, 0:1, c("Control","Treated")), read_comp)) +
  geom_jitter(data = ~filter(.x, panel == "Every student, jittered"),
              width = .10, alpha = .08, size = .45, color = INK) +
  geom_boxplot(data = ~filter(.x, panel == "The same data as a boxplot"),
               fill = "#E8EEF4", color = INK, width = .5, outlier.size = .4) +
  facet_wrap(~panel) +
  labs(x = "Treatment status", y = "EOY reading score",
       title = "EOY reading scores by treatment status")
ggsave(paste0(OUT,"fig05_treatment.png"), p5, width = 8.28, height = 3.41, dpi = 200, bg = "white")

## fig06 -- the treatment estimate four ways
se_df <- tibble(model = factor(c("OLS","OLS + MAP","MLM","MLM + MAP"),
                               levels = c("OLS","OLS + MAP","MLM","MLM + MAP")),
                est = c(coef(lm(read_comp~treat,jebs))[2],
                        coef(lm(read_comp~treat+std_baseline,jebs))[2],
                        fixef(lmer(read_comp~treat+(1|cluster_id),jebs))[2],
                        fixef(lmer(read_comp~treat+std_baseline+(1|cluster_id),jebs))[2]),
                se  = c(summary(lm(read_comp~treat,jebs))$coef[2,2],
                        summary(lm(read_comp~treat+std_baseline,jebs))$coef[2,2],
                        sqrt(diag(vcov(lmer(read_comp~treat+(1|cluster_id),jebs))))[2],
                        sqrt(diag(vcov(lmer(read_comp~treat+std_baseline+(1|cluster_id),jebs))))[2]))
p6 <- ggplot(se_df, aes(model, est)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = MUT) +
  geom_errorbar(aes(ymin = est-1.96*se, ymax = est+1.96*se), width = .14, linewidth = .6, color = BLU) +
  geom_point(size = 2.8, color = RED) +
  geom_text(aes(y = est + 1.96*se + 0.22, label = sprintf("SE %.3f", se)), size = 3.5, color = MUT) +
  scale_y_continuous(limits = c(-1.6, 3.1)) +
  labs(x = NULL, y = "Estimated treatment effect (EOY points)",
       title = "The treatment estimate and its 95% interval, four ways",
       subtitle = "OLS treats 7,797 students as independent; the MLM does not.")
ggsave(paste0(OUT,"fig06_treatment_se.png"), p6, width = 6.81, height = 3.86, dpi = 200, bg = "white")

## fig07 -- cross-level interaction
g7 <- expand_grid(std_baseline = seq(-3,3,.25), treat = c(0,1)) |> predband(model = m1) |>
      mutate(grp = factor(treat, 0:1, c("Control","Treated")))
p7 <- ggplot(g7, aes(std_baseline, predicted, color = grp, fill = grp)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .18, color = NA) +
  geom_line(linewidth = .95) +
  scale_color_manual(values = c(RED, BLU), name = NULL) +
  scale_fill_manual(values = c(RED, BLU), name = NULL) +
  labs(x = "BOY MAP score (SD units)", y = "Predicted EOY reading score",
       title = "Cross-level interaction: treatment by BOY MAP",
       subtitle = "Control slope 6.95, treated slope 7.03. The gap between the lines barely widens.")
ggsave(paste0(OUT,"fig07_crosslevel.png"), p7, width = 6.81, height = 3.95, dpi = 200, bg = "white")

## fig08 -- interacted Mundlak
g8 <- expand_grid(std_baseline = seq(-3,3,.25), sch_map = c(-1,0,1)) |> predband(model = mi) |>
      mutate(grp = factor(sch_map, c(-1,0,1), c("1 unit below average","average","1 unit above average")))
p8 <- ggplot(g8, aes(std_baseline, predicted, color = grp)) +
  geom_line(linewidth = .95) +
  scale_color_manual(values = c(RED, BLU, GRN), name = "School mean BOY MAP") +
  guides(color = guide_legend(nrow = 2)) +
  labs(x = "A student's own BOY MAP score (SD units)", y = "Predicted EOY reading score",
       title = "The within-school slope depends on the school",
       subtitle = "Slope 6.21 one unit below average on mean BOY MAP; 7.71 one unit above.")
ggsave(paste0(OUT,"fig08_interacted_mundlak.png"), p8, width = 6.81, height = 3.95, dpi = 200, bg = "white")

## fig09 -- splitting one student's MAP score into a between part and a within part
s13 <- jebs |> filter(cluster_id == 13)
stu <- s13 |> filter(id == 422)
mn  <- unique(s13$sch_map)
p9 <- ggplot(s13, aes(std_baseline, y = 0)) +
  geom_vline(xintercept = 0, color = MUT, linetype = "dashed", linewidth = .45) +
  geom_point(color = MUT, alpha = .35, size = 1.7) +
  geom_vline(xintercept = mn, color = RED, linewidth = .8) +
  geom_point(data = stu, aes(std_baseline, 0), color = BLU, size = 3.4) +
  annotate("segment", x = 0, xend = mn, y = .20, yend = .20, color = RED, linewidth = .7,
           arrow = arrow(length = unit(.16,"cm"), ends = "both", type = "closed")) +
  annotate("text", x = mn + 0.07, y = .32, hjust = 0, label = "between part: 0.459\n(where school 13 sits)",
           color = RED, size = 3.9, lineheight = .95) +
  annotate("segment", x = mn, xend = stu$std_baseline, y = -.20, yend = -.20, color = BLU, linewidth = .7,
           arrow = arrow(length = unit(.16,"cm"), ends = "both", type = "closed")) +
  annotate("text", x = mn + 0.07, y = -.30, hjust = 0,
           label = "within part: 0.958\n(how far above her classmates)", color = BLU, size = 3.9, lineheight = .95) +
  annotate("text", x = -0.06, y = .46, hjust = 1, label = "overall mean, 0", color = MUT, size = 3.7) +
  annotate("text", x = mn + 0.08, y = -.46, hjust = 0, label = "school 13's mean", color = RED, size = 3.7) +
  annotate("text", x = stu$std_baseline + 0.08, y = .46, hjust = 0, label = "student 422, own MAP 1.417", color = BLU, size = 3.7) +
  scale_y_continuous(limits = c(-.56,.56), breaks = NULL) +
  coord_cartesian(xlim = c(-2.7, 3.1)) +
  labs(x = "BOY MAP score (SD units)", y = NULL,
       title = "Student 422's BOY MAP score, split into two parts",
       subtitle = "Each gray dot is one of the 104 students in school 13, placed by their own BOY MAP score.") +
  theme(panel.grid.major.y = element_blank())
ggsave(paste0(OUT,"fig09_decomposition.png"), p9, width = 7.54, height = 3.31, dpi = 200, bg = "white")

cat("figures written\n")
cat(sprintf("precision weights: within %.4f between %.4f\n", prec$wt[1], prec$wt[2]))


## ---------------------------------------------------------------------------
## fig10 / fig11 / fig12 -- the Raudenbush & Bryk pictures (pp. 137 and 140), redrawn.
## fig10 is a schematic with exaggerated slopes; fig11 and fig12 use this lecture's numbers.
## ---------------------------------------------------------------------------
bW <- 6.9569; bB <- 8.7164; g0 <- 49.7406; bC <- bB - bW
cloud <- function(cx, cy, a, b, slope, n = 200) {      # a student cloud whose long axis follows the within line
  t <- seq(0, 2*pi, length.out = n)
  tibble(x = cx + a*cos(t), y = cy + slope*a*cos(t) + b*sin(t))
}
noaxes <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

## fig10 -- schematic: within lines, between line, and the single-slope compromise (R&B p. 137)
sW <- 4.5; sB <- 9; sM <- 5.6
mk <- c(-1, 0, 1); ck <- 50 + sB*mk
cl10 <- bind_rows(lapply(1:3, function(k) cloud(mk[k], ck[k], a = 1.0, b = 3.0, slope = sW) |> mutate(g = k)))
wl10 <- bind_rows(lapply(1:3, function(k) tibble(x = mk[k] + c(-1.15, 1.15), y = ck[k] + sW*c(-1.15, 1.15), g = k)))
lines10 <- bind_rows(
  wl10 |> mutate(what = "Within-school line (one per school, same slope)"),
  tibble(x = c(-2.2, 2.2), y = 50 + sB*c(-2.2, 2.2), g = 9, what = "Between-school line (through the school averages)"),
  tibble(x = c(-2.2, 2.2), y = 50 + sM*c(-2.2, 2.2), g = 8, what = "One slope forced to serve both comparisons"))
p10 <- ggplot() +
  geom_path(data = cl10, aes(x, y, group = g), color = MUT, linewidth = .5) +
  geom_line(data = lines10, aes(x, y, group = g, color = what, linetype = what), linewidth = 1.0) +
  geom_point(data = tibble(x = mk, y = ck), aes(x, y), color = RED, size = 3) +
  annotate("text", x = mk + 1.28, y = ck + sW*1.28 - 1.6, hjust = 0, label = c("School A", "School B", "School C"), color = MUT, size = 3.5) +
  scale_color_manual(values = c("Within-school line (one per school, same slope)" = GRN,
                                "Between-school line (through the school averages)" = RED,
                                "One slope forced to serve both comparisons" = INK), name = NULL) +
  scale_linetype_manual(values = c("Within-school line (one per school, same slope)" = "solid",
                                   "Between-school line (through the school averages)" = "solid",
                                   "One slope forced to serve both comparisons" = "dashed"), name = NULL) +
  guides(color = guide_legend(ncol = 1)) +
  coord_cartesian(xlim = c(-2.3, 2.9), ylim = c(28, 72)) +
  labs(x = "BOY MAP score", y = "EOY reading score",
       title = "Three schools, two kinds of comparison",
       subtitle = "Red dots: each school's average student. Slopes exaggerated to make the picture clear.") +
  noaxes + theme(axis.text = element_blank(), axis.ticks = element_blank(), legend.text = element_text(size = 10))
ggsave(paste0(OUT,"fig10_rb137.png"), p10, width = 6.9, height = 5.0, dpi = 200, bg = "white")

## fig11 -- beta_B split into beta_W + beta_C over one unit of school mean (R&B p. 140), our numbers
cA <- c(0, g0); cB <- c(1, g0 + bB)
yA1 <- cA[2] + bW                                     # school B's within line, one unit to the right of its average
wl11 <- bind_rows(tibble(x = cA[1] + c(-1.1, 1.3), y = cA[2] + bW*c(-1.1, 1.3), g = 1),
                  tibble(x = cB[1] + c(-1.1, 1.3), y = cB[2] + bW*c(-1.1, 1.3), g = 2))
p11 <- ggplot() +
  geom_line(data = wl11, aes(x, y, group = g), color = GRN, linewidth = 1.0) +
  geom_abline(intercept = g0, slope = bB, color = RED, linewidth = 1.0) +
  geom_segment(aes(x = 0, xend = 2.05, y = cA[2], yend = cA[2]), color = MUT, linetype = "dotted") +
  geom_segment(aes(x = 1, xend = 2.05, y = cB[2], yend = cB[2]), color = MUT, linetype = "dotted") +
  geom_segment(aes(x = 1, xend = 1, y = cA[2], yend = yA1), color = GRN, linewidth = 1.1) +
  geom_segment(aes(x = 1, xend = 1, y = yA1, yend = cB[2]), color = BLU, linewidth = 1.1) +
  annotate("segment", x = 2.05, xend = 2.05, y = cA[2], yend = cB[2], color = RED, linewidth = .8,
           arrow = arrow(length = unit(.14,"cm"), ends = "both", type = "closed")) +
  geom_point(data = tibble(x = c(0, 1), y = c(cA[2], cB[2])), aes(x, y), color = RED, size = 3) +
  geom_point(aes(x = 1, y = yA1), color = BLU, size = 3.2) +
  annotate("text", x = 1.06, y = cA[2] + bW/2, label = sprintf("beta[W] == %.2f", bW), parse = TRUE, hjust = 0, color = GRN, size = 4) +
  annotate("text", x = 1.34, y = yA1 + bC/2 + 2.6, label = sprintf("beta[C] == %.2f", bC), parse = TRUE, hjust = 0, color = BLU, size = 4) +
  annotate("segment", x = 1.02, xend = 1.32, y = yA1 + bC/2, yend = yA1 + bC/2 + 2.3, color = BLU, linewidth = .35) +
  annotate("text", x = 2.12, y = (cA[2] + cB[2])/2, label = sprintf("beta[B] == %.2f", bB), parse = TRUE, hjust = 0, color = RED, size = 4) +
  annotate("text", x = 0, y = cA[2] - 7.0, label = "School B\naverage BOY MAP = 0", color = MUT, size = 3.3, lineheight = .9) +
  annotate("text", x = 0.45, y = cB[2] + 3.0, label = "School C\naverage BOY MAP = 1", color = MUT, size = 3.3, lineheight = .9) +
  annotate("text", x = 1.0, y = cA[2] - 3.6, label = "a School B student\nwhose BOY MAP is 1", color = BLU, size = 3.2, lineheight = .9) +
  annotate("segment", x = 1.0, xend = 1.0, y = cA[2] - 2.3, yend = cA[2] - 0.3, color = BLU, linewidth = .35) +
  coord_cartesian(xlim = c(-1.2, 2.75), ylim = c(41, 65)) +
  labs(x = "BOY MAP score", y = "EOY reading score",
       title = "The contextual effect is the gap between the two slopes",
       subtitle = "One unit along a within line: 6.96.  One unit between school averages: 8.72.  Gap: 1.76.") +
  noaxes + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave(paste0(OUT,"fig11_rb140.png"), p11, width = 6.9, height = 4.8, dpi = 200, bg = "white")

## fig12 -- two real students, 1436 (school 1) and 2202 (school 21), both BOY MAP 0.488
m1 <- unique(jebs$sch_map[jebs$cluster_id == 1]); m21 <- unique(jebs$sch_map[jebs$cluster_id == 21])
x0 <- 0.4882
l1  <- function(x) g0 + bW*x + bC*m1
l21 <- function(x) g0 + bW*x + bC*m21
p12 <- ggplot() +
  geom_function(fun = l1,  color = BLU, linewidth = 1.0, xlim = c(-2, 1.6)) +
  geom_function(fun = l21, color = RED, linewidth = 1.0, xlim = c(-2, 1.6)) +
  geom_vline(xintercept = x0, color = MUT, linetype = "dotted") +
  geom_segment(aes(x = x0, xend = x0, y = l1(x0), yend = l21(x0)), color = INK, linewidth = .8,
               arrow = arrow(length = unit(.15,"cm"), ends = "both", type = "closed")) +
  geom_point(aes(x = x0, y = l1(x0)),  color = BLU, size = 3.4) +
  geom_point(aes(x = x0, y = l21(x0)), color = RED, size = 3.4) +
  annotate("text", x = x0 + 0.32, y = (l1(x0)+l21(x0))/2, label = sprintf("%.2f points", l21(x0)-l1(x0)), hjust = 0, size = 3.9, color = INK) +
  annotate("text", x = 1.68, y = l1(1.6),  label = "school 1 (mean BOY MAP -0.41)",  hjust = 0, color = BLU, size = 3.5) +
  annotate("text", x = 1.68, y = l21(1.6) + 0.3, label = "school 21 (mean BOY MAP +0.59)", hjust = 0, color = RED, size = 3.5) +
  annotate("text", x = x0 + 0.32, y = l1(x0) - 2.3,  label = "student 1436, school 1",  hjust = 0, color = BLU, size = 3.6) +
  annotate("text", x = x0 - 0.12, y = l21(x0) + 2.0, label = "student 2202, school 21", hjust = 1, color = RED, size = 3.6) +
  scale_x_continuous(breaks = c(-2, -1, 0, x0, 1, 2), labels = c("-2", "-1", "0", "0.488", "1", "2")) +
  coord_cartesian(xlim = c(-2.05, 3.6)) +
  labs(x = "A student's own BOY MAP score", y = "Predicted EOY reading score",
       title = "Same BOY MAP score, different school",
       subtitle = "Both scored 0.488 at BOY; their schools' averages differ by one unit.") +
  noaxes
ggsave(paste0(OUT,"fig12_two_students.png"), p12, width = 6.9, height = 4.3, dpi = 200, bg = "white")
cat(sprintf("fig12 check: pred1 %.3f pred21 %.3f diff %.3f\n", l1(x0), l21(x0), l21(x0)-l1(x0)))
