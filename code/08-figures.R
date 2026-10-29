# Figures for Chapter 8: Cross-Classified Models, Process Data, and Double Centering.
# Written by Bill Young. Run from the root of the book repository; set DATA
# to wherever your copy of vocab.csv lives. Figures are written to OUT.
suppressMessages({library(tidyverse); library(lme4); library(broom.mixed); library(patchwork)})
DATA <- "data/vocab.csv"
OUT  <- "chapters/images/08/"; dir.create(OUT, showWarnings = FALSE)

# The data, read exactly as the lecture reads it
vocab <- read_csv(DATA, show_col_types = FALSE) |>
  select(id, item, resp, rt, cov_lep, cov_ses_high, cov_ses_med) |>
  mutate(ln_rt = log(rt),
         c_ln_rt = ln_rt - mean(ln_rt)) |>
  group_by(id) |>
  mutate(p_ln_rt = mean(c_ln_rt)) |>
  group_by(item) |>
  mutate(i_ln_rt = mean(c_ln_rt)) |>
  ungroup() |>
  mutate(d_rt = c_ln_rt - p_ln_rt - i_ln_rt,
         item_num = parse_number(item))

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
pct <- function(x, d = 1) paste0(formatC(100 * x, format = "f", digits = d), "%")
FAC <- c(fig01 = 0.92, fig02 = 0.92, fig03 = 0.84, fig04 = 0.87, fig05 = 0.87, fig06 = 0.92, fig07 = 0.92, fig08 = 0.92, fig09 = 0.92,
         fig10 = 0.92, fig11 = 0.92, fig12 = 0.92, fig13 = 0.92, fig14 = 0.92, fig15 = 0.80, fig16 = 0.92, fig17 = 0.71, fig18 = 0.92)
save <- function(p, f, w, h) {
  k <- names(FAC)[startsWith(f, names(FAC))]; sc <- if (length(k)) FAC[[k[length(k)]]] else 1
  cat("###", f, round(sc, 3), "\n")
  ggsave(paste0(OUT, f), p, width = w * sc, height = h * sc, dpi = 200, bg = "white") }

# ---- the lecture's four models ----------------------------------------------------------
m1 <- glmer(resp ~ (1|id) + (1|item), vocab, family = binomial)
m2 <- glmer(resp ~ cov_lep + cov_ses_med + cov_ses_high + (1|id) + (1|item), vocab, family = binomial)
m3 <- glmer(resp ~ c_ln_rt + (1|id) + (1|item), vocab, family = binomial)
m4 <- glmer(resp ~ d_rt + p_ln_rt + i_ln_rt + (1|id) + (1|item), vocab, family = binomial)
g0 <- fixef(m1)[[1]]
tj <- attr(VarCorr(m1)$id, "stddev")[[1]]; ti <- attr(VarCorr(m1)$item, "stddev")[[1]]
cat("m1: g0", g0, " tau_j", tj, " tau_i", ti, "\n")
re1 <- tidy(m1, effects = "ran_vals")          # EB estimates with correct SEs (not augment())
pe <- re1 |> filter(group == "id"); ie <- re1 |> filter(group == "item")

## fig01 --- nested vs crossed, schematic (lecture PDF 10) ---------------------------------
sch <- tibble(s = c("School A", "School B", "School C"), x = c(2, 5, 8), y = 2)
stu <- tibble(x = 1:9, y = 1, s = rep(sch$s, each = 3))
seg <- stu |> left_join(sch |> rename(sx = x, sy = y), by = "s")
p1a <- ggplot() +
  geom_segment(data = seg, aes(x = x, y = y + .08, xend = sx, yend = sy - .12), color = MUT, linewidth = .7) +
  geom_label(data = sch, aes(x, y, label = s), fill = "grey94", label.size = 0, size = 4.2, fontface = "bold") +
  geom_point(data = stu, aes(x, y), size = 4, color = BLU) +
  annotate("text", x = 5, y = .62, label = "students (level 1)", color = BLU, size = 4.2) +
  annotate("text", x = 5, y = 2.42, label = "schools (level 2)", color = INK, size = 4.2) +
  scale_x_continuous(limits = c(.4, 9.6)) + scale_y_continuous(limits = c(.45, 2.6)) +
  labs(title = "Nested: students in schools", subtitle = "Each student belongs to one school") +
  theme_void(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 14, color = INK),
        plot.subtitle = element_text(size = 11.5, color = MUT))
grid <- expand_grid(stud = 1:4, itm = 1:5) |> mutate(hl = stud == 2 & itm == 3)
p1b <- ggplot(grid, aes(itm, stud)) +
  annotate("rect", xmin = .5, xmax = 5.5, ymin = 1.6, ymax = 2.4, fill = BLU, alpha = .12) +
  annotate("rect", xmin = 2.6, xmax = 3.4, ymin = .5, ymax = 4.5, fill = ORANGE, alpha = .18) +
  geom_point(aes(color = hl, size = hl)) +
  scale_color_manual(values = c(`FALSE` = MUT, `TRUE` = RED), guide = "none") +
  scale_size_manual(values = c(`FALSE` = 3, `TRUE` = 5), guide = "none") +
  scale_x_continuous(breaks = 1:5, labels = paste("item", 1:5), position = "top") +
  scale_y_reverse(breaks = 1:4, labels = paste("student", 1:4)) +
  annotate("text", x = 3, y = 4.85, label = "red: student 2's response to item 3", color = RED, size = 3.9) +
  labs(x = NULL, y = NULL, title = "Crossed: students and items",
       subtitle = "Each response belongs to one student and one item") +
  theme(panel.grid = element_blank(), axis.ticks = element_blank(),
        panel.border = element_blank(), axis.text = element_text(size = 11.5, color = INK)) +
  coord_cartesian(ylim = c(4.95, .5), clip = "off")
save(p1a + p1b + plot_layout(widths = c(1, 1)), "fig01_nested_crossed.png", 10.6, 4.0)

## fig02 --- the real responses: 30 students by 36 items ------------------------------------
set.seed(808)
complete_ids <- vocab |> count(id) |> filter(n == 36) |> pull(id)
ids30 <- sample(complete_ids, 30)
item_order <- vocab |> group_by(item) |> summarise(pc = mean(resp)) |> arrange(desc(pc)) |> pull(item)
g2 <- vocab |> filter(id %in% ids30) |>
  group_by(id) |> mutate(nc = sum(resp)) |> ungroup() |>
  mutate(item = factor(item, levels = item_order),
         idf = fct_reorder(factor(id), nc),
         r = factor(resp, levels = c(1, 0), labels = c("correct (1)", "incorrect (0)")))
cat("fig02: number correct of the 30 students range", range(g2$nc), "\n")
p2 <- ggplot(g2, aes(item, idf, fill = r)) +
  geom_tile(color = "white", linewidth = .4) +
  scale_fill_manual(values = c("correct (1)" = BLU, "incorrect (0)" = "grey86"), name = NULL) +
  labs(x = "item (easiest on the left, hardest on the right)",
       y = "student (most correct at the top)",
       title = "Thirty students' responses to the 36 items") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())
save(p2, "fig02_response_grid.png", 8.6, 5.4)

## fig03 --- one curve per item (the cross-classified logit model as item curves) -------------
pick <- ie |> filter(level %in% c("vocab_q28", "vocab_q20", "vocab_q34")) |>
  left_join(vocab |> group_by(item) |> summarise(pc = mean(resp)), by = c("level" = "item")) |>
  mutate(lab = case_when(level == "vocab_q28" ~ "easiest item",
                         level == "vocab_q20" ~ "an item near average easiness",
                         TRUE ~ "hardest item"),
         lab = paste0(lab, " (item ", parse_number(level), ")"))
curves <- pick |> select(lab, estimate) |> expand_grid(uj = seq(-3.6, 3.6, length.out = 300)) |>
  mutate(p = plogis(g0 + uj + estimate))
avg <- pick |> mutate(p = plogis(g0 + estimate))
print(avg |> select(level, estimate, pc, p))
cols3 <- setNames(c(GRN, BLU, RED), sort(unique(pick$lab))[c(2, 1, 3)])
cols3 <- c(setNames(GRN, pick$lab[pick$level == "vocab_q28"]), setNames(BLU, pick$lab[pick$level == "vocab_q20"]),
           setNames(RED, pick$lab[pick$level == "vocab_q34"]))
lev3 <- pick$lab[match(c("vocab_q28", "vocab_q20", "vocab_q34"), pick$level)]
curves$lab <- factor(curves$lab, levels = lev3); avg$lab <- factor(avg$lab, levels = lev3)
p3 <- ggplot(curves, aes(uj, p, color = lab)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = MUT) +
  geom_line(linewidth = 1.2) +
  geom_point(data = avg, aes(x = 0, y = p), size = 3) +
  geom_text(data = avg, aes(x = 0, y = p, label = pct(p)), hjust = -.25, vjust = 1.1, size = 4, show.legend = FALSE) +
  scale_color_manual(values = cols3, name = NULL) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  guides(color = guide_legend(nrow = 3)) +
  labs(x = expression("student's "*u[j]*" (log-odds above or below a student of average proficiency)"),
       y = "probability of a correct response",
       title = "One curve per item, all with the same shape",
       subtitle = expression("Dashed line: a student of average proficiency ("*u[j] == 0*")"))
save(p3, "fig03_item_curves.png", 7.6, 5.6)

## fig04 --- the two sets of random effects (lecture PDF 14-15) -----------------------------
cat("SD person EB", sd(pe$estimate), " SD item EB", sd(ie$estimate), "\n")
p4a <- ggplot(pe, aes(estimate)) + geom_histogram(bins = 30, fill = BLU, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = RED) +
  labs(x = expression("student's estimated "*u[j]*" (log-odds)"), y = "number of students",
       title = "1,914 students", subtitle = sprintf("SD of the estimates %.2f", sd(pe$estimate)))
p4b <- ggplot(ie, aes(estimate)) + geom_histogram(bins = 15, fill = ORANGE, color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed", color = RED) +
  labs(x = expression("item's estimated "*u[i]*" (easiness, log-odds)"), y = "number of items",
       title = "36 items", subtitle = sprintf("SD of the estimates %.2f", sd(ie$estimate)))
save(p4a + p4b, "fig04_random_effects.png", 10.6, 4.4)

## fig05 --- estimates against the raw proportion correct ------------------------------------
st <- vocab |> group_by(id) |> summarise(nc = sum(resp), n = n(), pc = mean(resp)) |>
  mutate(level = as.character(id)) |> left_join(pe |> select(level, estimate, std.error), by = "level")
st36 <- st |> filter(n == 36) |> group_by(nc) |> summarise(estimate = first(estimate), k = n(), pc = first(pc))
stfew <- st |> filter(n < 36)
cat("fig05: complete students", sum(st36$k), " fewer than 36:", nrow(stfew), "\n")
print(st36 |> filter(nc %in% c(0, 35)))
it <- vocab |> group_by(item) |> summarise(pc = mean(resp), n = n()) |> left_join(ie |> select(level, estimate, std.error), by = c("item" = "level"))
p5a <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = MUT) +
  geom_point(data = st36, aes(pc, estimate, size = k), color = BLU, alpha = .85) +
  geom_point(data = stfew, aes(pc, estimate), shape = 21, color = RED, fill = "white", size = 2.6, stroke = 1.1) +
  scale_size_area(max_size = 5.5, name = "students", breaks = c(10, 50, 100)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "percent of the student's responses correct", y = expression("student's estimated "*u[j]),
       title = "Students", subtitle = "Blue: 36 responses; red: fewer than 36")
p5b <- ggplot(it, aes(pc, estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = MUT) +
  geom_point(color = ORANGE, size = 2.8) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = "percent of the item's responses correct", y = expression("item's estimated "*u[i]),
       title = "Items", subtitle = "About 1,900 responses each")
save(p5a + p5b, "fig05_estimates_vs_percent.png", 10.6, 4.6)

## fig06 --- the latent response (Section 4.1) -------------------------------------------------
lat_items <- ie |> filter(level %in% c("vocab_q20", "vocab_q34")) |>
  mutate(center = g0 + estimate, p = plogis(center),
         lab = ifelse(level == "vocab_q20", "Item 20 (near average easiness)", "Item 34 (hardest)"))
print(lat_items |> select(level, center, p))
lat <- lat_items |> select(lab, center, p) |> expand_grid(x = seq(-9, 9, length.out = 600)) |>
  mutate(dens = dlogis(x, location = center, scale = 1), above = x > 0)
lat_lab <- lat_items |> mutate(txt = paste0(pct(p), " of the area\nis above 0"))
p6 <- ggplot(lat, aes(x, dens)) +
  geom_area(data = filter(lat, above), fill = BLU, alpha = .35) +
  geom_line(color = INK, linewidth = .9) +
  geom_vline(xintercept = 0, linetype = "dashed", color = RED) +
  geom_vline(data = lat_items, aes(xintercept = center), linetype = "dotted", color = MUT) +
  geom_text(data = lat_lab, aes(x = 4.6, y = .2, label = txt), size = 3.9, color = BLU, lineheight = .95) +
  facet_wrap(~ lab, ncol = 2) +
  scale_x_continuous(breaks = seq(-8, 8, 4)) +
  labs(x = "latent response (log-odds scale)", y = "how common each value is",
       title = "Latent response of a student of average proficiency to two items",
       subtitle = "Correct (1) if the latent response is above 0 (dashed red line)")
save(p6, "fig06_latent_response.png", 10.6, 4.2)

## fig07 --- from the latent response to the probability curve (adapted from RH&S Vol. 2, Fig. 10.3)
# Item of average easiness (u_i = 0). Each student's latent response is a logistic curve, drawn on its side,
# centered at g0 + u_j; the shaded part above the cutoff 0 is the probability plotted in the top panel.
xs <- -3:3; kk <- 3                                   # kk scales the sideways densities
dens <- map_dfr(xs, function(s0) { cc <- g0 + s0; yy <- seq(cc - 5.5, cc + 5.5, length.out = 241)
  xx <- s0 + kk * dlogis(yy, cc); tibble(x0 = s0, y = yy, x = xx) })
shade <- map_dfr(xs, function(s0) { cc <- g0 + s0; yy <- seq(0, cc + 5.5, length.out = 200)
  xx <- c(s0 + kk * dlogis(yy, cc), rep(s0, length(yy))); tibble(x0 = s0, y = c(yy, rev(yy)), x = xx) })
pts <- tibble(x0 = xs, p = plogis(g0 + xs))
print(pts)
p7b <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = RED) +
  geom_vline(xintercept = -g0, linetype = "dotted", color = MUT) +
  geom_polygon(data = shade, aes(x, y, group = x0), fill = BLU, alpha = .35) +
  geom_segment(data = tibble(x0 = xs), aes(x = x0, xend = x0, y = g0 + x0 - 5.5, yend = g0 + x0 + 5.5), color = MUT, linewidth = .3) +
  geom_path(data = dens, aes(x, y, group = x0), color = INK, linewidth = .6) +
  geom_abline(intercept = g0, slope = 1, color = INK, linewidth = .9) +
  annotate("text", x = 3.9, y = -.5, label = "cutoff 0", color = RED, size = 3.6, hjust = 1, vjust = 1) +
  scale_x_continuous(breaks = -3:3, limits = c(-3.6, 3.9)) +
  coord_cartesian(ylim = c(-6.5, 9.2)) +
  labs(x = expression("student's "*u[j]*" (above or below average proficiency, log-odds)"),
       y = "latent response (log-odds)", title = "Latent response",
       subtitle = "One curve, on its side, per student; shaded: above 0")
p7a <- ggplot() +
  geom_hline(yintercept = .5, linetype = "dotted", color = MUT) +
  geom_vline(xintercept = -g0, linetype = "dotted", color = MUT) +
  geom_line(data = tibble(x = seq(-3.6, 3.6, length.out = 300)) |> mutate(p = plogis(g0 + x)), aes(x, p), color = BLU, linewidth = 1.1) +
  geom_point(data = pts, aes(x0, p), color = BLU, size = 2.6) +
  geom_text(data = pts, aes(x0, p, label = pct(p)), vjust = -1, size = 3.6, color = INK) +
  scale_x_continuous(breaks = -3:3, limits = c(-3.6, 3.9)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1.15), breaks = c(0, .25, .5, .75, 1)) +
  labs(x = NULL, y = "probability correct", title = "Probability of a correct response",
       subtitle = "Each point: the shaded share of the curve below it")
save(p7a / p7b + plot_layout(heights = c(1, 1.35)), "fig07_latent_curve.png", 7.6, 5.6)

## fig08 --- reliability as agreement between two versions of the test (simulation) ---------
# Version 1 is the 36 items of the data; version 2 is 36 new items drawn from N(0, tau_i^2).
# Seed 2 is the first seed for which both simulated correlations equal the formula values
# (0.30 and 0.94) to two decimals.
set.seed(2)
J <- 1914; u_sim <- rnorm(J, 0, tj)
uiA <- ie$estimate[ie$level == "vocab_q20"]; uiB <- ie$estimate[ie$level == "vocab_q28"]
oneA <- g0 + u_sim + uiA + rlogis(J); oneB <- g0 + u_sim + uiB + rlogis(J)
newi <- rnorm(36, 0, ti)
avgA <- g0 + u_sim + mean(ie$estimate) + rowMeans(matrix(rlogis(J * 36), J))
avgB <- g0 + u_sim + mean(newi) + rowMeans(matrix(rlogis(J * 36), J))
r1 <- cor(oneA, oneB); r36 <- cor(avgA, avgB); cat("fig08: r one item", r1, " r 36 items", r36, "\n")
p7a <- ggplot(tibble(a = oneA, b = oneB), aes(a, b)) + geom_point(alpha = .18, size = .7, color = BLU) +
  labs(x = "latent response to item A", y = "latent response to item B",
       title = "One item", subtitle = sprintf("Correlation across students: %.2f", r1))
p7b <- ggplot(tibble(a = avgA, b = avgB), aes(a, b)) + geom_point(alpha = .18, size = .7, color = BLU) +
  labs(x = "average over the 36 items of version 1", y = "average over 36 other items",
       title = "Average of 36 items", subtitle = sprintf("Correlation across students: %.2f", r36))
save(p7a + p7b, "fig08_two_versions.png", 10.6, 4.8)

## fig09 --- why the item variance is not in the reliability formula -------------------------
comp_ids <- vocab |> count(id) |> filter(n == 36) |> pull(id)
pe36 <- pe |> filter(as.numeric(level) %in% comp_ids)
qs <- quantile(pe36$estimate, c(.1, .3, .5, .7, .9))
five <- tibble(q = c("10th", "30th", "50th", "70th", "90th"), uj = as.numeric(qs)) |>
  mutate(lab = paste0(q, " percentile"), lab = factor(lab, levels = rev(lab)))
items_ord <- ie |> arrange(desc(estimate)) |> mutate(rank = row_number())
rk <- five |> expand_grid(items_ord |> select(rank, estimate)) |> mutate(lo = g0 + uj + estimate)
cols5 <- setNames(c(RED, ORANGE, MUT, GRN, BLU), levels(five$lab))
p8a <- ggplot(rk, aes(rank, lo, color = lab)) + geom_line(linewidth = 1.1) + geom_point(size = 1.3) +
  scale_color_manual(values = cols5, name = "student at the") +
  labs(x = "item rank (1 = easiest, 36 = hardest)", y = "predicted log-odds of a correct response",
       title = "Every item, five students", subtitle = "The lines never cross")
halves <- five |> mutate(easy = sapply(uj, function(u) mean(plogis(g0 + u + items_ord$estimate[1:18]))),
                         hard = sapply(uj, function(u) mean(plogis(g0 + u + items_ord$estimate[19:36])))) |>
  pivot_longer(c(easy, hard), names_to = "half", values_to = "pc") |>
  mutate(half = factor(half, levels = c("easy", "hard"), labels = c("18 easiest items", "18 hardest items")))
print(halves |> select(lab, half, pc))
p8b <- ggplot(halves, aes(half, pc, group = lab, color = lab)) + geom_line(linewidth = 1.1) + geom_point(size = 2.6) +
  scale_color_manual(values = cols5, guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = NULL, y = "predicted percent correct", title = "Easy half and hard half",
       subtitle = "Everyone lower on the hard half; same order")
save(p8a + p8b + plot_layout(widths = c(1.35, 1), guides = "collect"), "fig09_ranking.png", 10.6, 4.8)

## fig10 --- reliability as the number of items / students changes (lecture PDF 18-19) -------
relI <- expand_grid(tau = c(tj, 0.6), I = 1:60) |>
  mutate(rel = tau^2 / (tau^2 + 3.29 / I),
         lab = factor(ifelse(tau > 1, "ours", "alike"), levels = c("ours", "alike")))
lab6 <- c(ours = expression(hat(tau)[j] == 1.19~"(our students)"), alike = expression(tau[j] == "0.60"~"(students more alike)"))
ptsI <- tibble(I = c(1, 10, 36)) |> mutate(rel = tj^2 / (tj^2 + 3.29 / I))
print(ptsI); cat("tau .6 at 36:", .36/(.36 + 3.29/36), " at 10:", .36/(.36+.329), "\n")
p6a <- ggplot(relI, aes(I, rel, color = lab, linetype = lab)) +
  geom_hline(yintercept = .8, linetype = "dotted", color = MUT) +
  geom_line(linewidth = 1.1) +
  geom_point(data = ptsI, aes(I, rel), inherit.aes = FALSE, color = BLU, size = 3) +
  geom_text(data = ptsI, aes(I + c(5, 4, 0), rel + c(-.15, -.08, -.08), label = sprintf("%.2f", rel)), inherit.aes = FALSE,
            color = BLU, size = 4) +
  scale_color_manual(values = c(ours = BLU, alike = ORANGE), labels = lab6, name = NULL) +
  scale_linetype_manual(values = c(ours = "solid", alike = "dashed"), labels = lab6, name = NULL) +
  scale_y_continuous(limits = c(0, 1)) +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  labs(x = "number of items, I", y = "reliability", title = "Students' scores",
       subtitle = "Dotted line: reliability 0.80")
relP <- tibble(P = exp(seq(log(2), log(3000), length.out = 300))) |> mutate(rel = ti^2 / (ti^2 + 3.29 / P))
ptsP <- tibble(P = c(10, 30, 1914)) |> mutate(rel = ti^2 / (ti^2 + 3.29 / P))
print(ptsP)
p6b <- ggplot(relP, aes(P, rel)) +
  geom_hline(yintercept = .8, linetype = "dotted", color = MUT) +
  geom_line(linewidth = 1.1, color = ORANGE) +
  geom_point(data = ptsP, color = ORANGE, size = 3) +
  geom_text(data = ptsP, aes(label = c("0.79", "0.92", "0.9986")), color = ORANGE, hjust = c(-.3, -.3, .5), vjust = c(1.2, 1.2, 1.9), size = 4) +
  scale_x_log10(breaks = c(2, 10, 30, 100, 300, 1914), labels = c("2", "10", "30", "100", "300", "1,914")) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "number of students, P (log scale)", y = "reliability", title = "Items' easiness",
       subtitle = expression(hat(tau)[i] == 1.11))
save(p6a + p6b, "fig10_reliability.png", 10.6, 4.9)

## fig08 --- one LEP difference in log-odds, different differences in probability -------------
b2 <- fixef(m2)
ie2 <- tidy(m2, effects = "ran_vals") |> filter(group == "item")
d7 <- ie2 |> transmute(item = level, ui = estimate) |>
  mutate(item = fct_reorder(item, ui),
         lo0 = b2[1] + ui, lo1 = b2[1] + b2[2] + ui,
         p0 = plogis(lo0), p1 = plogis(lo1), gap = p0 - p1)
cat("fig07: gap in points range", range(d7$gap), " largest for item with p0", d7$p0[which.max(d7$gap)], "\n")
long7 <- d7 |> pivot_longer(c(lo0, lo1, p0, p1), names_to = "k", values_to = "v") |>
  mutate(scale = ifelse(startsWith(k, "lo"), "log-odds", "probability"),
         who = ifelse(endsWith(k, "0"), "not an English learner", "English learner"))
mk7 <- function(sc, ttl, xl) {
  dd <- long7 |> filter(scale == sc)
  w <- dd |> select(item, who, v) |> pivot_wider(names_from = who, values_from = v)
  ggplot(dd, aes(v, item)) +
    geom_segment(data = w, aes(x = `English learner`, xend = `not an English learner`, y = item, yend = item),
                 color = "grey70", linewidth = .8, inherit.aes = FALSE) +
    geom_point(aes(color = who), size = 2.1) +
    scale_color_manual(values = c("not an English learner" = BLU, "English learner" = ORANGE), name = NULL) +
    labs(x = xl, y = "item (easiest at the top)", title = ttl) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank())
}
p7a <- mk7("log-odds", "Log-odds", "predicted log-odds of a correct response")
p7b <- mk7("probability", "Probability", "predicted probability of a correct response") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1))
save(p7a + p7b + plot_layout(guides = "collect"), "fig11_lep_items.png", 10.6, 5.6)

## fig09 --- RT and log(RT) (lecture PDF 27-28) ----------------------------------------------
n200 <- sum(vocab$rt > 200); cat("fig08: responses over 200 s:", n200, " median", median(vocab$rt), "\n")
p8a <- ggplot(filter(vocab, rt <= 200), aes(rt)) + geom_histogram(binwidth = 4, boundary = 0, fill = BLU, color = "white") +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "response time (seconds)", y = "number of responses", title = "RT",
       subtitle = paste0("Right-skewed; ", scales::comma(n200), " responses over 200 s not shown"))
p8b <- ggplot(vocab, aes(ln_rt)) + geom_histogram(bins = 40, fill = BLU, color = "white") +
  scale_x_continuous(sec.axis = sec_axis(~ exp(.), name = "seconds", breaks = c(1, 3, 10, 30, 100, 300, 1000), labels = c("1", "3", "10", "30", "100", "300", "1,000"))) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = "log(RT)", y = "number of responses", title = "log(RT)", subtitle = "Much closer to symmetric")
save(p8a + p8b, "fig12_rt_hist.png", 10.6, 4.5)

## fig10 --- log(RT) by student and by item (lecture PDF 29-30) ------------------------------
p9a <- vocab |> filter(id <= 12) |>
  ggplot(aes(factor(id), ln_rt)) + geom_boxplot(fill = "grey94", outlier.size = 1, color = BLU) +
  labs(x = "student (ID 1 to 12)", y = "log(RT)", title = "By student")
p9b <- vocab |> filter(item_num <= 11) |>
  ggplot(aes(factor(item_num), ln_rt)) + geom_boxplot(fill = "grey94", outlier.size = .6, color = ORANGE) +
  labs(x = "item (0 to 11)", y = "log(RT)", title = "By item")
save(p9a + p9b, "fig13_rt_by_student_item.png", 10.6, 4.4)

## fig11 --- RT and accuracy (lecture PDF 31) and the naive model (PDF 32-34) ------------------
set.seed(31)
p10a <- ggplot(vocab, aes(c_ln_rt, resp)) +
  geom_jitter(alpha = .03, height = .03, width = 0, size = .5, color = MUT) +
  geom_smooth(color = BLU, linewidth = 1.1, method = "gam", formula = y ~ s(x, bs = "cs")) +
  geom_smooth(method = glm, method.args = list(family = binomial), formula = y ~ x, color = RED, linewidth = 1.1) +
  scale_y_continuous(breaks = c(0, .25, .5, .75, 1)) +
  labs(x = "centered log(RT)", y = "1 = correct response", title = "Responses and accuracy",
       subtitle = "Blue: no assumed shape; red: one logistic curve")
b3 <- fixef(m3); V3 <- vcov(m3)
nd <- tibble(x = seq(-3, 4, length.out = 200)) |>
  mutate(lo = b3[1] + b3[2] * x,
         se = sqrt(V3[1, 1] + x^2 * V3[2, 2] + 2 * x * V3[1, 2]),
         p = plogis(lo), lw = plogis(lo - 1.96 * se), up = plogis(lo + 1.96 * se))
p10b <- ggplot(nd, aes(x, p)) +
  geom_ribbon(aes(ymin = lw, ymax = up), fill = BLU, alpha = .18) +
  geom_line(color = BLU, linewidth = 1.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(x = "centered log(RT)", y = "predicted probability correct", title = "m3: one slope for log(RT)",
       subtitle = "Student of average proficiency,\nitem of average easiness; 95% band")
save(p10a, "fig14_rt_accuracy.png", 6.2, 4.6)
save(p10b, "fig15_m3_prediction.png", 5.8, 4.4)
cat("glm slope:", coef(glm(resp ~ c_ln_rt, vocab, family = binomial)), "\n")

## fig12 --- the three pieces of log(RT) (lecture PDF 40-41) -----------------------------------
pp <- vocab |> distinct(id, p_ln_rt); ii <- vocab |> distinct(item, i_ln_rt)
cat("SDs: d", sd(vocab$d_rt), " person (students)", sd(pp$p_ln_rt), " item (items)", sd(ii$i_ln_rt), "\n")
cat("d outside [-3,3]:", sum(abs(vocab$d_rt) > 3), " person means > 3:", sum(pp$p_ln_rt > 3), "\n")
p11a <- ggplot(filter(vocab, abs(d_rt) <= 3), aes(d_rt)) + geom_histogram(bins = 40, fill = BLU, color = "white") +
  scale_y_continuous(labels = scales::comma) + coord_cartesian(xlim = c(-3, 3)) +
  labs(x = expression(drt[ij]), y = "number of responses", title = "Double-centered",
       subtitle = sprintf("68,538 responses; SD %.2f", sd(vocab$d_rt)))
p11b <- ggplot(filter(pp, p_ln_rt <= 3), aes(p_ln_rt)) + geom_histogram(bins = 40, fill = GRN, color = "white") +
  coord_cartesian(xlim = c(-3, 3)) +
  labs(x = expression(bar(RT)[j]), y = "number of students", title = "Student means",
       subtitle = sprintf("1,914 students; SD %.2f", sd(pp$p_ln_rt)))
p11c <- ggplot(ii, aes(i_ln_rt)) + geom_histogram(binwidth = .1, fill = ORANGE, color = "white") +
  annotate("text", x = 2.12, y = 2.2, label = "item 0", size = 4, color = INK) +
  coord_cartesian(xlim = c(-3, 3)) +
  labs(x = expression(bar(RT)[i]), y = "number of items", title = "Item means",
       subtitle = sprintf("36 items; SD %.2f", sd(ii$i_ln_rt)))
set.seed(41); smp <- slice_sample(vocab, n = 4000)
p11d <- ggplot(smp, aes(p_ln_rt, d_rt)) + geom_point(alpha = .15, size = .8, color = MUT) +
  geom_smooth(method = lm, formula = y ~ x, color = RED, se = FALSE) +
  coord_cartesian(xlim = c(-2.5, 3), ylim = c(-3.5, 3.5)) +
  labs(x = expression(bar(RT)[j]), y = expression(drt[ij]), title = "Against the student mean",
       subtitle = sprintf("4,000 responses; correlation %.3f", cor(vocab$d_rt, vocab$p_ln_rt)))
save((p11a + p11b) / (p11c + p11d), "fig17_three_pieces.png", 10.6, 8.0)

## fig13 --- double centering for one student (the first rows of the lecture's table, PDF 39) ---
s909 <- vocab |> filter(id == 909) |> mutate(expected = p_ln_rt + i_ln_rt)
cat("student 909: n", nrow(s909), " p", unique(s909$p_ln_rt), " correct", sum(s909$resp), "\n")
p12 <- ggplot(s909, aes(item_num)) +
  geom_hline(yintercept = unique(s909$p_ln_rt), linetype = "dashed", color = GRN) +
  geom_hline(yintercept = 0, color = MUT, linewidth = .4) +
  geom_segment(aes(xend = item_num, y = expected, yend = c_ln_rt, color = d_rt > 0), linewidth = .9) +
  geom_point(aes(y = expected), shape = 95, size = 6, color = INK) +
  geom_point(aes(y = c_ln_rt, fill = factor(resp)), shape = 21, size = 2.8, color = INK) +
  scale_color_manual(values = c(`TRUE` = RED, `FALSE` = BLU), labels = c(`TRUE` = "slower than expected", `FALSE` = "faster than expected"), name = NULL) +
  scale_fill_manual(values = c(`1` = INK, `0` = "white"), labels = c(`1` = "correct", `0` = "incorrect"), name = NULL) +
  scale_x_continuous(breaks = seq(0, 35, 5)) +
  guides(color = guide_legend(nrow = 2), fill = guide_legend(nrow = 2)) +
  labs(x = "item", y = "centered log(RT)",
       title = "Student 909: each response against its expected value",
       subtitle = "Black dash: student mean + item mean; green dashed: the student's mean (0.1955)")
save(p12, "fig16_one_student.png", 9.4, 5.6)

## fig14 --- opposite signs within and between students (lecture PDF 42-43) -------------------
b4 <- fixef(m4); sdp <- sd(pp$p_ln_rt)
cat("m4:", b4, " person SD", sdp, "\n")
studs <- tibble(p = c(-1, 0, 1) * sdp,
                lab = c("faster (1 SD below average)", "average", "slower (1 SD above average)"),
                lab2 = c("faster student", "average-speed student", "slower student"))
wl <- studs |> expand_grid(dd = seq(-1, 1, length.out = 20)) |>
  mutate(x = p + dd, y = b4[1] + b4[3] * p + b4[2] * dd)
bl <- tibble(x = seq(-1.6, 1.6, length.out = 20)) |> mutate(y = b4[1] + b4[3] * x)
p13a <- ggplot() +
  geom_line(data = bl, aes(x, y), color = RED, linewidth = 1.1) +
  geom_line(data = wl, aes(x, y, group = lab), color = GRN, linewidth = 1.2) +
  geom_point(data = studs, aes(p, b4[1] + b4[3] * p), color = RED, size = 3.6) +
  geom_text(data = studs, aes(p - 1.05, b4[1] + b4[3] * p + b4[2] * (-1), label = lab2), hjust = 1, size = 3.8, color = INK) +
  scale_x_continuous(limits = c(-2.75, 1.7)) +
  labs(x = "centered log(RT)", y = "predicted log-odds of a correct response",
       title = "Log-odds",
       subtitle = "Green: within (slope -0.17); red: between students (slope 0.16)")
pr <- studs |> expand_grid(d = seq(-2, 2, length.out = 100)) |>
  mutate(pc = plogis(b4[1] + b4[2] * d + b4[3] * p),
         lab = factor(lab, levels = studs$lab))
p13b <- ggplot(pr, aes(d, pc, color = lab)) + geom_line(linewidth = 1.2) +
  geom_text(data = filter(pr, d == max(d)), aes(label = lab2), hjust = -0.08, size = 3.8, show.legend = FALSE) +
  scale_color_manual(values = c(BLU, MUT, ORANGE) |> setNames(studs$lab), guide = "none") +
  scale_x_continuous(limits = c(-2, 5.2), breaks = -2:2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = expression(drt[ij]*" (double-centered log(RT))"), y = "predicted probability correct",
       title = "Probability", subtitle = expression("Item of average easiness and mean time, "*u[j] == 0))
save(p13a + p13b + plot_layout(widths = c(1.15, 1)), "fig18_simpson.png", 10.6, 4.4)
cat("done\n")
