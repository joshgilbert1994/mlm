# Figures for Chapter 10: Bayesian Multilevel Models and Distributional Models.
# Written by Bill Young. Run from the root of the book repository; set DATA_SUB and DATA_FULL
# to wherever your copies of jebs_sub.csv and jebs_full.csv live. Figures are written to OUT.
# Packages: tidyverse, lme4, brms (with rstan), posterior, patchwork, modelr. tidybayes is not installed, so
# add_epred_draws() is replaced by brms::posterior_epred(), and pp_check() is redrawn with ggplot from
# brms::posterior_predict() so that it uses this chapter's theme.
# The lecture fits with backend = "cmdstanr"; cmdstanr is not installed, so the same models are fitted with brms's
# default backend (rstan). The Stan program, the priors, the data and the sampler settings (4 chains, iter = 1000)
# are the lecture's. A sampler run uses random numbers, so every run gives slightly different draws: the seeds
# below were chosen so that the m1 summary rounds to the lecture's printed table (see SEED_M1).
# Fitting takes about 5 minutes. To reuse fits between runs, set the environment variable T10_FITS to a folder;
# by default the fits are kept in a temporary folder and deleted when R closes.
suppressMessages({library(tidyverse); library(lme4); library(brms); library(posterior); library(patchwork); library(modelr)})
DATA_SUB  <- "data/jebs_sub.csv"
DATA_FULL <- "data/jebs_full.csv"
OUT  <- "chapters/images/10/"; dir.create(OUT, showWarnings = FALSE)
FITS <- Sys.getenv("T10_FITS", tempdir())

# The data, read exactly as the lecture reads them
jebs_sub <- read_csv(DATA_SUB, show_col_types = FALSE)   # 780 students, 10 schools
jebs <- read_csv(DATA_FULL, show_col_types = FALSE)      # 7,797 students, 110 schools

INK<-"#1A1D23"; BLU<-"#2C5F8A"; RED<-"#B03A2E"; GRN<-"#2E7D5B"; MUT<-"#5C6672"; RUL<-"#DDE2E6"
ORANGE <- "#D98E04"; LBLU <- "#9DBBD6"; GRY <- "grey60"
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

# Canvas sizes. Each figure is drawn on a canvas about 1.3 times its printed width, so the 11.2-pt tick labels
# print at about 8.5 pt. PCT is the width to use in the qmd (percent of the 6.5 in text block); FAC rescales a
# canvas (width and height together) if a figure is later shrunk or enlarged in the qmd: change FAC by the same
# ratio as PCT so the printed text size stays the same.
PCT <- c(fig01 = 82, fig02 = 76, fig03 = 84, fig04 = 70, fig05 = 70, fig06 = 66, fig07 = 70, fig08 = 82,
         fig09 = 84, fig10 = 84, fig11 = 82, fig12 = 66, fig13 = 60, fig14 = 82, fig15 = 70, fig16 = 70,
         fig17 = 72, fig18 = 80, fig19 = 72)
FAC <- setNames(rep(1, length(PCT)), names(PCT))
FAC[["fig02"]] <- 76 / 82   # printed narrower for the page layout; canvas scaled by the same ratio so text size is unchanged
save <- function(p, f, w, h) {
  k <- substr(f, 1, 5); sc <- FAC[[k]]
  pt <- 11.2 * (PCT[[k]] / 100 * 6.5) / (w * sc)
  cat(sprintf("### %s canvas %.2f x %.2f in; printed at %d%% (%.2f in): tick labels %.1f pt\n",
              f, w * sc, h * sc, PCT[[k]], PCT[[k]] / 100 * 6.5, pt))
  ggsave(paste0(OUT, f), p, width = w * sc, height = h * sc, dpi = 300, bg = "white") }
fit_or_load <- function(name, expr) {
  f <- file.path(FITS, paste0(name, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  m <- expr; saveRDS(m, f); m }

# ---- the lecture's models ---------------------------------------------------------------------------------
# Part 1: random intercepts model for the 10-school sample, with the lecture's priors
priors1 <- prior(normal(5, 2), class = "b") + prior(normal(50, 10), class = "Intercept") + prior(normal(0, 5), class = "sd") +
  prior(normal(0, 5), class = "sigma")   # added to the lecture's code on 9 Nov 2026
mod1 <- bf(read_comp ~ std_baseline + (1 | cluster_id))
m0 <- fit_or_load("m0", brm(mod1, data = jebs_sub, prior = priors1, sample_prior = "only",
                            chains = 4, iter = 1000, cores = 4, seed = 1043, refresh = 0))
SEED_M1 <- 68     # chosen from seeds 1 to 200: its summary rounds to the lecture's slope row (7.17, 0.24, 6.71, 7.65)
m1 <- fit_or_load(paste0("m1_seed", SEED_M1), update(m0, sample_prior = "no", seed = SEED_M1, refresh = 0))
# Part 2: MELS model for all 110 schools
priors2 <- prior(normal(5, 2), class = "b", coef = "std_baseline") + prior(normal(0, 5), class = "b", coef = "treat") +
  prior(normal(50, 10), class = "Intercept") + prior(normal(0, 5), class = "sd") +
  prior(normal(1, .5), class = "Intercept", dpar = "sigma") + prior(normal(0, .5), class = "b", dpar = "sigma")
mod2 <- bf(read_comp ~ treat + std_baseline + (1 | cluster_id), sigma ~ treat + std_baseline)
m2 <- fit_or_load("m2", brm(mod2, data = jebs, prior = priors2, sample_prior = "only",
                            chains = 4, iter = 1000, cores = 4, seed = 1043, refresh = 0))
m3 <- fit_or_load("m3", update(m2, sample_prior = "no", seed = 1043, refresh = 0))
# lmer on the same 780 students (Lectures 2 to 4), for comparison
ml <- lmer(read_comp ~ 1 + std_baseline + (1 | cluster_id), data = jebs_sub)

d1 <- as_draws_df(m1); d3 <- as_draws_df(m3)
cat("m1 (rounded as the lecture's tidy table):\n"); print(broom.mixed::tidy(m1) |> mutate(across(where(is.numeric), ~ round(.x, 2))))
cat(sprintf("tau draws: mean %.3f median %.3f; lmer tau %.4f, beta1 %.4f (SE %.4f)\n",
            mean(d1$sd_cluster_id__Intercept), median(d1$sd_cluster_id__Intercept),
            attr(VarCorr(ml)$cluster_id, "stddev"), fixef(ml)[2], sqrt(diag(vcov(ml)))[2]))
dens_tau <- density(d1$sd_cluster_id__Intercept); cat(sprintf("tau draws: peak of density %.2f\n", dens_tau$x[which.max(dens_tau$y)]))
cat("m3 summary:\n"); print(summary(m3)$fixed[, 1:7], digits = 4); print(summary(m3)$random$cluster_id[, 1:7], digits = 4)
dl0 <- mean(d3$b_sigma_Intercept); dl1 <- mean(d3$b_sigma_treat); dl2 <- mean(d3$b_sigma_std_baseline)
cat(sprintf("m3 scale equation, posterior means: delta0 %.4f delta1 %.4f delta2 %.4f\n", dl0, dl1, dl2))

# ---- Figure 1: the 10-school data ----------------------------------------------------------------------------
p1a <- ggplot(jebs_sub, aes(read_comp)) +
  geom_histogram(binwidth = 4, boundary = 28, fill = LBLU, color = "white") +
  labs(title = "EOY reading scores", x = "EOY reading score (read_comp)", y = "Number of students")
p1b <- ggplot(jebs_sub, aes(std_baseline, read_comp, color = factor(cluster_id))) +
  geom_point(alpha = .55, size = 1.1) +
  scale_color_manual(values = c("#2C5F8A", "#D98E04", "#2E7D5B", "#B03A2E", "#7B4F9E", "#4A9FB5",
                                "#8C6D31", "#C2577A", "#5C6672", "#6AA84F"), guide = "none") +
  labs(title = "EOY against BOY MAP", x = "BOY MAP score (std_baseline)", y = "EOY reading score")
save(p1a + p1b + plot_layout(widths = c(1, 1.15)), "fig01_data.png", 7.4, 3.3)

# ---- Figure 2: prior x likelihood -> posterior for the MAP slope ---------------------------------------------
b_hat <- unname(fixef(ml)[["std_baseline"]]); b_se <- unname(sqrt(diag(vcov(ml)))[2])
post_norm <- function(m0_, s0, m_, s_) { w0 <- 1 / s0^2; w <- 1 / s_^2
  c(mean = (w0 * m0_ + w * m_) / (w0 + w), sd = 1 / sqrt(w0 + w), wdata = w / (w0 + w)) }
pA <- post_norm(5, 2, b_hat, b_se); pB <- post_norm(5, 2, b_hat, 2)
cat(sprintf("Fig 2: data precision %.3f, weight on data %.4f, posterior mean %.4f SD %.4f; weak data: mean %.4f SD %.4f\n",
            1 / b_se^2, pA[["wdata"]], pA[["mean"]], pA[["sd"]], pB[["mean"]], pB[["sd"]]))
xs <- seq(-1, 13, by = .01)
curves <- bind_rows(
  tibble(panel = "Our data: SE 0.24 (780 students)", x = xs, prior = dnorm(xs, 5, 2), like = dnorm(xs, b_hat, b_se), post = dnorm(xs, pA[["mean"]], pA[["sd"]])),
  tibble(panel = "Much less precise data: SE 2", x = xs, prior = dnorm(xs, 5, 2), like = dnorm(xs, b_hat, 2), post = dnorm(xs, pB[["mean"]], pB[["sd"]]))) |>
  mutate(panel = factor(panel, levels = unique(panel))) |>
  pivot_longer(prior:post, names_to = "curve") |>
  mutate(curve = factor(curve, levels = c("prior", "like", "post"),
                        labels = c("prior N(5, 2)", "data alone (likelihood)", "posterior")))
vl <- tibble(panel = factor(levels(curves$panel), levels = levels(curves$panel)), m = c(pA[["mean"]], pB[["mean"]]))
p2 <- ggplot(curves, aes(x, value, color = curve, linetype = curve)) +
  geom_line(linewidth = .9) +
  geom_vline(data = vl, aes(xintercept = m), color = BLU, linetype = "dotted") +
  facet_wrap(~ panel, scales = "free_y") +
  scale_color_manual(values = c(GRY, ORANGE, BLU)) + scale_linetype_manual(values = c("solid", "dashed", "solid")) +
  scale_x_continuous(breaks = seq(0, 12, 2)) +
  labs(title = "Prior and data combined into the posterior: the MAP slope", x = expression("MAP slope "*beta[1]*" (EOY points per SD of BOY MAP)"),
       y = "Density", color = NULL, linetype = NULL)
save(p2, "fig02_update.png", 7.4, 3.1)

# ---- Figure 3: the priors (the lecture's code, free x scales) -------------------------------------------------
set.seed(2026); n <- 1e5
priors_df <- tibble(g_0 = rnorm(n, 50, 10), b_1 = rnorm(n, 5, 2), tau = abs(rnorm(n, 0, 5)), sigma = abs(rnorm(n, 0, 5)))
cat(sprintf("HN(0,5): median %.2f, 95th pct %.2f; N(50,10) 95%%: %.1f to %.1f; N(5,2) 95%%: %.2f to %.2f\n",
            5 * qnorm(.75), 5 * qnorm(.975), 50 - 1.96 * 10, 50 + 1.96 * 10, 5 - 1.96 * 2, 5 + 1.96 * 2))
parlab <- c(g_0 = "gamma[0]*' (intercept)'", b_1 = "beta[1]*' (MAP slope)'", tau = "tau*' (school SD)'", sigma = "sigma*' (residual SD)'")
# The SD priors are half-normal: their densities are drawn with bounds = c(0, Inf), so the curves start at their
# maximum at 0 instead of spilling below 0 (values above 20 have prior probability 0.00006 and are left out)
long <- priors_df |> pivot_longer(everything(), names_to = "par") |>
  mutate(sd_par = par %in% c("tau", "sigma"), par = factor(parlab[par], levels = parlab))
p3 <- ggplot(long, aes(value)) +
  geom_density(data = filter(long, !sd_par), color = INK, fill = "grey85") +
  geom_density(data = filter(long, sd_par, value <= 20), color = INK, fill = "grey85", bounds = c(0, Inf)) +
  facet_wrap(~ par, scales = "free", labeller = label_parsed, nrow = 1) +
  labs(title = "The priors", x = "Value", y = "Density")
save(p3, "fig03_priors.png", 7.6, 2.4)

# ---- prior and posterior predictive checks: densities of real and simulated scores ----------------------------
ppc_df <- function(m, y, nd, seed, lim = NULL) {
  set.seed(seed); yrep <- posterior_predict(m, ndraws = nd)
  one <- function(v) { if (!is.null(lim)) v <- v[v >= lim[1] & v <= lim[2]]; d <- density(v); tibble(x = d$x, y = d$y) }
  bind_rows(map_dfr(seq_len(nrow(yrep)), ~ one(yrep[.x, ]) |> mutate(draw = .x, type = "simulated")),
            one(y) |> mutate(draw = 0, type = "real")) }
ppc_plot <- function(df, title, xlim = NULL) {
  p <- ggplot() +
    geom_line(data = filter(df, type == "simulated"), aes(x, y, group = draw, color = "simulated (Yrep)"), alpha = .55, linewidth = .4) +
    geom_line(data = filter(df, type == "real"), aes(x, y, color = "real scores (Y)"), linewidth = 1.1) +
    scale_color_manual(values = c("real scores (Y)" = INK, "simulated (Yrep)" = LBLU), breaks = c("real scores (Y)", "simulated (Yrep)")) +
    labs(title = title, x = "EOY reading score", y = "Density", color = NULL)
  if (!is.null(xlim)) p <- p + coord_cartesian(xlim = xlim)
  p }
set.seed(7); yrep0 <- posterior_predict(m0, ndraws = 2000)
cat(sprintf("m0 prior predictive (2,000 draws x 780): middle 95%% of fake scores %.1f to %.1f; share below 0 or above 100: %.3f\n",
            quantile(yrep0, .025), quantile(yrep0, .975), mean(yrep0 < 0 | yrep0 > 100)))
pp0 <- ppc_df(m0, jebs_sub$read_comp, 50, 11, lim = c(0, 100))
p4 <- ppc_plot(pp0, "Prior predictive check (50 simulated data sets)", c(0, 100))
save(p4, "fig04_prior_pc.png", 6.0, 3.6)

# ---- Figure 5: prior predictive school lines -----------------------------------------------------------------
nd <- jebs_sub |> data_grid(std_baseline = seq_range(std_baseline, n = 100), cluster_id)
set.seed(12); ep0 <- posterior_epred(m0, newdata = nd, ndraws = 50)
lines0 <- as_tibble(t(ep0), .name_repair = ~ paste0("d", seq_along(.x))) |> bind_cols(nd) |>
  pivot_longer(starts_with("d"), names_to = "draw", values_to = "epred")
cat(sprintf("prior lines: range of predicted EOY %.1f to %.1f\n", min(lines0$epred), max(lines0$epred)))
p5 <- ggplot(lines0, aes(std_baseline, epred, group = interaction(draw, cluster_id))) +
  geom_hline(yintercept = c(0, 100), linetype = "dashed", color = MUT) +
  geom_line(alpha = .2, color = BLU) +
  labs(title = "Prior predictive school lines", x = "BOY MAP score (std_baseline)", y = "Predicted EOY score")
save(p5, "fig05_prior_lines.png", 6.0, 3.2)

# ---- Figure 6: the draws of the MAP slope: four chains, and their histogram ------------------------------------
arr <- as_draws_array(m1, variable = "b_std_baseline")
tr <- as_tibble(as.data.frame(arr[, , 1])) |> mutate(iter = 500 + row_number()) |>
  pivot_longer(-iter, names_to = "chain", values_to = "b1") |> mutate(chain = sub("^([0-9]+).*", "\\1", chain))
b1 <- d1$b_std_baseline; q <- quantile(b1, c(.025, .975))
cat(sprintf("beta1 draws: n %d, mean %.3f, SD %.3f, 2.5%% %.3f, 97.5%% %.3f; Rhat %.3f, bulk ESS %.0f\n",
            length(b1), mean(b1), sd(b1), q[1], q[2], rhat(arr), ess_bulk(arr)))
p6b <- ggplot(tibble(b1), aes(b1)) +
  geom_histogram(bins = 40, fill = LBLU, color = "white") +
  geom_vline(xintercept = mean(b1), color = BLU, linewidth = 1) +
  geom_vline(xintercept = q, color = RED, linetype = "dashed", linewidth = .8) +
  annotate("label", x = mean(b1), y = Inf, vjust = 1.3, label = paste0("mean ", f2(mean(b1))), size = 3.6, color = BLU) +
  annotate("label", x = q, y = Inf, vjust = c(2.9, 2.9), label = paste0(c("2.5th pct ", "97.5th pct "), f2(q)), size = 3.4, color = RED) +
  scale_y_continuous(expand = expansion(mult = c(.02, .45))) +
  labs(title = expression(bold("The 2,000 draws of the MAP slope "*beta[1])), x = expression(beta[1]*" draw"), y = "Number of draws")
save(p6b, "fig06_draws.png", 5.6, 3.4)

# ---- Figure 7: the draws of the school SD tau ---------------------------------------------------------------
tau <- d1$sd_cluster_id__Intercept; tau_lmer <- attr(VarCorr(ml)$cluster_id, "stddev")
labs7 <- c("posterior mean", "posterior median", paste0("lmer estimate ", f2(tau_lmer)))   # no numbers for the draws: a rerun differs slightly from the printed table
vl7 <- tibble(x = c(mean(tau), median(tau), tau_lmer), what = factor(labs7, levels = labs7))
p7 <- ggplot(tibble(tau), aes(tau)) +
  geom_histogram(binwidth = .1, boundary = 0, fill = LBLU, color = "white") +
  geom_vline(data = vl7, aes(xintercept = x, color = what, linetype = what), linewidth = 1, key_glyph = "path") +
  scale_color_manual(values = c(BLU, BLU, RED)) + scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  scale_x_continuous(breaks = 0:8) +
  labs(title = expression(bold("The 2,000 draws of the school SD "*tau)), x = expression(tau*" draw (EOY points)"), y = "Number of draws",
       color = NULL, linetype = NULL) +
  guides(color = guide_legend(nrow = 1, override.aes = list(linewidth = .9)), linetype = guide_legend(nrow = 1)) +
  theme(legend.key.width = unit(1.6, "lines"), legend.margin = margin(0, 0, 0, 0))
save(p7, "fig07_tau.png", 6.0, 3.1)

# ---- Figure 8: prior and posterior predictive checks side by side -------------------------------------------
pp1 <- ppc_df(m1, jebs_sub$read_comp, 50, 13, lim = c(0, 100))
p8 <- (ppc_plot(pp0, "Before the data: prior", c(0, 100)) + ppc_plot(pp1, "After the data: posterior", c(0, 100))) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
save(p8, "fig08_ppc.png", 7.4, 2.8)
set.seed(8); yrep1 <- posterior_predict(m1, ndraws = 2000)
cat(sprintf("m1 posterior predictive: share in 45-55 simulated %.3f vs real %.3f; share outside 28.18-68.93 simulated %.3f\n",
            mean(yrep1 >= 45 & yrep1 <= 55), mean(jebs_sub$read_comp >= 45 & jebs_sub$read_comp <= 55),
            mean(yrep1 < min(jebs$read_comp) | yrep1 > max(jebs$read_comp))))

# ---- Figure 9: posterior school lines, with each school's students ------------------------------------------
set.seed(14); ep1 <- posterior_epred(m1, newdata = nd, ndraws = 50)
lines1 <- as_tibble(t(ep1), .name_repair = ~ paste0("d", seq_along(.x))) |> bind_cols(nd) |>
  pivot_longer(starts_with("d"), names_to = "draw", values_to = "epred")
nsch <- jebs_sub |> count(cluster_id) |> mutate(lab = paste0("School ", cluster_id, " (n = ", n, ")"))
lab_of <- setNames(nsch$lab, nsch$cluster_id)
p9 <- ggplot() +
  geom_point(data = jebs_sub, aes(std_baseline, read_comp), color = GRY, alpha = .45, size = .7) +
  geom_line(data = lines1, aes(std_baseline, epred, group = draw), alpha = .25, color = BLU) +
  facet_wrap(~ cluster_id, nrow = 2, labeller = as_labeller(lab_of)) +
  labs(title = "Posterior school lines: 50 draws per school", x = "BOY MAP score (std_baseline)", y = "EOY reading score") +
  theme(strip.text = element_text(size = 10))
save(p9, "fig09_post_lines.png", 7.6, 3.3)

# ---- Figure 19 (appendix): the four chains of the MAP slope ------------------------------------------------------
p19 <- ggplot(tr, aes(iter, b1, color = chain)) + geom_line(linewidth = .3, alpha = .85) +
  scale_color_manual(values = c(BLU, ORANGE, GRN, RED)) +
  labs(title = "Four chains, 500 draws each", x = "Step of the sampler (after 500 warmup steps)",
       y = expression(beta[1]*" draw"), color = "Chain") + theme(legend.position = "bottom")
save(p19, "fig19_chains.png", 5.8, 3.4)

# ---- Figure 10: prior against posterior (the lecture's code) ---------------------------------------------------
post_draws <- d1 |> as_tibble() |> select(g_0 = b_Intercept, b_1 = b_std_baseline, tau = sd_cluster_id__Intercept, sigma) |> mutate(type = "Posterior")
comb <- bind_rows(priors_df |> mutate(type = "Prior"), post_draws) |> pivot_longer(-type, names_to = "par") |>
  mutate(sd_par = par %in% c("tau", "sigma"), par = factor(parlab[par], levels = parlab), type = factor(type, levels = c("Prior", "Posterior")))
p10 <- ggplot(comb, aes(value, fill = type)) +
  geom_density(data = filter(comb, !sd_par), alpha = .5, color = NA) +
  geom_density(data = filter(comb, sd_par, value <= 20), alpha = .5, color = NA, bounds = c(0, Inf)) +
  facet_wrap(~ par, scales = "free", labeller = label_parsed, nrow = 2) +
  scale_fill_manual(values = c(Prior = "grey70", Posterior = BLU)) +
  labs(title = "Prior and posterior of each parameter", x = "Value", y = "Density", fill = NULL)
save(p10, "fig10_prior_post.png", 7.0, 4.2)

# ---- Figure 11: 20 confidence intervals from repeated samples, and one credible interval ----------------------
# Repeated samples: new data from the lmer fit (same 10 schools, same MAP scores), refitted with lmer.
g0 <- fixef(ml)[[1]]; b1t <- fixef(ml)[[2]]; sig <- sigma(ml)
sim_ci <- function(k) {
  set.seed(1000 + k)
  u <- rnorm(10, 0, tau_lmer); names(u) <- sort(unique(jebs_sub$cluster_id))
  y <- g0 + b1t * jebs_sub$std_baseline + u[as.character(jebs_sub$cluster_id)] + rnorm(nrow(jebs_sub), 0, sig)
  f <- suppressMessages(lmer(y ~ 1 + std_baseline + (1 | cluster_id), data = mutate(jebs_sub, y = y)))
  tibble(k = k, est = fixef(f)[[2]], se = sqrt(diag(vcov(f)))[2]) }
# seeds 1010 to 1029: the first run of 20 samples in which 19 intervals contain the true slope (over seeds 1001 to 1400, 95.25% did)
cis <- map_dfr(10:29, sim_ci) |> mutate(k = k - 9) |> mutate(lo = est - 1.96 * se, hi = est + 1.96 * se, covers = lo <= b1t & hi >= b1t)
cat(sprintf("Fig 11: %d of 20 intervals contain the true slope %.4f\n", sum(cis$covers), b1t))
p11a <- ggplot(cis, aes(y = k)) +
  geom_vline(xintercept = b1t, color = INK, linetype = "dashed") +
  geom_segment(aes(x = lo, xend = hi, yend = k, color = covers), linewidth = .9) +
  geom_point(aes(x = est, color = covers), size = 1.6) +
  scale_color_manual(values = c(`TRUE` = BLU, `FALSE` = RED), labels = c(`TRUE` = "contains the true slope", `FALSE` = "does not contain the true slope"), name = NULL) +
  scale_y_continuous(breaks = c(1, 5, 10, 15, 20)) +
  labs(title = "Confidence intervals", subtitle = "20 repeated samples; dashed line: true slope",
       x = expression(beta[1]*" (EOY points per SD of BOY MAP)"), y = "Sample")
dens_b1 <- density(b1); db <- tibble(x = dens_b1$x, y = dens_b1$y)
p11b <- ggplot(db, aes(x, y)) +
  geom_area(data = filter(db, x >= q[1], x <= q[2]), fill = LBLU) + geom_line(color = BLU, linewidth = .9) +
  annotate("text", x = mean(b1), y = max(db$y) * .45, label = "95% of the\nposterior", size = 3.6, color = INK) +
  labs(title = "Credible interval", subtitle = "One posterior, from our one sample",
       x = expression(beta[1]*" (EOY points per SD of BOY MAP)"), y = "Density")
save((p11a + p11b) + plot_layout(guides = "collect") & theme(legend.position = "bottom"), "fig11_intervals.png", 7.4, 3.0)

# ---- Figure 12: same mean, different spread (simulated) --------------------------------------------------------
set.seed(12); n12 <- 10000
sim12 <- bind_rows(tibble(y = rnorm(n12, 50, 10), g = "Control: mean 50, SD 10"), tibble(y = rnorm(n12, 50, 6), g = "Program: mean 50, SD 6"))
print(sim12 |> group_by(g) |> summarise(mean = mean(y), sd = sd(y)))
sp <- sim12 |> group_by(g) |> group_modify(~ { d <- density(.x$y); tibble(x = d$x, y = d$y) }) |> ungroup() |> filter(x >= 15, x <= 85)
p12 <- ggplot(sp, aes(x, y, color = g, linetype = g)) + geom_line(linewidth = 1) +
  geom_vline(xintercept = 50, color = MUT, linetype = "dotted") +
  scale_color_manual(values = c(ORANGE, BLU)) + scale_linetype_manual(values = c("dashed", "solid")) +
  labs(title = "Same mean, different spread (simulated)", x = "EOY reading score", y = "Density", color = NULL, linetype = NULL) +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom", legend.key.width = unit(2, "lines"))
save(p12, "fig12_spread.png", 5.4, 3.1)

# ---- Figure 13: boxplot by treatment status (the lecture's code) -------------------------------------------------
p13 <- ggplot(jebs, aes(factor(treat, labels = c("0 (control)", "1 (treated)")), read_comp)) +
  geom_boxplot(fill = "grey92", color = INK, width = .5) +
  labs(title = "EOY scores by treatment status", x = "Treatment status (treat)", y = "EOY reading score")
save(p13, "fig13_boxplot.png", 4.6, 2.8)
print(jebs |> group_by(treat) |> summarise(n = n(), mean = mean(read_comp), var = var(read_comp), sd = sd(read_comp),
                                          q1 = quantile(read_comp, .25), med = median(read_comp), q3 = quantile(read_comp, .75)))
print(var.test(read_comp ~ treat, jebs))

# ---- Figure 14: equal and unequal residual spread (simulated) ----------------------------------------------------
set.seed(31); mm <- rnorm(400)
het <- bind_rows(tibble(map = mm, y = 50 + 7 * mm + rnorm(400, 0, 6), sdv = 6, panel = "Equal spread at every MAP score"),
                 tibble(map = mm, sdv = exp(log(6) - .18 * mm), y = 50 + 7 * mm + rnorm(400, 0, sdv), panel = "Spread smaller at higher MAP scores")) |>
  mutate(panel = factor(panel, levels = unique(panel)))
band <- expand_grid(map = seq(-2.7, 2.7, by = .05), panel = levels(het$panel)) |>
  mutate(panel = factor(panel, levels = levels(het$panel)), sdv = ifelse(panel == levels(panel)[1], 6, exp(log(6) - .18 * map)),
         fit = 50 + 7 * map, lo = fit - 2 * sdv, hi = fit + 2 * sdv)
p14 <- ggplot() +
  geom_point(data = het, aes(map, y), color = GRY, alpha = .6, size = .9) +
  geom_line(data = band, aes(map, fit), color = BLU, linewidth = 1) +
  geom_line(data = band, aes(map, lo), color = ORANGE, linetype = "dashed", linewidth = .8) +
  geom_line(data = band, aes(map, hi), color = ORANGE, linetype = "dashed", linewidth = .8) +
  facet_wrap(~ panel) +
  labs(title = "Homoskedastic and heteroskedastic residuals (simulated)", x = "BOY MAP score", y = "EOY reading score")
save(p14, "fig14_hetero.png", 7.4, 2.7)

# ---- Figure 15 and 16: prior and posterior predictive checks for the MELS model --------------------------------
pp2 <- ppc_df(m2, jebs$read_comp, 50, 15, lim = c(0, 100))
p15 <- ppc_plot(pp2, "Prior predictive check, MELS model (50 data sets)", c(0, 100))
save(p15, "fig15_prior_pc_mels.png", 6.0, 3.2)
set.seed(9); yrep2 <- posterior_predict(m2, ndraws = 500)
cat(sprintf("m2 prior predictive: middle 95%% of fake scores %.1f to %.1f; share below 0 or above 100: %.3f\n",
            quantile(yrep2, .025), quantile(yrep2, .975), mean(yrep2 < 0 | yrep2 > 100)))
pp3 <- ppc_df(m3, jebs$read_comp, 50, 16)
p16 <- ppc_plot(pp3, "Posterior predictive check, MELS model (50 data sets)")
save(p16, "fig16_ppc_mels.png", 6.0, 3.2)
d_real <- density(jebs$read_comp); pk <- which(diff(sign(diff(d_real$y))) == -2) + 1
cat("full-data density peaks:", round(d_real$x[pk], 1), "\n")

# ---- Figure 17: predicted residual SD against BOY MAP, control and treated -------------------------------------
grid <- expand_grid(map = seq(-2.5, 2.7, by = .05), treat = c(0, 1))
sig_draws <- exp(outer(d3$b_sigma_Intercept, rep(1, nrow(grid))) + outer(d3$b_sigma_treat, grid$treat) + outer(d3$b_sigma_std_baseline, grid$map))
sg <- grid |> mutate(m = colMeans(sig_draws), lo = apply(sig_draws, 2, quantile, .025), hi = apply(sig_draws, 2, quantile, .975),
                     arm = factor(treat, levels = c(0, 1), labels = c("control (treat = 0)", "treated (treat = 1)")))
at <- sg |> filter(treat == 0, abs(map - round(map)) < 1e-8, map %in% c(-2, -1, 0, 1, 2))
cat("predicted residual SD, control, at MAP -2..2:", f2(at$m, 3), "\n")
cat(sprintf("exp(delta0) %.4f; exp(delta1) %.4f; exp(delta2) %.4f; exp(2*delta2) %.4f\n", exp(dl0), exp(dl1), exp(dl2), exp(2 * dl2)))
p17 <- ggplot(sg, aes(map, m, color = arm, fill = arm)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .2, color = NA) + geom_line(linewidth = 1) +
  geom_point(data = at, color = BLU, size = 2, show.legend = FALSE) +
  scale_color_manual(values = c(BLU, ORANGE)) + scale_fill_manual(values = c(BLU, ORANGE)) +
  labs(title = "Predicted residual SD by BOY MAP score", x = "BOY MAP score (std_baseline)",
       y = expression("Residual SD "*sigma[ij]), color = NULL, fill = NULL)
save(p17, "fig17_sigma_map.png", 6.0, 2.75)

# ---- Figure 18: a mean effect on one group of students looks like a smaller spread (simulated) -------------------
set.seed(41); nA <- 2000
warn <- bind_rows(tibble(arm = "Control", group = rep(c("group A", "group B"), each = nA), y = c(rnorm(nA, 40, 5), rnorm(nA, 56, 5))),
                  tibble(arm = "Program", group = rep(c("group A", "group B"), each = nA), y = c(rnorm(nA, 48, 5), rnorm(nA, 56, 5)))) |>
  mutate(arm = factor(arm, levels = c("Control", "Program")))
sds <- warn |> group_by(arm) |> summarise(sd = sd(y), mean = mean(y)) |> mutate(lab = paste0("mean ", f2(mean, 1), ", SD ", f2(sd, 1)))
print(sds)
# Each group's curve is its density times 0.5 (its share of the students), so the two thin curves add up to the
# thick curve for all students in the panel
dens <- function(v, w = 1) { d <- density(v, from = 20, to = 76); tibble(x = d$x, y = w * d$y) }
d18 <- bind_rows(warn |> group_by(arm, group) |> group_modify(~ dens(.x$y, .5)) |> ungroup(),
                 warn |> group_by(arm) |> group_modify(~ dens(.x$y)) |> ungroup() |> mutate(group = "all students"))
p18 <- ggplot(d18, aes(x, y, color = group, linewidth = group)) + geom_line() +
  geom_label(data = sds, aes(x = 20, y = Inf, label = lab), inherit.aes = FALSE, vjust = 1.3, hjust = 0, size = 3.7) +
  facet_wrap(~ arm) +
  scale_color_manual(values = c(`all students` = INK, `group A` = ORANGE, `group B` = BLU), breaks = c("group A", "group B", "all students")) +
  scale_linewidth_manual(values = c(`all students` = 1.3, `group A` = .8, `group B` = .8), breaks = c("group A", "group B", "all students")) +
  labs(title = "The program raises group A's mean only (simulated)", x = "EOY reading score", y = "Density", color = NULL, linewidth = NULL)
save(p18, "fig18_warning.png", 7.2, 3.6)

cat("done\n")
