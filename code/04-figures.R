# Figures for Chapter 4: Random Slopes Models.
# Written by Bill Young. Run from the root of the book repository; set DATA
# to wherever your copy of jebs_full.csv lives. Figures are written to OUT.
suppressMessages({library(tidyverse); library(lme4); library(patchwork); library(MASS); select <- dplyr::select})
DATA <- "data/jebs_full.csv"
OUT  <- "chapters/images/04/"; dir.create(OUT, showWarnings = FALSE)
set.seed(2026)

jebs <- read_csv(DATA, show_col_types = FALSE) |>
  mutate(cluster_id = factor(cluster_id))

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

ri <- lmer(read_comp ~ std_baseline + (1|cluster_id), jebs)
rs <- lmer(read_comp ~ std_baseline + (std_baseline|cluster_id), jebs)
ols <- lm(read_comp ~ std_baseline*cluster_id, jebs)
v <- as.data.frame(VarCorr(rs))
t0<-v$sdcor[1]; t1<-v$sdcor[2]; rho<-v$sdcor[3]; cv<-v$vcov[3]
g0<-fixef(rs)[1]; g1<-fixef(rs)[2]
re <- ranef(rs)$cluster_id; names(re)<-c("u0","u1")
re$cid <- rownames(re)

jebs <- jebs |> mutate(ri_fit = fitted(ri), rs_fit = fitted(rs), ols_fit = fitted(ols))

## fig01 --- three models side by side -----------------------------------
sp <- function(yv, ttl, sub, showy) {
  g <- ggplot(jebs, aes(std_baseline, .data[[yv]], group = cluster_id)) +
    geom_line(alpha=.42, linewidth=.32, color=INK) +
    coord_cartesian(ylim=c(26,80)) +
    labs(x="BOY MAP score", y=if (showy) "predicted EOY reading" else NULL,
         title=ttl, subtitle=sub) +
    theme(plot.title=element_text(size=13), plot.subtitle=element_text(size=10.5))
  if (!showy) g <- g + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())
  g
}
cat("### fig01\n")
ggsave(paste0(OUT,"fig01_three_models.png"),
  (sp("ri_fit","Random intercepts","110 intercepts, one slope", TRUE) |
   sp("rs_fit","Random slopes","110 intercepts and 110 slopes", FALSE) |
   sp("ols_fit","Separate regressions","each school fit on its own", FALSE)),
  width=10.4, height=3.9, dpi=200, bg="white")

## fig02 --- anatomy of one school ---------------------------------------
r1 <- re |> filter(cid=="1")
xs <- seq(-2.2,2.0,.01)
i1 <- g0+r1$u0; s1 <- g1+r1$u1
an <- bind_rows(
  tibble(x=xs, y=g0+g1*xs, what="the average school"),
  tibble(x=xs, y=i1+s1*xs,  what="school 1"))
x0 <- -1.0; x1 <- 0.0                       # a one-unit run, left of the crossing
p2 <- ggplot(an, aes(x,y,color=what,linetype=what)) +
  geom_line(linewidth=1.3) +
  scale_color_manual(values=c("the average school"=MUT,"school 1"=BLU), name=NULL) +
  scale_linetype_manual(values=c("the average school"="22","school 1"="solid"), name=NULL) +
  # the intercept gap at BOY = 0
  annotate("segment", x=0, xend=0, y=g0, yend=i1, color=RED, linewidth=1.7,
           arrow=arrow(length=unit(.18,"cm"), ends="both", type="closed")) +
  annotate("point", x=c(0,0), y=c(g0,i1), color=c(MUT,BLU), size=2.8) +
  annotate("text", x=-0.15, y=(g0+i1)/2+3.2, hjust=1, color=RED, size=4.0,
           label="3.11 points apart\nat a BOY score of 0", lineheight=.95) +
  # one unit of run on each line, in separate windows so nothing overlaps
  annotate("segment", x=-2.0, xend=-1.0, y=g0+g1*-2.0, yend=g0+g1*-2.0, color=MUT, linewidth=.6, linetype="21") +
  annotate("segment", x=-1.0, xend=-1.0, y=g0+g1*-2.0, yend=g0+g1*-1.0, color=MUT, linewidth=1.5,
           arrow=arrow(length=unit(.15,"cm"), ends="both", type="closed")) +
  annotate("text", x=-0.88, y=g0+g1*-1.55, hjust=0, color=INK, size=3.8,
           label="one unit across,\n6.95 points up", lineheight=.95) +
  annotate("segment", x=0.6, xend=1.6, y=i1+s1*0.6, yend=i1+s1*0.6, color=BLU, linewidth=.6, linetype="21") +
  annotate("segment", x=1.6, xend=1.6, y=i1+s1*0.6, yend=i1+s1*1.6, color=BLU, linewidth=1.5,
           arrow=arrow(length=unit(.15,"cm"), ends="both", type="closed")) +
  annotate("text", x=1.35, y=i1+s1*0.6-2.9, hjust=.5, color=BLU, size=3.8,
           label="one unit across,\n4.91 points up", lineheight=.95) +
  coord_cartesian(xlim=c(-2.4,2.1), ylim=c(33,68)) +
  labs(x="BOY MAP score", y="predicted EOY score",
       title="What the two departures mean, for one school",
       subtitle="School 1 has a higher intercept than the average school and a flatter slope.")
cat("### fig02\n"); ggsave(paste0(OUT,"fig02_anatomy.png"), p2, width=8.6, height=5.4, dpi=200, bg="white")

## fig03 --- trellis, small to large --------------------------------------
pick <- c("2","96","50","26","30","52")
lab <- jebs |> count(cluster_id) |> filter(cluster_id %in% pick) |>
  mutate(lb = paste0("school ", cluster_id, "  (n = ", n, ")"))
tre <- jebs |> filter(cluster_id %in% pick) |> left_join(lab, by="cluster_id") |>
  mutate(lb = factor(lb, levels = lab$lb[order(lab$n)]))
p3 <- ggplot(tre, aes(std_baseline, read_comp)) +
  geom_point(alpha=.35, size=1.1, color=MUT) +
  geom_line(aes(y=ols_fit, color="separate regression", linetype="separate regression"), linewidth=.75) +
  geom_line(aes(y=ri_fit, color="random intercepts", linetype="random intercepts"), linewidth=.75) +
  geom_line(aes(y=rs_fit, color="random slopes", linetype="random slopes"), linewidth=.9) +
  scale_color_manual(values=c("separate regression"=INK,"random intercepts"=BLU,"random slopes"=RED), name=NULL) +
  scale_linetype_manual(values=c("separate regression"="42","random intercepts"="12","random slopes"="solid"), name=NULL) +
  facet_wrap(~lb, nrow=2) +
  labs(x="BOY MAP score", y="EOY score",
       title="The same three models, school by school",
       subtitle="Top row: small schools, where the separate regression is unstable and the random slopes line is\nshrunk almost to the average slope. Bottom row: large schools, where the separate regression and the\nrandom slopes line are nearly identical.")
cat("### fig03
"); ggsave(paste0(OUT,"fig03_trellis.png"), p3, width=10.5, height=6.4, dpi=200, bg="white")

## fig04 --- distribution of school slopes + PI ---------------------------
sl <- tibble(slope = g1 + re$u1)
p4 <- ggplot(sl, aes(slope)) +
  geom_histogram(bins=22, fill="grey80", color="white") +
  geom_vline(xintercept=g1, color=INK, linewidth=.9) +
  geom_vline(xintercept=c(g1-1.96*t1, g1+1.96*t1), color=RED, linetype="42", linewidth=.8) +
  annotate("text", x=g1+0.10, y=17.5, label=sprintf("average slope %.2f", g1), color=INK, size=4, hjust=0, vjust=0) +
  annotate("text", x=g1-1.96*t1, y=9, label="4.74", color=RED, size=4, hjust=1.15) +
  annotate("text", x=g1+1.96*t1, y=9, label="9.16", color=RED, size=4, hjust=-0.15) +
  coord_cartesian(xlim=c(1.7,10.2)) +
  labs(x="school's own slope on BOY MAP", y="number of schools",
       title="Every school gets its own slope",
       subtitle="Bars are the 110 school slopes the model estimates \u2014 already shrunk toward the average.\nDashed lines mark the prediction interval, which covers about 95% of the TRUE school slopes.")
cat("### fig04
"); ggsave(paste0(OUT,"fig04_slope_spread.png"), p4, width=8.2, height=4.8, dpi=200, bg="white")

## fig05 --- intercept-slope scatter --------------------------------------
p5 <- ggplot(re, aes(u0,u1)) +
  geom_hline(yintercept=0, color=MUT, linetype="dotted") +
  geom_vline(xintercept=0, color=MUT, linetype="dotted") +
  geom_point(alpha=.6, size=1.9, color=BLU) +
  geom_smooth(method="lm", se=FALSE, color=RED, linewidth=.9, formula=y~x) +
  labs(x=expression("intercept departure  "*hat(u)[0*j]), y=expression("slope departure  "*hat(u)[1*j]),
       title="Schools with higher intercepts tend to have flatter slopes",
       subtitle="Each dot is one school.")
cat("### fig05
"); ggsave(paste0(OUT,"fig05_rho_scatter.png"), p5, width=7.6, height=5.0, dpi=200, bg="white")

## fig06 --- what rho looks like ------------------------------------------
# tau0 = 3, tau1 = 1.4 puts the convergence point at about x = 1.9, inside the frame
Z <- MASS::mvrnorm(22, c(0,0), diag(2), empirical = TRUE)   # one standardized sample, reused
mk <- function(r, lab) {
  S <- matrix(c(9, r*3*1.4, r*3*1.4, 1.96), 2)
  u <- Z %*% chol(S)
  xs <- seq(-2.5,2.5,.05)
  map_dfr(1:22, ~tibble(g=.x, x=xs, y=(50+u[.x,1]) + (7+u[.x,2])*xs)) |> mutate(panel=lab)
}
L1 <- "negative correlation: the lines fan in"
L2 <- "no correlation: no pattern"
L3 <- "positive correlation: the lines fan out"
sim <- bind_rows(mk(-0.9,L1), mk(0,L2), mk(0.9,L3)) |> mutate(panel=factor(panel, levels=c(L1,L2,L3)))
p6 <- ggplot(sim, aes(x,y,group=g)) + geom_line(alpha=.5, linewidth=.42, color=INK) +
  geom_vline(xintercept=0, color=RED, linetype="42", linewidth=.55) +
  facet_wrap(~panel) + coord_cartesian(ylim=c(24,76)) +
  labs(x="predictor (the dashed line marks 0)", y="outcome",
       title="What the intercept-slope correlation looks like",
       subtitle="Simulated schools. The three panels have the same intercept variance and slope variance; only the correlation differs.")
cat("### fig06\n"); ggsave(paste0(OUT,"fig06_rho_meaning.png"), p6, width=10.4, height=3.95, dpi=200, bg="white")

## fig07 --- scale dependence ---------------------------------------------
xs7 <- seq(-2.2,2.2,.02)
lines7 <- map_dfr(1:110, ~tibble(g=.x, x=xs7, y=(g0+re$u0[.x]) + (g1+re$u1[.x])*xs7))
band <- tibble(x=xs7, mid=g0+g1*xs7, sd=sqrt(t0^2 + 2*cv*xs7 + t1^2*xs7^2)) |>
  mutate(lo=mid-sd, hi=mid+sd)
marks <- tibble(m=c(-1,0,1),
                lab=c("\u00b13.35", "\u00b12.80", "\u00b12.64"),
                rl =c("\u03c1 = -0.62", "\u03c1 = -0.34", "\u03c1 = +0.07"))
marks <- marks |> mutate(ytop = g0+g1*m + sqrt(t0^2+2*cv*m+t1^2*m^2))
p7 <- ggplot() +
  geom_line(data=lines7, aes(x,y,group=g), alpha=.16, linewidth=.28, color=INK) +
  geom_ribbon(data=band, aes(x=x, ymin=lo, ymax=hi), fill=BLU, alpha=.22) +
  geom_line(data=band, aes(x,mid), color=BLU, linewidth=1.1) +
  geom_segment(data=marks, aes(x=m, xend=m, y=g0+g1*m - sqrt(t0^2+2*cv*m+t1^2*m^2), yend=ytop),
               color=RED, linewidth=1.5, inherit.aes=FALSE) +
  geom_text(data=marks, aes(x=m, y=ytop+1.4, label=lab), color=RED, size=4.0, fontface="bold", inherit.aes=FALSE) +
  geom_text(data=marks, aes(x=m, y=ytop+3.0, label=rl),  color=RED, size=3.6, inherit.aes=FALSE) +
  coord_cartesian(ylim=c(31,69)) +
  labs(x="BOY MAP score", y="predicted EOY score",
       title="How far apart the schools are depends on where along the predictor it is measured",
       subtitle="Grey lines are the 110 schools. The band is one SD either side of the average line;
the red bars measure it at three BOY scores.")
cat("### fig07\n"); ggsave(paste0(OUT,"fig07_scale_dependence.png"), p7, width=9.4, height=5.5, dpi=200, bg="white")

## fig08 --- the variance parabola ----------------------------------------
xx <- seq(-2.2,2.2,.01)
par8 <- tibble(x=xx, sd=sqrt(t0^2 + 2*cv*xx + t1^2*xx^2))
pts <- tibble(x=c(-1,0,1), sd=sqrt(t0^2 + 2*cv*c(-1,0,1) + t1^2*c(-1,0,1)^2))
p8 <- ggplot(par8, aes(x,sd)) + geom_line(color=BLU, linewidth=1.2) +
  geom_point(data=pts, color=RED, size=2.8) +
  geom_text(data=pts, aes(label=sprintf("%.2f", sd)), color=RED, size=4,
            hjust=c(-0.35, 0.5, 0.5), vjust=c(0.35, -1.1, -1.1)) +
  geom_vline(xintercept=-cv/t1^2, color=GRN, linetype="42", linewidth=.7) +
  annotate("text", x=-cv/t1^2+.08, y=4.24, hjust=0, color=GRN, size=3.7,
           label=sprintf("smallest at %.2f", -cv/t1^2)) +
  coord_cartesian(ylim=c(2.5,4.45)) +
  labs(x="BOY MAP score", y="SD of the school departures at that score",
       title="How far apart the schools are, as a function of the predictor",
       subtitle="In a random intercepts model this would be a flat line. Once slopes vary it is a curve.")
cat("### fig08
"); ggsave(paste0(OUT,"fig08_variance_curve.png"), p8, width=8.2, height=4.6, dpi=200, bg="white")

## fig09 --- Mundlak vs REWB, as a forest plot ------------------------------
rim <- lmer(read_comp ~ std_baseline + sch_map + (1|cluster_id), jebs)
rir <- lmer(read_comp ~ demeaned_map + sch_map + (1|cluster_id), jebs)
rsm <- lmer(read_comp ~ std_baseline + sch_map + (std_baseline|cluster_id), jebs)
rsr <- lmer(read_comp ~ demeaned_map + sch_map + (demeaned_map|cluster_id), jebs)

# within slope = the coefficient on the student-level score in each model.
# between slope = the coefficient on sch_map under REWB, and the SUM of the two
# under Mundlak, whose SE needs the covariance between them.
est_ci <- function(m, which) {
  V <- as.matrix(vcov(m)); b <- fixef(m)
  if (which == "within") { e <- b[2]; se <- sqrt(V[2,2]) }
  else if (which == "between_rewb") { e <- b[3]; se <- sqrt(V[3,3]) }
  else { e <- b[2] + b[3]; se <- sqrt(V[2,2] + V[3,3] + 2*V[2,3]) }
  e <- unname(e); se <- unname(se)
  list(est = e, lo = e - 1.96*se, hi = e + 1.96*se)
}
cmp <- bind_rows(
  tibble(model="random intercepts", src="Mundlak", param="within slope",  !!!est_ci(rim,"within")),
  tibble(model="random intercepts", src="REWB",    param="within slope",  !!!est_ci(rir,"within")),
  tibble(model="random intercepts", src="Mundlak", param="between slope", !!!est_ci(rim,"between_mund")),
  tibble(model="random intercepts", src="REWB",    param="between slope", !!!est_ci(rir,"between_rewb")),
  tibble(model="random slopes",     src="Mundlak", param="within slope",  !!!est_ci(rsm,"within")),
  tibble(model="random slopes",     src="REWB",    param="within slope",  !!!est_ci(rsr,"within")),
  tibble(model="random slopes",     src="Mundlak", param="between slope", !!!est_ci(rsm,"between_mund")),
  tibble(model="random slopes",     src="REWB",    param="between slope", !!!est_ci(rsr,"between_rewb"))
) |>
  mutate(param = factor(param, levels = c("within slope","between slope")),
         src   = factor(src, levels = c("REWB","Mundlak")))

row9 <- function(which, ttl, sub) {
  ggplot(filter(cmp, param == which), aes(est, src, color = src)) +
    geom_linerange(aes(xmin = lo, xmax = hi), linewidth = .9) +
    geom_point(size = 2.9) +
    geom_text(aes(label = sprintf("%.3f", est)), vjust = -1.2, size = 3.9,
              color = INK, show.legend = FALSE) +
    facet_wrap(~ model, nrow = 1) +
    scale_color_manual(values = c("Mundlak" = BLU, "REWB" = RED), guide = "none") +
    scale_y_discrete(expand = expansion(add = .8)) +
    labs(x = NULL, y = NULL, title = ttl, subtitle = sub) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 11, color = MUT))
}
p9 <- row9("within slope", "Within slope",
           "Identical under random intercepts; still agreeing to 0.001 under random slopes.") /
      (row9("between slope", "Between slope",
            "Identical under random intercepts. Under random slopes they differ by 0.249 \u2014 small next to either interval.") +
         labs(x = "estimate, with 95% confidence interval")) +
  plot_annotation(
    title = "The two parameterizations agree, until the slope is allowed to vary",
    theme = theme(plot.title = element_text(face = "bold", size = 15, color = INK),
                  plot.background = element_rect(fill = "white", color = NA)))
cat("### fig09\n"); ggsave(paste0(OUT,"fig09_mundlak_rewb.png"), p9, width=9.4, height=5.4, dpi=200, bg="white")

## fig10 --- appendix: treated vs control spread --------------------------
rs_c <- lmer(read_comp ~ std_baseline + (std_baseline|cluster_id), filter(jebs, treat==0))
rs_t <- lmer(read_comp ~ std_baseline + (std_baseline|cluster_id), filter(jebs, treat==1))
sd_c <- as.data.frame(VarCorr(rs_c))$sdcor[1]; sd_t <- as.data.frame(VarCorr(rs_t))$sdcor[1]
hd <- tibble(grp=c("control schools","treated schools"),
             one=c(2.998368, 2.305683), two=c(sd_c, sd_t)) |>
  pivot_longer(c(one,two), names_to="how", values_to="sd") |>
  mutate(how = recode(how, one="from the single model with a random slope on treatment",
                           two="from fitting the two groups separately"))
p10 <- ggplot(hd, aes(grp, sd, fill=how)) +
  geom_col(position=position_dodge(.7), width=.6) +
  geom_text(aes(label=sprintf("%.2f", sd)), position=position_dodge(.7), vjust=-0.5, size=4, color=INK) +
  scale_fill_manual(values=c(BLU,MUT), name=NULL) +
  scale_y_continuous(expand=expansion(mult=c(0,.18))) +
  guides(fill=guide_legend(nrow=2)) +
  labs(x=NULL, y="SD of the school departures",
       title="Treated schools are more alike than control schools",
       subtitle="Two ways of estimating the same comparison. Both say the spread between schools\nis smaller where the program ran.")
cat("### fig10
"); ggsave(paste0(OUT,"fig10_treat_spread.png"), p10, width=8.0, height=5.2, dpi=200, bg="white")

## fig11 --- how far each slope moves, against school size ----------------
f <- jebs |> group_by(cluster_id) |> filter(n()>2, sd(std_baseline)>0) |>
  group_modify(~tibble(b=coef(lm(read_comp~std_baseline,.x))[2], nj=nrow(.x))) |> ungroup() |>
  mutate(cid=as.character(cluster_id)) |> left_join(re, by="cid") |> mutate(bm=g1+u1, move=bm-b)
p11 <- ggplot(f, aes(nj, move)) +
  geom_hline(yintercept=0, color=MUT, linetype="dotted") +
  geom_point(alpha=.65, size=2, color=BLU) +
  annotate("text", x=8, y=11.6, hjust=0, color=RED, size=3.7,
           label="school 2: three students, moved 11.6 points") +
  annotate("segment", x=7, xend=3.6, y=11.5, yend=11.6, color=RED, linewidth=.35) +
  labs(x="students in the school", y="model slope minus separate-regression slope",
       title="The model moves a slope further the fewer students the school has",
       subtitle="Small schools are shrunk heavily toward the average slope. Large schools hardly move at all.")
cat("### fig11
"); ggsave(paste0(OUT,"fig11_movement.png"), p11, width=8.2, height=4.8, dpi=200, bg="white")

cat("done\n")
