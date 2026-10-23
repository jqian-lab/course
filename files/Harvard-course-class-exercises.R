#' ---
#' title: "Truncated time-to-event data: in-class R exercises"
#' subtitle: "Harvard Catalyst Biostatistics Short Course, October 23, 2026"
#' author: "Jing Qian (qian@umass.edu)"
#' output:
#'   html_document:
#'     toc: true
#'     toc_float: true
#' ---
#'
#' Short examples to run during the class, one block per part. Each block runs on its
#' own in under a minute. Run the code, compare the output with the slide, then do the
#' **Try it** step below the block. The full course code is in `Harvard-course-R-code.R`.

#' # Setup ----
#' Needs survival, boot and gss (Overview, Parts I, II and IV) and tranSurv 1.2.3 or
#' later (Part III). coxrt and seqTrun are optional extras.
#+ setup, message = FALSE
pkgs <- c("survival", "boot", "gss", "tranSurv")
new_pkgs <- setdiff(pkgs, rownames(installed.packages()))
if (length(new_pkgs) > 0) install.packages(new_pkgs)
library(survival)


#' # Overview: who is in the risk set? ----
#' Five people; A dies at time 15. E only enters the study at time 17, so a death of E
#' at time 15 could never have been seen.
#+ overview-toy
toy <- data.frame(id    = c("A", "B", "C", "D", "E"),
                  entry = c(0, 3, 6, 9, 17),
                  exit  = c(15, 20, 24, 30, 26),
                  event = c(1, 0, 1, 0, 1))
naive <- survfit(Surv(exit, event) ~ 1, data = toy)          # ignores entry
fixed <- survfit(Surv(entry, exit, event) ~ 1, data = toy)   # delayed entry
summary(naive, times = 15)$n.risk   # 5 at risk: hazard 1/5
summary(fixed, times = 15)$n.risk   # 4 at risk: hazard 1/4

#' **Try it.** Change E's entry time from 17 to 12. How many are at risk at time 15 in
#' each analysis now, and why?
# Your code here:


#' ## Overview: GENIE BPC colorectal cancer (needs internet) ----
#' Stage IV colorectal cancer; time zero is diagnosis, and patients enter the database
#' when their tumor is sequenced (Brown et al., JAMA Oncol 2022; public data).
#+ overview-genie
url <- paste0("https://raw.githubusercontent.com/slb2240/",
              "delayed_entry_clin_genom_studies/main/crc_stage_iv_os_dx_data.csv")
crc <- read.csv(url)
crc <- with(crc, data.frame(dx_to_seq   = tt_cpt_report_mos,   # entry (months)
                            dx_to_death = tt_os_dx_mos,        # death or last contact
                            died        = os_status_dx))
survfit(Surv(dx_to_death, died) ~ 1, data = crc)              # naive: median 38.1
survfit(Surv(dx_to_seq, dx_to_death, died) ~ 1, data = crc)   # delayed entry: median 25.1

#' **Try it.** Move time zero to the day of sequencing: analyze
#' `seq_to_death = dx_to_death - dx_to_seq`. Is there still delayed entry? What
#' question does this median answer?
# Your code here:



#' # Part I: Channing House (left truncation) ----
#' Ages in months (converted to years). Time zero is birth; residents are observable
#' once they moved in alive. Despite its name, `cens = 1` means **died**.
#+ part1-channing
library(boot)   # channing
ch <- transform(channing, entry = entry/12, exit = exit/12)
fw <- subset(ch, entry < exit & sex == "Female")
naive <- survfit(Surv(exit, cens) ~ 1, data = fw)
fixed <- survfit(Surv(entry, exit, cens) ~ 1, data = fw)
naive   # median age at death 88.0
fixed   # median age at death 84.8
plot(naive, conf.int = FALSE, lty = 2, xlim = c(60, 102),
     xlab = "Age (years)", ylab = "Survival")
lines(fixed, conf.int = FALSE, col = "red", lwd = 2)
legend("bottomleft", c("Naive", "Delayed entry"), lty = c(2, 1), col = c(1, 2), bty = "n")

#' **Try it.** (1) Repeat the delayed-entry analysis for **men**. Look at the first rows of
#' `summary()`: how many men are at risk at the first two deaths, and what happens to the
#' curve? (2) Refit with `start.time = 68`. What does the median now estimate?
# Your code here:


#' ## Part I: right truncation, transfusion AIDS ----
#' `incu`: infection to AIDS ($X^*$); `infe`: infection to the July 1986 cutoff ($T^*$).
#' A case is in the data only if $X^* \le T^*$. Reverse time, then use delayed entry.
#+ part1-aids
library(gss); data(aids)   # months
shift <- max(aids$infe)    # reversal constant c (any c >= 99)
fit <- survfit(Surv(shift - infe - 0.5, shift - incu, rep(1, 295)) ~ 1, data = aids)
## Back to the original scale: F(x | X* <= 99) = S_X((shift - x)-)
F_rt <- function(x) sapply(x, function(u) {
  k <- which(fit$time < shift - u)
  if (length(k) == 0) 1 else fit$surv[max(k)]
})
median(aids$incu)                         # naive median incubation: 29 months
(0:99)[min(which(F_rt(0:99) >= 0.5))]     # corrected: 68 months (given <= 99 months)

#' **Try it.** Use `shift <- 150` instead of `max(aids$infe)` and rerun the block. Does
#' the corrected median change? Why or why not?
# Your code here:



#' # Part II: Cox regression with delayed entry (flchain) ----
#' Residents aged 50+ enter at the age of their blood draw. Exposure: top decile of
#' free light chain (FLC). Time scale: age.
#+ part2-flchain
fl <- subset(flchain, futime > 0)             # drop entry >= exit (3 rows)
fl$age_exit <- fl$age + fl$futime / 365.25    # age at death or last contact
fl$high <- as.numeric(fl$flc.grp == 10)       # top FLC decile
m_naive <- coxph(Surv(age_exit, death) ~ high + sex, data = fl)        # HR 1.33
m_age   <- coxph(Surv(age, age_exit, death) ~ high + sex, data = fl)   # HR 2.27
m_study <- coxph(Surv(futime, death) ~ high + sex + age, data = fl)    # HR 2.13
round(rbind(naive = summary(m_naive)$conf.int["high", c(1, 3, 4)],
            delayed_entry = summary(m_age)$conf.int["high", c(1, 3, 4)],
            time_on_study = summary(m_study)$conf.int["high", c(1, 3, 4)]), 2)

#' **Try it.** (1) Why does ignoring entry age shrink the hazard ratio? Compare the
#' entry ages of the two FLC groups. (2) Using `m_age`, estimate the median age at death
#' for women **alive at 65**, by FLC group: `survfit(m_age, newdata = ..., start.time = 65)`.
# Your code here:


#' ## Part II: right truncation, IPW Cox (transfusion AIDS) ----
#' Do recipients aged $\le$ 59 develop AIDS faster? Keep the 116 people infected in 1983
#' or later. Each case is weighted by $1/\Pr(R > x)$, its inverse chance of being seen.
#+ part2-ipw
data(aids, package = "gss")
d <- data.frame(X = pmax(aids$incu, 0.5),    # infection -> AIDS (months)
                R = aids$infe,               # infection -> June 30, 1986
                young = as.numeric(aids$age <= 59))
d <- subset(d, 102 - R > 60)                 # infected in 1983 or later: n = 116
fR <- survfit(Surv(X, R, rep(1, nrow(d))) ~ 1, data = d, stype = 2)  # R left-truncated by X
d$w <- 1 / stepfun(fR$time, c(1, fR$surv))(d$X)                      # w = 1 / Pr(R > X)
coxph(Surv(X) ~ young, data = d)                    # naive: HR 1.83
coxph(Surv(X) ~ young + offset(log(w)), data = d)   # IPW: HR 3.24 (printed SE too small)

#' **Try it.** (1) Inspect the weights: what are the median and the largest weight, and
#' what share of the total weight do the top 5% of cases carry? (2) If `coxrt` is
#' installed, get the correct SE with `coxrt::coxph.RT(X ~ young, right = R, data = d)`.
# Your code here:



#' # Part III: is entry quasi-independent of death? ----
#' Channing House again, ages in **months**. tranSurv must be version 1.2.3 or later.
#+ part3-tranSurv
library(tranSurv)
packageVersion("tranSurv")
data(channing, package = "boot")
chan <- subset(channing, entry < exit)                     # n = 457, 175 deaths
with(chan, cKendall(entry, exit, cens))                    # tau 0.085, p = 0.074
with(chan, cKendall(entry, exit, cens, method = "IPW2"))   # tau 0.444, p < 0.001 (prints 0)
coxph(Surv(entry, exit, cens) ~ sex, data = chan)          # Cox QI: male 0.322
trReg(Surv(entry, exit, cens) ~ sex, data = chan)          # CK: a = -0.80, male 0.489

#' **Try it.** (1) Test quasi-independence separately for women and men. Do the pooled
#' and within-sex tests agree? (2) Fit the Cox2 version (about 15 seconds):
#' `trReg(..., method = "adjust", control = list(Q = 2))`.
# Your code here:



#' # Part IV: sequential truncation, the fix as weighting ----
#' Setup 2 of Betensky, Qian and Hou (2023): $L$ = onset age, $X = L + U$ = entry age,
#' $R$ = death age, observed only if $L < X < R$; $\Pr(X < R) = 0.15$. Weight each
#' person by $1/\Pr(R > x_j)$ (NPMLE) and compare with the two shortcuts.
#+ part4-weighting
sim_seq <- function(n, r.shape) {     # keep sampling until n people satisfy X < R
  out <- NULL
  while (NROW(out) < n) {
    l <- rgamma(n, 26, 0.5); x <- l + runif(n, 24, 34); r <- rgamma(n, r.shape, 0.5)
    out <- rbind(out, cbind(l, x, r)[x < r, , drop = FALSE])
  }
  out[1:n, ]
}
set.seed(2026)
n <- 200
dat <- sim_seq(n, r.shape = 32.4)     # Pr(X < R) = 0.15
t.a <- 60.43                          # largest L we expect to see; estimate F_L given L <= t.a
l <- dat[, "l"]; x <- dat[, "x"]; r <- dat[, "r"]
S_R <- function(fit) stepfun(fit$time, c(1, fit$surv))
w_x <- 1 / S_R(survfit(Surv(x, r, rep(1, n)) ~ 1))(x)   # NPMLE: 1 / Pr(R > x_j)
w_l <- 1 / S_R(survfit(Surv(l, r, rep(1, n)) ~ 1))(l)   # KM(L | L < R): 1 / Pr(R > l_j)
F_w <- function(u, w) {                                  # weighted CDF of L, given L <= t.a
  F <- sapply(c(u, t.a), function(v) sum(w * (l <= v)) / sum(w))
  F[seq_along(u)] / F[length(u) + 1]
}
u <- seq(30, t.a, by = 0.5)
truth <- pgamma(u, 26, 0.5) / pgamma(t.a, 26, 0.5)
plot(u, truth, type = "l", lwd = 3, col = "grey60", ylim = c(0, 1),
     xlab = "Age at onset, L", ylab = "CDF given L <= t.a")
lines(u, F_w(u, rep(1, n)), lty = 3, col = "blue")      # ECDF
lines(u, F_w(u, w_l), lty = 4, col = "purple")          # KM(L | L < R)
lines(u, F_w(u, w_x), lty = 2, lwd = 2, col = "red")    # NPMLE
legend("topleft", c("Truth", "ECDF", "KM(L | L < R)", "NPMLE"), bty = "n",
       lty = c(1, 3, 4, 2), lwd = c(3, 1, 1, 2), col = c("grey60", "blue", "purple", "red"))
## CDF at the true median age at onset (truth 0.50): ECDF 0.76, KM 0.76, NPMLE 0.43
round(c(ECDF = F_w(u, rep(1, n))[which.min(abs(truth - 0.5))],
        KM_LR = F_w(u, w_l)[which.min(abs(truth - 0.5))],
        NPMLE = F_w(u, w_x)[which.min(abs(truth - 0.5))]), 2)

#' **Try it.** Change to light truncation, `r.shape = 47.5` ($\Pr(X < R) = 0.80$) with
#' `t.a <- 69.78`, and rerun the block. How do the shortcuts compare with the NPMLE now?
#' (If any weight is `Inf`, the KM for $R$ hit 0; rerun with another seed.)
# Your code here:

