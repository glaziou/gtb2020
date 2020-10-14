#' ---
#' title: TB prevalence and incidence in Mozambique
#' author: Philippe Glaziou
#' date: 07/1072020
#' output:
#'    html_document:
#'      mode: selfcontained
#'      toc: true
#'      toc_depth: 3
#'      toc_float: true
#'      number_sections: true
#'      theme: flatly
#'      highlight: zenburn
#'      df_print: paged
#'      code_folding: hide
#' ---

#' Last updated: `r Sys.Date()`
#'
#'
#' # Load libraries and data
#'
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(propagate))
suppressPackageStartupMessages(library(rjags))
suppressPackageStartupMessages(library(knitr))
library(here)

load(here('Rdata/est.Rdata'))    # WHO estimates
load(here('Rdata/cty.Rdata'))    # iso3 
load(here('Rdata/pop.Rdata'))    # pop estimates from UN Pop Division
load('../../../gtb2019/Rdata/incs.Rdata') # WHO incidence splits by age
load(here('Rdata/tb.Rdata'))     # case notifications
pshiv <-
  fread(here('input/psHIV.csv')) # HIV prevalence in notified vs prevalent
yr <- 2019                       # year of the prevalence survey


#' ## Prevalence of pulmonary TB in adults (b+ pulmonary), received from Irwin on 8/7/2020: 
#' 3.34 (2.52 - 4.16)/1000 - including post-stratification adjustment
#' 
pa <- c(334, (416 - 252) / 3.92) / 1e5

#' Received from Irwin on 6/7/2020:
#'
#' Number of survey cases = 89
#' Number of survey cases on current tx = 7
#' 
#' Number of survey cases tested for HIV = 43 (actual tested)
#' Number or survey cases that are HIV+ (tested) = 5 (out of 43 tested)
#' Number or survey cases that are HIV+ (tested positive or verbally acknowledged) = 5+15=20 (out of 64 documented)
#' 
#' Number of participants that tested for HIV = 27,052 eligible to test of whom 15,651 accepted to test
#' Number of participants that tested for HIV that are HIV+= 804
#' Number of participants that are HIV+ (tested positive or verbally acknowledged)=2966
#' 

btreated <- 7
cases <- 89

#' ## Prevalence of HIV among prevalent TB
#' 
phiv <- 20 / 64
phiv.sd <- binom.test(20, 89)$conf.int
phiv.sd <- (phiv.sd[2] - phiv.sd[1]) / 3.92



#' sensitivity of xpert pooled (BGD, KEN, PHL, VNM) 
#' 
#' (source: GDG meeting, Dec 2019, Geneva)
#'
#' Pooled sensitivity and standard deviation against culture, using one xpert test:
(se <- c(0.73, (0.82 - 0.62) / 3.92))

#' Assume 7% increase in overall sensitivity after repeating the test in a second sample
#' 
#' based on observed 7% increase in Kenya using Xpert MTB-Rif 
#' 
(se <- se * 1.07)



#' ## Functions

#' ### Ratio of two random variables using Taylor expansion
#'
divXY <- function(X, Y, varX, varY, covXY = 0) {
  eXY <- X / Y - covXY / Y ^ 2 + X * varY / Y ^ 3
  varXY <-
    varX / Y ^ 2 - 2 * X * covXY / Y ^ 3 + X ^ 2 * varY / Y ^ 4
  return(list("E(X/Y)" = eXY, "Var(X/Y)" = varXY))
}

#' ### Beta parameters using the method of moments
#'
get.beta <- function(ev, sd) {
  #' @param ev expected value
  #' @param sd standard deviation
  #' @export
  stopifnot(ev > 0 & ev < 1)
  stopifnot(sd > 0)
  
  S = (ev * (1 - ev) / sd ^ 2) - 1
  if (S < 0)
    stop('Not distributed Beta: sd^2 >= ev*(1-ev)')
  
  a = S * ev
  b = S * (1 - ev)
  return(c(a = a, b = b))
}

#' ### Generate low and high bounds assuming Beta distribution
#'
lohi <- function(ev, sd) {
  #' @param ev expected value
  #' @param sd standard deviation
  #' @export
  stopifnot(ev > 0 & ev < 1)
  stopifnot(sd > 0)
  
  par <- get.beta(ev, sd)
  lo <- qbeta(0.025, par[1], par[2])
  hi <- qbeta(0.975, par[1], par[2])
  return(c(lo = lo, hi = hi))
}

vlohi <- Vectorize(lohi, c('ev', 'sd'))


#' ### Generate posterior prevalence (accounting for sensitivity)
#'
aprev <- function(tpos, n, se, se.sd) {
  #' bayesian post-hoc prevalence accounting for sensitivity
  #'
  #' @param tpos number testing positive
  #' @param n number tested
  #' @param se sensitivity
  #' @param se.sd SD of se
  #' @export
  require(rjags)
  
  stopifnot(tpos >= 0 | n > 0)
  stopifnot(!is.na(se) | !is.na(se.sd))
  stopifnot(se > 0 & se < 1)
  stopifnot(se.sd > 0 & se.sd < 1)
  
  se.par <- get.beta(se, se.sd)
  
  # model specs
  m <- "model {
  tpos ~ dbin(theta, n)
  theta <- se*phi
  se ~ dbeta(se1, se2)
  phi ~ dbeta(1, 1)
}"
  
  out <- jags.model(
    textConnection(m),
    data = list(
      tpos = tpos,
      n = n,
      se1 = se.par[1],
      se2 = se.par[2]
    ),
    n.chains = 3
  )
  
  update(out, 10000)
  
  mcmc <-
    coda.samples(out, variable.names = c("phi"), n.iter = 20000)
  omcmc <- summary(mcmc)
  post <- omcmc$statistics[1]
  post.lo <- omcmc$quantiles[1]
  post.hi <- omcmc$quantiles[5]
  post.sd <- omcmc$statistics[2]
  
  return(list(
    mcmc = mcmc,
    res = c(
      post = post,
      post.lo = post.lo,
      post.hi = post.hi,
      post.sd = post.sd
    )
  ))
}


#' ### Statistical ensemble
#'
#' The rate $R$ obtained using method i is assumed
#' distributed Beta with shape and scale parameters $\alpha_i+1$ and
#' $\beta_i+1$, respectively, and determined using the method of moments:
#' $R_i \sim B(\alpha_i+1,\beta_i+1)$ so that
#'
#'
#' $\textrm{Prob}(x = \textrm{TB}) = \int_{0}^{1} x B(\alpha_i, \beta_i)\, dx = \frac{\alpha_i + 1}{\alpha_i + \beta_i + 2}$
#'
#' The combined probability is then expressed as
#'
#' $\textrm{Prob}(x = \textrm{TB}) = \frac{\sum{\alpha_i} + 1}{\sum{\alpha_i} + \sum{\beta_i} + 2}$
#'
#' $\textrm{Var} = \frac{(\sum{\alpha_i+1})(\sum{\beta_i+1})}{(\sum{\alpha}+\sum{\beta+2})^2(\sum{\alpha}+\sum{\beta+3})}$
#'
#'
ensbeta <- function(xi, xi.sd) {
  stopifnot(xi < 1 & xi.sd < 1)
  stopifnot(xi > 0 & xi.sd > 0)
  vget.beta <- Vectorize(get.beta, c('ev', 'sd'))
  
  w <- vget.beta(xi, xi.sd) - 1
  a <- sum(w[1, ])
  b <- sum(w[2, ])
  pw <- list(c = a + 1, d = b + 1)
  k <- pw$c + pw$d
  post.ev <- pw$c / k
  post.lo <- qbeta(0.025, pw$c, pw$d)
  post.hi <- qbeta(0.975, pw$c, pw$d)
  post.sd <- sqrt(pw$c * pw$d / (k ^ 2 * (k + 1)))
  return(
    list(
      post.param = c(shape = pw$c, scale = pw$d),
      post.ev = post.ev,
      post.lo = post.lo,
      post.hi = post.hi,
      post.sd = post.sd
    )
  )
}



#' ### Generate incidence from prevalence
#'
prev2inc <- function(prev,
                     prev.sd,
                     prevk = NA,
                     newinc,
                     tbhiv,
                     tbhiv.sd,
                     rtbhiv = NA,
                     rtbhiv.sd = NA)
{
  #' function to derive incidence from prevalence
  #'
  #' @param prev prevalence per capita
  #' @param prev.sd standard deviation of prevalence
  #' @param prevk prevalence on tx (known cases)
  #' @param newinc detection rate (new+relapse) per capita
  #' @param tbhiv proportion HIV+ among incident cases
  #' @param tbhiv.sd SD of proportion HIV+
  #' @param rtbhiv HIV+ rate ratio (prevalent/incident)
  #' @param rtbhiv.sd SD of HIV+ rate ratio
  #' @export
  require(propagate) # use second-order Taylor expansion about moments
  
  stopifnot(prev > 0 & prev.sd > 0)
  stopifnot((prevk > 0 &
               prevk <= prev) | is.na(prevk))
  stopifnot(newinc >= 0)
  stopifnot(tbhiv > 0 & tbhiv < 1)
  stopifnot(tbhiv.sd > 0 & tbhiv.sd < 1)
  stopifnot((rtbhiv > 0 &
               rtbhiv < 1 &
               rtbhiv.sd > 0 &
               rtbhiv.sd < 1) | is.na(rtbhiv)) & !is.na(rtbhiv.sd)
  
  # durations
  Du.nh <-
    c(2.5, sqrt(3 / 4))           # not detected, HIV-neg ~U(1,4)
  Dn.nh <-
    c(1.1, sqrt(1.8 ^ 2 / 12))    # detected, HIV-neg ~U(0.2,2)
  Du.h <-
    c(0.105, sqrt(0.19 ^ 2 / 12)) # not detected, HIV-pos ~U(0.01,0.2)
  Dn.h <-
    c(0.505, sqrt(0.99 ^ 2 / 12)) # detected, HIV-pos ~U(0.01,1)
  
  Pr <- c(prev, prev.sd)
  Prk <- c(prevk, prev.sd * prevk / prev)
  Ni <- c(newinc, 0)
  H <- c(tbhiv, tbhiv.sd)
  r <- c(rtbhiv, rtbhiv.sd)
  if (is.na(rtbhiv) | is.na(rtbhiv.sd))
    r <- c(1, 0)
  known <- !is.na(prevk)
  
  if (known) {
    # $I = Sum_{i=u,n} \frac{P_i}{D_i}$
    # P denotes prevalence; D duration; n detected; u undetected
    p <-
      c(0.025, sqrt(0.05 ^ 2 / 12)) # self cures or dies before tx
    DT <- cbind(Pr, Prk, r, H, p, Dn.nh, Dn.h)
    EXPR <-
      expression((Pr - Prk) / ((1 - p) * (r * H * Dn.h + (1 - r * H) * Dn.nh)))
  } else if (Pr[1] > Ni[1] * ((1 - H[1]) * Dn.nh[1] + H[1] * Dn.h[1])) {
    # I = P_u/D_u + newinc, P_u > 0
    DT <-
      cbind(Pr, r, H, Ni, Du.nh, Dn.nh, Du.h, Dn.h)
    EXPR <-
      expression((Pr - Ni * ((1 - H) * Dn.nh + H * Dn.h)) / (r * H * Du.h + (1 -
                                                                               r * H) * Du.nh) + Ni)
    
  } else
    stop('prev too low compared with newinc')
  
  out <-
    propagate(
      expr = EXPR,
      data = DT,
      type = 'stat',
      do.sim = F,
      second.order = T
    )
  return(out)
}


#' ### Plot incidence trends
#'
iplot <- function(iso,
                  dta = est,
                  hiv = FALSE,
                  ylog = TRUE) {
  #' @param iso iso3 country code
  #' @param dta dataset
  #' @param hiv hiv lines?
  #' @param ylog log scale on y-axis?
  #' @export
  if (missing(iso) || !is.character(iso) || nchar(iso) != 3)
    stop("Error: iso must be a valid ISO3 country code")
  if (is.na(match('cty', objects(
    envir = parent.frame(), all.names = TRUE
  ))))
    stop("Error: cty must be loaded")
  iso <- toupper(iso)
  isocty <- function(iso)
    cty[iso, country]
  
  p <- qplot(
    year,
    inc,
    data = subset(dta, iso3 == iso),
    geom = 'line',
    main = paste('TB incidence in', as.character(isocty(iso)))
  ) +
    geom_ribbon(
      aes(year, ymin = inc.lo, ymax = inc.hi),
      fill = I('blue'),
      alpha = I(0.3)
    ) +
    geom_line(aes(year, newinc)) +
    xlab('') + ylab('Rate per 100,000/year') +
    theme_bw(base_size = 16)
  
  q <- p + geom_line(aes(year, inc.h), colour = I('red')) +
    geom_ribbon(
      aes(year, ymin = inc.h.lo, ymax = inc.h.hi),
      fill = I('red'),
      alpha = I(0.3)
    )
  
  if (ylog)
    p <-
    p + coord_trans(y = "log10") + ylab('Rate per 100,000/year (log scale)')
  
  if (hiv)
    return(q)
  else
    return(p)
}



#'
#' <br>
#'


#' # Pulmonary TB prevalence in 15+ year old
#'

#' Effective sample size: assume beta distribution, using the method of moments;
#'
#' the first beta parameters is the number of cases,
#'
#' the sum of the two parameters is the effective sample size.
#'
#' both parameters are rounded to the nearest integer and inputed into the
#'
#' bayesian model of prevalence
#'
(size <- as.integer(get.beta(pa[1], pa[2])))


#' Draw posterior distribution of prevalence
#' adjusted for sensitivity of xpert
#'
#+ echo=FALSE, message=FALSE, results='hide'
out <- aprev(
  tpos = size[1],
  n = sum(size),
  se = se[1],
  se.sd = se[2]
)

#' The trace plot
#' indicates good mixing of the markov chains:
#'
#+ fig.width=8, fig.height=6
plot(out$mcmc)


#' Posterior prevalence per 100000 pop:
(res = round(out$res * 1e5))

#' Small loss in relative precision:
#' coeff of variation (CV) of xpert positive prevalence:
(pa[2] / pa[1])

#' CV of bact positive prevalence:
#'
(res[4] / res[1])

#' <br>
#'

#' # Prevalence all forms and all ages
#'
#' Mozambique estimates and notification data
#' 
(moz <- est['MOZ', .(
  iso3,
  year,
  pop,
  inc,
  inc.lo,
  inc.hi,
  inc.sd,
  tbhiv,
  tbhiv.lo,
  tbhiv.hi,
  tbhiv.sd,
  newinc,
  newrel.hivpos,
  newrel.hivtest
)])

#' UN pop estimates (UNPD)
#' 
(pm <- pop['MOZ'][year %in% 2000:2019])

mozs <- incs['MOZ']


#' Proportion of children:
#' 
(ch <-
  pm[year == yr, .(ch = (e.pop.m014 + e.pop.f014) / e.pop.num)]$ch)

ch <- c(ch, 0)


#' ## Estimated incidence rate ratio
#'
#' incidence rates by ageXsex in Global TB Report 2019 
#' 
#' (assume same rate ratio in 2019 as estimated for 2018)
#' 
(mozs)

out <-
  divXY(mozs$inc.014,
        mozs$inc.15plus,
        mozs$inc.014.sd ^ 2,
        mozs$inc.15plus.sd ^ 2)
(r <- c(out[[1]], sqrt(out[[2]])))



#' ## Prop extrapulm, based on notifications
#'
#' Case notification data (EP and children)
#' 
(nm <- tb['MOZ',.(iso3, year, newinc, ep, ch)])

#' use data post 2014
#' 
ep <-
  nm[year > 2014, .(ep * (1 - ch)), by = year]

#' mean EP excluding children:
#' 
(emean <- mean(ep$V1))

#' SD:
#' 
(esd <- sd(ep$V1))
e <- c(emean, esd)



#' ## Overall correction factor for prevalence adults to all ages
#'
#' $\text{cf} = \frac{1-c + c r}{1-e}$
#'
(cf <- propagate(
  expression((1 - ch + ch * r) / (1 - e)),
  data = cbind(r, e, ch),
  do.sim = F,
  second.order = T
)$prop)


#' ## Prevalence all forms and all ages
#'
pra <- res[c(1, 4)]
(out <-
    propagate(
      expression(pra * (1 - ch + ch * r) / (1 - e)),
      data = cbind(r, e, ch, pra),
      do.sim = F,
      second.order = T
    )$prop)
pr <- out[2]
pr.sd <- out[4]
pr.ci <- c(out[5], out[6])



#' # Incidence
#'
#'

#' ## HIV risk ratio (among prevalent / among notified)
#'
#' ratio and SD

out <-
  divXY(phiv, last(moz$tbhiv), phiv.sd ^ 2, last(moz$tbhiv.sd) ^ 2)
(rhiv <- c(out[[1]], sqrt(out[[2]])))


#' ## Model 1
#'
(
  inc18a <- prev2inc(
    prev = pr / 1e5,
    prev.sd = pr.sd / 1e5,
    newinc = last(moz$newinc) / 1e5,
    tbhiv = last(moz$tbhiv),
    tbhiv.sd = last(moz$tbhiv.sd),
    rtbhiv = rhiv[1],
    rtbhiv.sd = rhiv[2]
  )$prop
)

#' ## Model 2
#'
(
  inc18b <- prev2inc(
    prev = pr / 1e5,
    prev.sd = pr.sd / 1e5,
    prevk = btreated / cases * pr / 1e5,
    newinc = last(moz$newinc) / 1e5,
    tbhiv = last(moz$tbhiv),
    tbhiv.sd = last(moz$tbhiv.sd),
    rtbhiv = rhiv[1],
    rtbhiv.sd = rhiv[2]
  )$prop
)


#' ## Statistical ensemble of M1 and M2
#'
p1 <- get.beta(inc18a[2], inc18a[4])
p2 <- get.beta(inc18b[2], inc18b[4])
d1 <- rbeta(1e5, p1[1], p1[2]) * 1e5
d2 <- rbeta(1e5, p2[1], p2[2]) * 1e5

(inc18 <- ensbeta(c(inc18a[2], inc18b[2]), c(inc18a[4], inc18b[4])))
d3 <- rbeta(1e5, inc18$post.param[1], inc18$post.param[2]) * 1e5


#+ fig.width=8, fig.height=6
qplot(d1, geom = 'density', colour = I('red')) +
  geom_density(aes(d2), colour = I('blue')) +
  geom_density(aes(d3), colour = I('black')) +
  xlab('Incidence rate per 100 000 (2018)') + ylab('Density') +
  theme_bw(base_size = 18) +
  annotate(
    'text',
    x = 300,
    y = 0.0028,
    label = 'Model 1',
    colour = I('red')
  ) +
  annotate(
    'text',
    x = 1000,
    y = 0.001,
    label = 'Model 2',
    colour = I('blue')
  ) +
  annotate(
    'text',
    x = 750,
    y = 0.0031,
    label = 'Ensemble',
    colour = I('black')
  )



#' Same average between the statistical ensemble and model 1.
#' 
#' We keep model 1. Model 2 is consistent but based on small numbers, it inflates uncertainty in the 
#' statistical ensemble.
#' 


#' ## Rescaling
#'
rf <- inc18a[2] * 1e5 / last(moz$inc)
rf.sd <- inc18a[4] * 1e5 / last(moz$inc.sd)


#' ## Generate incidence time series
#'
#' Reuseing previously published trends.
#' 
moz2 <- data.table(
  iso3 = 'MOZ',
  year = 2000:2019,
  pop = moz$pop,
  inc = moz$inc * rf,
  inc.sd = moz$inc.sd * rf.sd,
  newinc = moz$newinc
)

#' (newinc = notification rate, new + relpase cases)

out <- vlohi(moz2$inc / 1e5, moz2$inc.sd / 1e5)
moz2$inc.lo <- out[1, ] * 1e5
moz2$inc.hi <- out[2, ] * 1e5

moz2[, inc.num := inc * pop / 1e5]
moz2[, inc.lo.num := inc.lo * pop / 1e5]
moz2[, inc.hi.num := inc.hi * pop / 1e5]

knitr::kable(moz2)

#+ fig.width=8, fig.height=6, fig.cap='Updated incidence (blue) with notification rates (solid line)'
iplot('moz', moz2, ylog = F) +
  ggtitle('Updated TB incidence estimates in Mozambique') +
  expand_limits(y = 0)


#+ fig.width=8, fig.height=6, fig.cap='Before (red) and updated (blue), the solid line represents the case notification rate'
qplot(year,
      inc,
      data = moz2,
      geom = 'line',
      colour = I('blue')) +
  geom_ribbon(aes(year, ymin = inc.lo, ymax = inc.hi),
              fill = I('blue'),
              alpha = I(.3)) +
  geom_line(
    aes(year, inc),
    data = moz,
    colour = I('red'),
    size = I(1),
    linetype = I(2)
  ) +
  geom_ribbon(
    aes(year, ymin = inc.lo, ymax = inc.hi),
    data = moz,
    fill = I('red'),
    alpha = I(.2)
  ) +
  geom_line(aes(year, newinc)) +
  scale_y_log10(breaks = c(100, 200, 300, 400, 500, 700, 1000, 1200, 1500, 2000, 2500)) +
  xlab('') + ylab('Rate per 100 000 / year') +
  theme_bw(base_size = 16)


# save
ggsave(file = 'moz_before_after.pdf',
       width = 8,
       height = 6)

fwrite(moz2, file = 'moz_updated.csv')
