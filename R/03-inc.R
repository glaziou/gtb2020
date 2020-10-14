#' ---
#' title: TB incidence
#' author: Philippe Glaziou
#' date: 2020-06-30
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

#' # Preamble
#' (Last updated: `r Sys.Date()`)
#'
#' Prevalence of HIV among notified cases
#'
#' # Load libraries and data
#'
suppressMessages(library(data.table))
suppressMessages(library(imputeTS))
suppressMessages(library(propagate))
library(here)

rm(list = ls())

source(here('R/fun.R'))

load(here('Rdata/tb.Rdata'))
load(here('Rdata/cty.Rdata'))
load(here('Rdata/unaids.Rdata'))
load(here('Rdata/tbhiv.Rdata'))
load(here('Rdata/pop.Rdata'))
load(here('Rdata/grpmbr.Rdata'))
load(here('Rdata/old.Rdata'))
load(here('Rdata/ihme.Rdata'))

vlohi <- Vectorize(lohi, c('ev', 'sd'))
yr <- 2019
m <- 1e5

#' update inc rules
#'
rinc <- fread(here('input/inc/rules_inc.csv'))
setkey(rinc, iso3)


#' income groups
#'
income <- grpmbr[group.type == 'g_income']
setnames(income, 'group.name', 'g.income')


#' new incidence series
#'
(dim(tbhiv))
(dim(old))
est <-
  merge(
    tbhiv,
    old[, list(iso3, year, inc, inc.sd)],
    by = c('iso3', 'year'),
    all.x = TRUE,
    all.y = FALSE
  )
(dim(est))
est <-
  merge(
    est,
    pop[, .(iso3, year, pop = e.pop.num)],
    by = c('iso3', 'year'),
    all.x = TRUE,
    all.y = FALSE
  )
(dim(est))

# check missing values
sum(is.na(est$inc) &
      est$year < yr) == 0   # TRUE: only 2018 inc values are missing

(sum(is.na(est$newinc)))
(sum(is.na(est$newinc) & est$year == yr))




#' # Outlier values in notification tseries
#'
#' check for outliers in the notification series, in countries with pop>1e5,
#' defined as notifications in any one year > 3 times the average over 2000-2015
#'
# sel <- est$iso3=='ESP' & est$year==2019
# if(est$c.newinc[sel] < 1000){
#   est[sel, newinc := NA]
#   est[sel, c.newinc := NA]
#   }


est[, imp.newinc := newinc] # will hold imputed values

#' check for outliers
#'
(est[pop > 1e5, .(outlier = sum(imp.newinc > 3 * mean(imp.newinc, na.rm =
                                                                  T)) > 0), by = iso3][outlier == T])
(est['STP', .(iso3, year, imp.newinc)])
sel <- est$iso3 == "STP" & est$year == 2003
est$imp.newinc[sel] <- NA # reset outlier to missing
(est['MDA', .(iso3, year, imp.newinc)])
sel <- est$iso3 == "MDA" & est$year == 2003
est$imp.newinc[sel] <- NA # reset outlier to missing
sel <- est$iso3 == "KGZ" & est$year == 2003
est$imp.newinc[sel] <- NA # reset outlier to missing

#' list outliers in the notification series in countries with pop>1e5,
#'
(est[pop > 1e5, .(outlier = sum(imp.newinc == 0)), by = iso3][outlier ==
                                                                          T])

#' KHM (Aug - NTP mentioned that 2014 peak due to 5000 false pos in children)
#'
(est['KHM', .(iso3, year, imp.newinc)])
sel <- est$iso3 == 'KHM' & est$year == 2014
est$imp.newinc[sel] <- NA # reset outlier to missing


#' check for effect of the global case chasing strategy
#' rapid increases (change > 10%)
#'
est[, chg.newinc := c(NA, diff(newinc)) / newinc, by = iso3]
sel1 <-
  est$year >= yr - 1 &
  est$newinc > 10 & est$chg.newinc > 0.1 & est$pop > 1e6
(est[sel1, .(iso3, year, newinc, imp.newinc, chg.newinc)])

# est[iso3 == 'MHL' &
#       year == 2018, imp.newinc := NA] # account for mass screening



#' # Interpolation of missing notifications
#'
#' using Kalman smoothing on structural TS, where possible
#'
tmp <- copy(est)
interp <- c('SMR', 'MSR', 'VGB')

est[iso3 %in% interp, imp.newinc := na_interpolation(imp.newinc), by = iso3]
est[iso3 %ni% interp, imp.newinc := na_kalman(imp.newinc, type='trend'), by = iso3]

est[, test.ispos(imp.newinc)]



#' visually check imputations
#'
wr <- c('AMR', 'AFR', 'EMR', 'EUR', 'SEA', 'WPR')

#+ fig.width=16, fig.height=12
for (i in wr) {
  p <-
    qplot(year, newinc, data = est[g.whoregion == i], geom = 'point') +
    geom_line(aes(year, imp.newinc), colour = I('red')) +
    facet_wrap( ~ iso3, scales = 'free_y')
  suppressWarnings(print(p))
  suppressWarnings(ggsave(here(
    paste('output/checks/imputations', i, '_newinc.pdf', sep = '')
  ),
  width = 14,
  height = 8))
  
}




#' add GBD regions
#'
gbd <- grpmbr[group.type == 'g_gbd']
setnames(gbd, 'group.name', 'g.gbd')
est <- merge(est, gbd[, .(iso3, g.gbd)], by = 'iso3', all.x = TRUE)

#' mdr regions
#'
gest <- grpmbr[group.type == 'g_est']
eeur <- gest[group.name == 'EEUR', iso3]
est[, g.mdr := g.whoregion]
est[iso3 %in% eeur, g.mdr := 'EEU']

#' add income groups
#'
est <-
  merge(est, income[, .(iso3, g.income)], by = 'iso3', all.x = TRUE)
est$g.income <- as.character(est$g.income)

#' add HBCs
#'
hbc <- as.character(grpmbr[group.type == 'g_hb_tb']$iso3)
est$g.hbc <- est$iso3 %in% hbc
hbmdr <- as.character(grpmbr[group.type == 'g_hb_mdr']$iso3)
est$g.hbmdr <- est$iso3 %in% hbmdr
hbtbhiv <- as.character(grpmbr[group.type == 'g_hb_tbhiv']$iso3)
est$g.hbtbhiv <- est$iso3 %in% hbtbhiv





#' # Standard adjustment
#'
std <-
  function(x,
           ISO3,
           f.lo = 1,
           f.hi,
           ystart = 2000,
           h = 1,
           source = 'Standard adjustment',
           smooth = FALSE) {
    #' $I = f N$
    #'
    #' @param x datatable
    #' @param country iso3
    #' @param f low bound
    #' @param f.hi high bound
    #' @param h 1 - over diagnosis
    #' @param smooth MA with exponential weighting
    #' @export
    f <- mean(c(f.hi, f.lo))
    sel <- x$iso3 == ISO3 & x$year >= ystart
    
    if (smooth == FALSE) {
      x[sel, inc := imp.newinc * f * h]
    } else {
      x[sel, inc := na_ma(imp.newinc, k = 4, weighting = 'exponential') * f * h]
    }
    x[sel, inc.sd := inc * (f.hi - 1) / 3.92]
    x[sel, source.inc := source]
  }




lst <- rinc$iso3[!is.na(rinc$hi)]
tmp <- copy(est) # backup
for (i in lst) {
  est <-
    std(
      est,
      ISO3 = i,
      f.lo = rinc[i, lo],
      f.hi = rinc[i, hi],
      ystart = rinc[i, ystart],
      h = 1,
      source = 'Standard adjustment',
      smooth = ifelse(rinc[i, hi] >= 1.5 &
                        !is.na(rinc[i, hi]), TRUE, FALSE)
    )
}

est[!is.na(inc), test.ispos(inc)]




#' over reporting in reviewed countries of the former USSR
#'
sel <- est$iso3 == 'KAZ'
est$inc[sel] <-
  predict(loess(imp.newinc ~ year, data = est[sel], span = 0.6))
est$inc.sd[sel] <- est$inc[sel] * 0.2
est$inc.lo[sel] <- est$inc[sel] - 1.96 * est$inc.sd[sel]
est$inc.hi[sel] <- est$inc[sel] + 1.96 * est$inc.sd[sel]


sel <- est$iso3 == 'RUS'
est$inc[sel] <-
  predict(loess(imp.newinc ~ year, data = est[sel], span = 0.6))
est$inc.sd[sel] <- est$inc[sel] * 0.2
est$inc.lo[sel] <- est$inc[sel] - 1.96 * est$inc.sd[sel]
est$inc.hi[sel] <- est$inc[sel] + 1.96 * est$inc.sd[sel]


sel <- est$iso3 %in% c('RUS', 'KAZ')
est$source.inc[sel] <- 'Standard adjustment'





#' # Inventory studies and capture-recapture
#'

#' GBR underreporting 5-17% (2010 report) ave = 11%
#' 1/(1-0.11)=1.12; 1/(1-0.05)=1.05: 1/(1-0.17)=1.2
#' use non-smoothed series with imputed 2013
#' to capture the increase in incidence in 2011
#'
est$source.inc[est$iso3 %in% c('GBR', 'NLD')] <- 'Inventory'


#' past CR, add 2018 based on current trends
#'
lst <- c('IRQ', 'YEM', 'EGY', 'IDN')

for (i in lst) {
  sel <- est$iso3 == i
  est[sel, inc := na_kalman(inc, type='trend')]
  est[sel, inc.sd := na_kalman(inc.sd, type='trend')]
  est[sel, source.inc := 'Inventory']
}

#' CHN
#' 
#' correspondance received on 6/8/2018
#' 
#' underreporting 8.2% in 2015, changes in care delivery system
#' and reporting since then. NTP states that 2018 value <= 8.2%
#' 
#' inc_2018 = newinc /(1 - 0.0082)
#' 
sel <- est$iso3 == 'CHN' & est$year == 2019
est[sel, inc := newinc * 1.144]
est[sel, inc.sd := inc * .288/3.92]

# est[sel, inc := na_kalman(inc)]
# est[sel, inc.sd := na_kalman(inc.sd)]

sel2 <- est$iso3 == 'CHN' & est$year == 2018
f <- est[sel2, newinc / (1 - 0.082) / inc]
est[sel, inc := inc * f]
est[sel, inc.sd := inc.sd * f]
est[sel, source.inc := 'Inventory']


#' DEU, inventory study (report shared in Aug 2020)
#' 
sel <- est$iso3 == 'DEU'
est[sel, inc := newinc / 0.95]
est[sel, inc.sd := (newinc / 0.9 - newinc) / 3.92]
est[sel, source.inc := 'Inventory']





#' # Prevalence surveys
pr <- unique(old$iso3[old$source.inc == 'Prevalence'])
est[iso3 %in% pr & iso3 %ni% c('CHN'), source.inc := 'Prevalence']

#' new survey results
#' ZAF, NPL, SWZ, LSO, MOZ
#' (LSO and MOZ include 2019 data points)
#' 
lst <- c('ZAF', 'NPL', 'SWZ', 'LSO', 'MOZ')

zaf <- fread(here('input/inc/zaf_updated.csv'))
npl <- fread(here('input/inc/npl_updated.csv'))
swz <- fread(here('input/inc/swz_updated.csv'))
lso <- fread(here('input/inc/lso_updated.csv'))
moz <- fread(here('surveys/MOZ/moz_updated.csv'))

est['ZAF', inc := na_kalman(c(zaf$inc, NA))]
est['ZAF', inc.sd := na_kalman(c(zaf$inc.sd, NA))]
est['NPL', inc := na_kalman(c(npl$inc, NA))]
est['NPL', inc.sd := na_kalman(c(npl$inc.sd, NA))]
est['SWZ', inc := na_kalman(c(swz$inc, NA))]
est['SWZ', inc.sd := na_kalman(c(swz$inc.sd, NA))]

est['LSO', inc := lso$inc]
est['LSO', inc.sd := lso$inc.sd]
est['MOZ', inc := moz$inc]
est['MOZ', inc.sd := moz$inc.sd]

est[iso3 %in% lst, source.inc := 'Prevalence']


ps <- old[source.inc=='Prevalence survey', unique(iso3)]
est[iso3 %in% ps, source.inc := 'Prevalence']




#' # Expert opinion
#'
bak <- copy(est)
est[is.na(source.inc), source.inc := 'Expert opinion']


#' add last point based on current inc trends
#'
(table(est$year[is.na(est$inc)])) # should only be missing values where year==yr
lst <- as.character(unique(est$iso3[is.na(est$inc)]))

est[, chg.inc := c(NA, diff(log(inc))), by = iso3]
est[is.infinite(chg.inc), chg.inc := NA]
est[is.nan(chg.inc), chg.inc := NA]

est[iso3 %in% lst, chg.inc := na_interpolation(chg.inc), by =
      iso3]
est[iso3 %in% lst, shift.inc := shift(inc, n = 1L, fill = NA, type = 'lag'), by =
      iso3]
est[iso3 %in% lst, shift.inc.sd := shift(inc.sd,
                                         n = 1L,
                                         fill = NA,
                                         type = 'lag'), by = iso3]

est[iso3 %in% lst &
      year == yr, inc := exp(log(shift.inc) + chg.inc)]
est[iso3 %in% lst &
      year == yr, inc.sd := exp(log(shift.inc.sd) + chg.inc)]
est[iso3 %in% lst, chg.inc := c(NA, diff(log(inc))), by = iso3]


#' inc increasing in 2019
#' 
sel <- est$chg.inc > 0 & !is.na(est$chg.inc) & est$year == yr
incr19 <- as.character(unique(est$iso3[sel]))

#' increasing in 2018
#' 
sel <- est$chg.inc > 0 & !is.na(est$chg.inc) & est$year == yr - 1
incr18 <- as.character(unique(est$iso3[sel]))

incr <- intersect(lst, union(incr19, incr18))


#' reset inc where a recent increase is unlikely
#'
retrend <- function(x,
                    country,
                    ystart = 2000,
                    trend = F) {
  sel <- x$iso3 == country
  sel2 <- sel & x$year >= ystart
  sel3 <- sel & x$year == ystart - 1
  
  ninc <- x$inc[sel2][1]
  ninc.sd <- x$inc.sd[sel2][1]
  slope <- ifelse(trend, x$chg.inc[sel3], 0)
  
  x$inc[sel2] <- exp(log(ninc) + slope * (0:(dim(x[sel2])[1] - 1)))
  x$inc.sd[sel2] <-
    exp(log(ninc.sd) + slope * (0:(dim(x[sel2])[1] - 1)))
  return(x)
}

bak <- copy(est)
est <- retrend(est, 'GNQ', 2016, trend = F)

est[iso3 == 'LBY' & year == yr, inc := inc / 0.68]
est[iso3 == 'LBY' & year == yr, inc.sd := inc * .22]


# val <-
#   ihme['GNB'][measure_name == 'Incidence' &
#                 cause_name == 'Tuberculosis' &
#                 metric_name == 'Rate'][, val]
# est[iso3 == 'GNB', inc := na.kalman(c(val, NA))] # from IHME
# est[iso3 == 'GNB', inc.sd := inc * 0.2]



#' changes in inc based on chg.newinc
#'
retrend2 <- function(x, country, ystart = 2000) {
  sel <- x$iso3 == country
  sel2 <- sel & x$year >= ystart
  cdr <- x$newinc[sel2][1] / x$inc[sel2][1]
  rsd <- x$inc.sd[sel2][1] / x$inc[sel2][1]
  x$inc[sel2] <- x$newinc[sel2] / cdr
  x$inc[sel] <- predict(loess(x$inc[sel] ~ x$year[sel]), span = 0.6)
  x$inc.sd[sel2] <- x$inc[sel2] * rsd
  x$chg.newinc[sel] <- c(NA, diff(log(x$inc[sel])))
  return(x)
}

est <- retrend2(est, 'KHM', 2016)


#' bug fixes from last year
#'
est[iso3 == 'KHM' & year ==2018, inc := NA]
est[iso3 == 'KHM', inc := na_interpolation(inc)]
est[iso3 == 'KHM', inc.sd := inc * .2]

#' SDN bug fix
#' 
est[iso3 == 'SDN' & year < 2011, inc := newinc / .5]
est[iso3 == 'SDN' & year < 2011, inc.sd := inc * .27]

#' MWI
#' 
est[iso3 == 'MWI' & year > 2005, inc := NA]
est[iso3 == 'MWI' & year > 2005, inc.sd := NA]
est[iso3 == 'MWI' & year==2013, inc := 261.4]
est[iso3 == 'MWI' & year==2015, inc := newinc/.477]
est[iso3 == 'MWI' & year==2013, inc.sd := 71.82]
est[iso3 == 'MWI' & year==2019, inc := newinc/.621]
est[iso3 == 'MWI', inc := na_kalman(inc)]
est[iso3 == 'MWI', inc.sd := inc * 0.2748]
# 
# (est['MWI',.(year,newinc,inc,inc.sd, newinc/inc)])
# iplot('mwi')


#' checks
#'
est[, test.ispos(inc)]
est[, test.ispos(inc.sd)]


# sel <- est$iso3 %ni% c('KAZ', 'RUS')
# est[sel & year == yr][, test.AgeB(inc, imp.newinc)]
# (est[sel &
#                  inc < imp.newinc, .(iso3, year, inc, inc.sd, newinc, imp.newinc)])
# 
# #+ fig.width=8, fig.height=6
# print(iplot('ERI'))
# print(iplot('GEO'))
# print(iplot('UKR'))


#' adj 2019 ---
#' 
est[iso3 == 'SSD', inc := newinc[9] / .65]
est[iso3 == 'SSD', inc.sd := inc * .2]
est <- retrend2(est, 'GMB', 2013)
est[iso3 == 'TGO' & year==yr-1, inc := newinc/.81]
est[iso3 == 'TGO' & year==yr, inc := newinc/.86]
est[iso3 == 'TGO' & year >= (yr - 1), inc.sd := inc*.1025]
est[iso3 == 'PER' & year >= (yr - 1), inc := NA]
est[iso3 == 'PER', inc := na_interpolation(inc)]
est[iso3 == 'PER' & year >= yr - 1, inc.sd := inc  * .1276]


#' update change in inc
#'
est[, chg.inc := c(NA, diff(log(inc))), by = iso3]
est[, shift.inc := NULL]
est[, shift.inc.sd := NULL]



#' # Comparison plots with last year's global TB report
#'
#+ fig.width=16, fig.height=12
for (i in wr) {
  p <- qplot(
    year,
    inc,
    data = subset(est, g.whoregion == i),
    geom = 'line',
    colour = I('grey90')
  ) +
    geom_ribbon(
      aes(
        year,
        ymin = inc - 1.96 * inc.sd,
        ymax = inc + 1.96 * inc.sd
      ),
      fill = I('blue'),
      alpha = I(.4)
    ) +
    geom_line(
      aes(year, inc),
      data = subset(old, g.whoregion == i),
      colour = I('red'),
      linetype = I(2)
    ) +
    geom_line(aes(year, newinc)) +
    facet_wrap(~ iso3, scales = 'free_y') + xlab('') + ylab('Incidence rate per 100k/yr')
  suppressWarnings(print(p))
  suppressWarnings(ggsave(here(
    paste('output/checks/inc', i, '_compare.pdf', sep = '')),
    width = 14,
    height = 8
  ))
}



#' # Comparison plots with focus on recent trends
#'
#+ fig.width=16, fig.height=12
for (i in wr) {
  p <- qplot(
    year,
    inc,
    data = subset(est, g.whoregion == i & year > 2013),
    geom = 'line',
    colour = I('grey90')
  ) +
    geom_ribbon(
      aes(
        year,
        ymin = inc - 1.96 * inc.sd,
        ymax = inc + 1.96 * inc.sd
      ),
      fill = I('blue'),
      alpha = I(.4)
    ) +
    geom_line(
      aes(year, inc),
      data = subset(old, g.whoregion == i & year > 2013),
      colour = I('red'),
      linetype = I(2)
    ) +
    geom_line(aes(year, newinc)) +
    facet_wrap(~ iso3, scales = 'free_y') + xlab('') + ylab('Incidence rate per 100k/yr')
  suppressWarnings(print(p))
  suppressWarnings(ggsave(here(
    paste('output/checks/inc', i, '_compare_last5years.pdf', sep = '')),
    width = 14,
    height = 8
  ))
}



#' # Indirect estimation of HIV prev in TB
#'
#' the rest of missing values is imputed using
#' a prediction model based on UNAIDS estimates of the
#' prevalence of HIV in the general population
#'
bak <- copy(est)
(dim(est))
est <-
  merge(
    est,
    unaids[, .(
      iso3,
      year,
      hiv.num,
      hiv.lo.num,
      hiv.hi.num,
      mort.hiv.num,
      mort.hiv.lo.num,
      mort.hiv.hi.num
    )],
    by = c('iso3', 'year'),
    all.x = T,
    all.y = F
  )
(dim(est))

est[, hiv := hiv.num / pop]
est[, hiv.lo := hiv.lo.num / pop]
est[, hiv.hi := hiv.hi.num / pop]
est[, hiv.sd := (hiv.hi - hiv.lo) / 3.92]

est[, mort.hiv := mort.hiv.num / pop * m]
est[, mort.hiv.lo := mort.hiv.lo.num / pop * m]
est[, mort.hiv.hi := mort.hiv.hi.num / pop * m]
est[, mort.hiv.sd := (mort.hiv.hi - mort.hiv.lo) / 3.92]



#' # IRR
#'
set.seed(121)

#' inc.h
#'
inc.h <- with(est, prodXY(inc, tbhiv, inc.sd ^ 2, tbhiv.sd ^ 2))
est$inc.h <- inc.h[[1]]
est$inc.h.sd <- sqrt(inc.h[[2]])

#' inc.nh
#'
inc.nh <-
  with(est, prodXY(inc, (1 - tbhiv), inc.sd ^ 2, tbhiv.sd ^ 2))
est$inc.nh <- inc.nh[[1]]
est$inc.nh.sd <- sqrt(inc.nh[[2]])

#' force of infection in HIV+
#'
fi.h <-
  with(est, divXY(inc.h / m, hiv, (inc.h.sd / m) ^ 2, hiv.sd ^ 2))
est$fi.h <- fi.h[[1]]
est$fi.h.sd <- sqrt(fi.h[[2]])

#' force of infection in HIV-
#'
fi.nh <-
  with(est, divXY(inc.nh / m, (1 - hiv), (inc.nh.sd / m) ^ 2, hiv.sd ^ 2))
est$fi.nh <- fi.nh[[1]]
est$fi.nh.sd <- sqrt(fi.nh[[2]])

#' incidence rate ratio, ignoring covariance
#'
irr <- with(est, divXY(fi.h, fi.nh, fi.h.sd ^ 2, fi.nh.sd ^ 2))
est$irr <- irr[[1]]
est$irr.sd <- sqrt(irr[[2]])
sel <- !is.na(est$irr) & est$irr > 1e3
table(sel)
est[sel, .(iso3, year, inc, tbhiv, inc.nh, inc.h, hiv, fi.h, fi.nh, irr)]
est$irr[sel] <- est$irr.sd[sel] <- NA

sel <-
  est$irr.sd == 0 & est$hiv > 0 &
  !is.na(est$irr) & !is.na(est$irr.sd)
table(sel)
est$irr[sel] <- NA
est$irr.sd[sel] <- NA
lst <- unique(as.character(est$iso3[sel]))

est[iso3 %in% lst, irr := na.interpolation(irr), by = iso3]
est[iso3 %in% lst, irr.sd := na.interpolation(irr.sd), by = iso3]


#' check IRR series
#'
#+ fig.width=16, fig.height=12
# for (i in wr) {
#   p <-
#     qplot(
#       year,
#       0,
#       data = subset(est, g.whoregion == i),
#       geom = 'line',
#       colour = I('grey90')
#     ) +
#     geom_line(aes(year, irr)) +
#     geom_ribbon(
#       aes(
#         year,
#         ymin = pmax(irr - 1.96 * irr.sd, 0),
#         ymax = irr + 1.96 * irr.sd
#       ),
#       fill = I('blue'),
#       alpha = I(.4)
#     ) +
#     facet_wrap( ~ iso3, scales = 'free_y') + xlab('') + ylab('IRR')
#   suppressWarnings(print(p))
# }




#' # tbhiv based on IRR
#'
(est[, .(median(irr, na.rm = T), mean(irr, na.rm = T)), by = year])
#' use last non-missing r = irr backwards in time
#'
#' $t = \frac{h r}{1 + h (r - 1)}$
#'
bak2 <- copy(est)
lcty <- unique(as.character(est$iso3))

est[!is.na(tbhiv), test.isbinom(tbhiv)]
est[!is.na(irr), test.ispos(irr)]



#' no irr, impute it
#'
est[, ghiv := hiv > 0.1]
est[is.na(ghiv), ghiv := F]
est[, hincome := g.income == 'HIC']
est[is.na(hincome), hincome := T]

out <-
  est[, .(n = .N, irr = weighted.mean(irr, w = pop, na.rm = TRUE)),
      by = list(year, ghiv, hincome)]

# out <-
#   est[, .(n = .N, irr = weighted.mean((inc * tbhiv / hiv) / (inc * (1 - tbhiv) /
#                                                                (1 - hiv)), w = pop, na.rm = TRUE),
#           irr.sd = sd(irr)),
#       by = list(year, ghiv, hincome)]
(out)

out[is.infinite(irr), irr := NA]
out[is.nan(irr), irr := NA]
# out[24, irr := NA]
out[, irr := na_interpolation(irr), by = list(ghiv, hincome)]

# mirr <- out[ghiv == F & hincome == T, mean(irr, na.rm = T)]
# out[is.na(irr), irr := mirr]
(out)

#' use g.income and generalized HIV as predictors of IRR
#'
(dim(est))
est <- merge(est,
             out[, .(year, ghiv, hincome, e.irr = irr)],
             by = c('year', 'ghiv', 'hincome'),
             all.x = T)
(dim(est))
setkey(est, iso3, year)


#' one IRR available
#' 
est[, n.irr := sum(!is.na(irr) & irr > 0), by=iso3]
est[n.irr==1, f.irr := irr / e.irr]
est[n.irr==1, f.irr.sd := irr.sd / irr]
est[n.irr==1, f.irr := max(f.irr, na.rm = TRUE), by = iso3]
est[n.irr==1, f.irr.sd := max(f.irr.sd, na.rm = TRUE), by = iso3]
est[n.irr==1, irr := e.irr * f.irr]
est[n.irr==1, irr.sd := irr * f.irr.sd]


#' multiple IRR available (todo: fix this)
#' 
est[n.irr>1, f.irr := mean(irr, na.rm=T) / mean(e.irr, na.rm=T)]
est[n.irr>1, f.irr.sd := mean(irr.sd, na.rm=T) / mean(irr, na.rm=T)]
est[n.irr>1, f.irr := max(f.irr, na.rm = TRUE), by = iso3]
est[n.irr>1, f.irr.sd := max(f.irr.sd, na.rm = TRUE), by = iso3]
est[n.irr>1, irr := e.irr * f.irr]
est[n.irr>1, irr.sd := irr * f.irr.sd]

#' no IRR available
#' 
est[n.irr==0, irr := e.irr]
est[n.irr==0, irr.sd := e.irr * 0.25]


est[!is.na(irr), test.ispos(irr)]

est[, n.irr := NULL]
est[, f.irr := NULL]
est[, f.irr.sd := NULL]
est[is.infinite(irr.sd), irr.sd := irr * .25]
est[is.na(hiv.sd) & !is.na(hiv), hiv.sd := hiv*.25]

#' e.tbhiv using imputed IRR
#'
sel <- !is.na(est$irr) & is.na(est$irr.sd)
table(sel)
# est$irr.sd[sel] <- .25 * est$irr[sel] # arbitrary
# est[irr.sd > irr * .25, irr.sd := irr * .25]   # as well

excl <- c('BLZ')
sel <-
  !is.na(est$irr) & !is.na(est$hiv) & est$hiv > 0 &
  est$iso3 %ni% excl
table(sel)

out <- est[sel, {
  tmp = h2t(hiv, hiv.sd, irr, irr.sd)$prop
  
  list(e.tbhiv = tmp[2],
       e.tbhiv.sd = tmp[4])
},
by = .(iso3, year)]


# out[is.infinite(e.tbhiv), e.tbhiv := NA]
# out[is.infinite(e.tbhiv.sd), e.tbhiv.sd := NA]

# smooth e.tbhiv series
out[, e.tbhiv := predict(loess(e.tbhiv~year, span=0.6)), by=iso3]
out[, e.tbhiv.sd := predict(loess(e.tbhiv.sd~year, span=0.6)), by=iso3]
out[e.tbhiv<0, e.tbhiv := 0]
out[e.tbhiv.sd<0, e.tbhiv.sd := 0]

out[!is.na(e.tbhiv), test.isbinom(e.tbhiv)]
out[!is.na(e.tbhiv), test.isbinom(e.tbhiv.sd)]

est$e.tbhiv[sel] <- out$e.tbhiv
est$e.tbhiv.sd[sel] <- out$e.tbhiv.sd





#' impute missing tbhiv with e.tbhiv
#'
sel <- is.na(est$tbhiv) & !is.na(est$e.tbhiv)
table(sel)
est$tbhiv[sel] <- est$e.tbhiv[sel]
est$tbhiv.sd[sel] <- est$e.tbhiv.sd[sel]

est[!is.na(tbhiv), test.isbinom(tbhiv)]
est[!is.na(tbhiv.sd), test.isbinom(tbhiv.sd)]

sum(is.na(est$tbhiv))
sum(is.na(est$tbhiv.sd))




#' fixes
#'
sel <- est$iso3 == 'KHM' & est$year %in% 2000:2009
sel2 <- old$iso3 == 'KHM' & old$year %in% 2000:2009
est$tbhiv[sel] <- old$tbhiv[sel2]
est$tbhiv.sd[sel] <- old$tbhiv.sd[sel2]

incl <- c('ZAF', 'ZMB', 'ZWE')
sel <-
  est$year < 2010 & is.na(est$tbhiv.routine.ok) & est$iso3 %in% incl
table(sel)
est$tbhiv[sel] <- est$e.tbhiv[sel]


est$source.inc[est$source.inc=='Prevalence'] <- 'Prevalence survey'
est$source.inc[est$iso3=='BGD'] <- 'Prevalence survey'
est$source.inc[est$source.inc %in% c('Standard adjustment')] <- 'Case notifications,\nStandard adjustment'
est$source.inc[est$source.inc=='Expert opinion'] <- 'Case notifications,\nExpert opinion'
est$source.inc[est$source.inc=='Inventory'] <- 'Inventory study'







#' check imputed series
#'

#+ fig.width=16, fig.height=12
for (i in wr) {
  p <-
    qplot(
      year,
      0,
      data = subset(est, g.whoregion == i),
      geom = 'line',
      colour = I('grey90')
    ) +
    geom_line(aes(year, tbhiv)) +
    geom_ribbon(
      aes(
        year,
        ymin = tbhiv - 1.96 * tbhiv.sd,
        ymax = tbhiv + 1.96 * tbhiv.sd
      ),
      fill = I('blue'),
      alpha = I(.4)
    ) +
    geom_point(aes(year, tbhiv.routine),
               colour = I('black'),
               shape = I(4)) +
    geom_point(aes(year, tbhiv.routine.ok),
               colour = I('blue'),
               shape = I(4)) +
    geom_point(aes(year, tbhiv.surv),
               colour = I('green'),
               shape = I(2)) +
    geom_point(aes(year, tbhiv.sentin),
               colour = I('red'),
               shape = I(3)) +
    facet_wrap(~ iso3, scales = 'free_y') + xlab('') + ylab('HIV prevalence in TB')
  suppressWarnings(print(p))
}

#' compare with last year
#'

#+ fig.width=16, fig.height=12
for (i in wr) {
  p <-
    qplot(
      year,
      0,
      data = subset(est, g.whoregion == i),
      geom = 'line',
      colour = I('grey90')
    ) +
    geom_ribbon(
      aes(
        year,
        ymin = tbhiv - 1.96 * tbhiv.sd,
        ymax = tbhiv + 1.96 * tbhiv.sd
      ),
      fill = I('blue'),
      alpha = I(.4)
    ) +
    geom_line(
      aes(year, tbhiv),
      data = subset(old, g.whoregion == i),
      colour = I('red'),
      linetype = I(2)
    ) +
    facet_wrap(~ iso3, scales = 'free_y') + xlab('') + ylab('HIV prevalence in TB')
  suppressWarnings(print(p))
}




#' global TBHIV
#'
(est[, weighted.mean(tbhiv, w = inc * pop / 1e5, na.rm = T), by = year])
#' # IRR
#'
set.seed(121)

#' inc.h
#'
inc.h <- with(est, prodXY(inc, tbhiv, inc.sd ^ 2, tbhiv.sd ^ 2))
est$inc.h <- inc.h[[1]]
est$inc.h.sd <- sqrt(inc.h[[2]])

#' inc.nh
#'
inc.nh <-
  with(est, prodXY(inc, (1 - tbhiv), inc.sd ^ 2, tbhiv.sd ^ 2))
est$inc.nh <- inc.nh[[1]]
est$inc.nh.sd <- sqrt(inc.nh[[2]])

#' force of infection in HIV+
#'
fi.h <-
  with(est, divXY(inc.h / m, hiv, (inc.h.sd / m) ^ 2, hiv.sd ^ 2))
est$fi.h <- fi.h[[1]]
est$fi.h.sd <- sqrt(fi.h[[2]])

#' force of infection in HIV-
#'
fi.nh <-
  with(est, divXY(inc.nh / m, (1 - hiv), (inc.nh.sd / m) ^ 2, hiv.sd ^ 2))
est$fi.nh <- fi.nh[[1]]
est$fi.nh.sd <- sqrt(fi.nh[[2]])

#' incidence rate ratio, ignoring covariance
#'
irr <- with(est, divXY(fi.h, fi.nh, fi.h.sd ^ 2, fi.nh.sd ^ 2))
est$irr <- irr[[1]]
est$irr.sd <- sqrt(irr[[2]])
sel <- !is.na(est$irr) & est$irr > 1e3
table(sel)

(est[sel, .(iso3, year, inc, tbhiv, inc.nh, inc.h, hiv, fi.h, fi.nh, irr)])
est$irr[sel] <- est$irr.sd[sel] <- NA




#' bounds
#'
sel1 <- est$inc > 0
table(sel1)
out1 <- vlohi(est$inc[sel1] / m, est$inc.sd[sel1] / m)

sel2 <- est$inc.nh > 0 & !is.na(est$inc.nh)
table(sel2)
out2 <- vlohi(est$inc.nh[sel2] / m, est$inc.nh.sd[sel2] / m)

sel3 <-
  est$inc.h > 0 &
  !is.na(est$inc.h) & est$inc.h.sd > 0 & !is.na(est$inc.h.sd)
table(sel3)
out3 <- vlohi(est$inc.h[sel3] / m, est$inc.h.sd[sel3] / m)

est$inc.lo[sel1] <- out1[1,] * m
est$inc.hi[sel1] <- out1[2,] * m
est$inc.lo[!sel1] <- est$inc.hi[!sel1] <- est$inc[!sel1]

est$inc.nh.lo[sel2] <- out2[1,] * m
est$inc.nh.hi[sel2] <- out2[2,] * m
est$inc.nh.lo[!sel2 & est$inc.nh == 0 & est$inc.nh.sd == 0] <- 0
est$inc.nh.hi[!sel2 & est$inc.nh == 0 & est$inc.nh.sd == 0] <- 0

est$inc.h.lo[sel3] <- out3[1,] * m
est$inc.h.hi[sel3] <- out3[2,] * m
est$inc.h.lo[!sel3 & est$inc.h == 0 & est$inc.h.sd == 0] <- 0
est$inc.h.hi[!sel3 & est$inc.h == 0 & est$inc.h.sd == 0] <- 0

sel4 <-
  (est$inc.h.lo > est$inc.h) |
  (est$inc.h.hi < est$inc.h) &
  (!is.na(est$inc.h) & !is.na(est$inc.h.lo) & !is.na(est$inc.h.hi))
table(sel4)
est[sel4, .(iso3, year, inc, inc.h, inc.h.sd, inc.h.lo, inc.h.hi)]
est[sel4, inc.h.lo := 0]
est[sel4, inc.h.hi := inc.h + 1.96 * inc.h.sd]

sel <-
  (est$inc.nh.lo > est$inc.nh) |
  (est$inc.nh.hi < est$inc.nh) &
  (!is.na(est$inc.nh) &
     !is.na(est$inc.nh.lo) & !is.na(est$inc.nh.hi))
table(sel)
est[sel, .(iso3, year, inc, inc.nh, inc.nh.sd, inc.nh.lo, inc.nh.hi)]
# est[sel, inc.nh.hi := inc.nh + 1.96 * inc.nh.sd]

est[is.na(inc.nh.lo), inc.nh.lo := 0]
est[is.na(inc.nh.hi), inc.nh.hi := inc.nh + 1.96 * inc.nh.sd]

sel <- !is.na(est$inc.h)
est[sel & is.na(inc.h.lo), inc.h.lo := 0]
est[sel & is.na(inc.h.hi), inc.h.hi := inc.h + 1.96 * inc.h.sd]


sel <- est$tbhiv > 0 & est$tbhiv < 1 & !is.na(est$tbhiv)
table(sel)
out <- vlohi(est$tbhiv[sel], est$tbhiv.sd[sel])
est$tbhiv.lo[sel] <- out[1,]
est$tbhiv.hi[sel] <- out[2,]


sel <- est$tbhiv == 0 & !is.na(est$tbhiv)
table(sel)
est$tbhiv.lo[sel] <- 0
est$tbhiv.hi[sel] <- est$tbhiv.sd[sel] * 1.96

sel <- est$tbhiv == 1 & !is.na(est$tbhiv)
table(sel)
est$tbhiv.hi[sel] <- 1
est$tbhiv.lo[sel] <- pmax(1 - est$tbhiv.sd[sel] * 1.96, 0)

est[tbhiv.hi<tbhiv,.(iso3,tbhiv,tbhiv.sd,tbhiv.lo,tbhiv.hi)]
est[tbhiv.hi<tbhiv, tbhiv.hi := tbhiv + 1.96 * tbhiv.sd]

#' checks
#'
est[, test.isbinom(inc / m)]
est[!is.na(inc.nh), test.isbinom(inc.nh / m)]
est[!is.na(tbhiv), test.isbinom(tbhiv)]
est[!is.na(inc.h), test.isbinom(inc.h / m)]

est[, test.bounds(inc, inc.lo, inc.hi)]
est[!is.na(inc.nh), test.bounds(inc.nh, inc.nh.lo, inc.nh.hi)]
est[!is.na(inc.h), test.bounds(inc.h, inc.h.lo, inc.h.hi)]
est[!is.na(tbhiv), test.bounds(tbhiv, tbhiv.lo, tbhiv.hi)]

est[!is.na(inc.h), sum(abs(inc.h + inc.nh - inc) > 1) == 0]
est[!is.na(inc.h), test.ispos(inc.h)]
est[!is.na(inc.h), test.ispos(inc.h.sd)]
est[!is.na(inc.h), test.ispos(inc.nh)]
est[!is.na(inc.h), test.ispos(inc.nh.sd)]
est[!is.na(irr), test.ispos(irr)]



#' # Comparison plots with focus on recent trends, HIV+ incidence
#'
#+ fig.width=16, fig.height=12
for (i in wr) {
  p <- qplot(
    year,
    inc.h,
    data = subset(est, g.whoregion == i & year > 2013),
    geom = 'line',
    colour = I('grey90')
  ) +
    geom_ribbon(
      aes(
        year,
        ymin = inc.h - 1.96 * inc.h.sd,
        ymax = inc.h + 1.96 * inc.h.sd
      ),
      fill = I('blue'),
      alpha = I(.4)
    ) +
    geom_line(
      aes(year, inc.h),
      data = subset(old, g.whoregion == i & year > 2013),
      colour = I('red'),
      linetype = I(2)
    ) +
    facet_wrap( ~ iso3, scales = 'free_y') + xlab('') + ylab('Incidence rate per 100k/yr')
  suppressWarnings(print(p))
}


#' # Comparison plots with last year's global TB report
#'
#+ fig.width=16, fig.height=12
for (i in wr) {
  p <- qplot(
    year,
    hiv,
    data = subset(est, g.whoregion == i),
    geom = 'line',
    colour = I('grey90')
  ) +
    geom_ribbon(
      aes(
        year,
        ymin = hiv - 1.96 * hiv.sd,
        ymax = hiv + 1.96 * hiv.sd
      ),
      fill = I('blue'),
      alpha = I(.4)
    ) +
    geom_line(
      aes(year, hiv),
      data = subset(old, g.whoregion == i),
      colour = I('red'),
      linetype = I(2)
    ) +
    facet_wrap(~ iso3, scales = 'free_y') + xlab('') + ylab('HIV prevalence')
  suppressWarnings(print(p))
  suppressWarnings(ggsave(here(
    paste('output/checks/hiv', i, '_compare.pdf', sep = '')),
    width = 14,
    height = 8
  ))
}



#' # Comparison plots with last year's global TB report
#'
#+ fig.width=16, fig.height=12
for (i in wr) {
  p <- qplot(
    year,
    tbhiv,
    data = subset(est, g.whoregion == i),
    geom = 'line',
    colour = I('grey90')
  ) +
    geom_ribbon(
      aes(
        year,
        ymin = tbhiv - 1.96 * tbhiv.sd,
        ymax = tbhiv + 1.96 * tbhiv.sd
      ),
      fill = I('blue'),
      alpha = I(.4)
    ) +
    geom_line(
      aes(year, tbhiv),
      data = subset(old, g.whoregion == i),
      colour = I('red'),
      linetype = I(2)
    ) +
    facet_wrap(~ iso3, scales = 'free_y') + xlab('') + ylab('HIV prevalence in TB')
  suppressWarnings(print(p))
  suppressWarnings(ggsave(here(
    paste('output/checks/tbhiv', i, '_compare.pdf', sep = '')),
    width = 14,
    height = 8
  ))
}



#' global aggregates
#'
(est[, .(inc.num = as.integer(sum(inc * pop / 1e5))), by = year])
(est[, .(inc.h.num = as.integer(sums(inc.h * pop / 1e5))), by =
                 year])



#' save
#'
save(est, file = here('Rdata/est.Rdata'))
fwrite(est, file = here(paste0('csv/est_s03_', Sys.Date(), '.csv')))
