#' ---
#' title: Mortality HIV-neg
#' author: Philippe Glaziou
#' date: 07/07/2020
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
#' Mortality HIV-neg
#'
#' # Load libraries and data
#'
library(data.table)
library(imputeTS)
library(propagate)
library(here)


load(here('Rdata/cty.Rdata'))
load(here('Rdata/est.Rdata'))
load(here('Rdata/ihme.Rdata'))
load(here('Rdata/ihme.tb.Rdata'))
load(here('Rdata/vr.Rdata'))
load(here('Rdata/old.Rdata'))
load(here('Rdata/gbd.Rdata'))

m <- 1e5
yr <- 2019

source(here('R/fun.R'))


#' vectorized lohi
#'
vlohi <- Vectorize(lohi, c('ev', 'sd'))



#' # Indirect mort.nh
#'
est$t <- est$tbhiv
est$t[is.na(est$tbhiv)] <- 0

est$U <- 0  # under-reporting of detected/treated cases
# est$U[est$iso3=='PAK'] <- .27
est$U[est$iso3 == 'IRQ'] <- .16
est$U[est$iso3 == 'IDN'] <- c(rep(.5, 17), .3, .12, .12)
est$U[est$iso3 == 'IND'] <- c(rep(.3, 18), .22, .1)

setkey(est, iso3, year)
out <-
  est[, {
    tmp = inc2mort(inc, inc.sd, imp.newinc / (1 - U), tbhiv, tbhiv.sd, noHIV =
                     T)$prop
    
    list(e.mort.nh = tmp[2],
         e.mort.nh.sd = tmp[4])
  },
  by = .(iso3, year)] # this takes a while

dim(est)
dim(out)
est <- merge(est, out, by = c('iso3', 'year'), all.x = T)
dim(est)


#' Import vr
est2 <-
  merge(est, vr[, .(
    iso3,
    year,
    vr.keep = keep.vr,
    vr.garbage = garbage,
    vr.coverage,
    vr.quality = codqual,
    vr.env = env,
    #                           ghe.env, ghe.env.lo, ghe.env.hi,
    vr.mort.nh = tb.adj * m / pop,
    vr.raw = tb * m / pop,
    vr.mort.nh.sd = tb.adj.sd * m / pop
  )],
  by = c('iso3', 'year'), all.x = TRUE)
dim(est)
dim(est2)
dim(vr)


#' import IHME estimates adjusted for WHO/IHME env ratio
#'
est3 <- merge(est2, gbd, by = c('iso3', 'year'), all.x = T)
dim(est2)
dim(est3)
est4 <- merge(est3,
              ihme.tb[, .(
                iso3,
                year,
                ihme.mort.nh.num = val,
                ihme.mort.nh.num.sd = (upper - lower) / 3.92
              )],
              by = c('iso3', 'year'),
              all.x = T)
dim(est3)
dim(est4)
est4[, ihme.mort.nh.num := ihme.mort.nh.num * env.ratio]
est4[, ihme.mort.nh.num.sd := ihme.mort.nh.num.sd * env.ratio]





#' # Set mort.nh to VR, IHME adjusted, or predicted
#'
#' VR
#'
est4[iso3 %in% c('ZAF', 'MNG') & !is.na(est4$vr.keep), vr.keep := FALSE]

sel <-
  est4$vr.keep == T & !is.na(est4$vr.keep) & !is.na(est4$vr.mort.nh)
table(sel) # use valid VR values

est4[sel, mort.nh := vr.mort.nh]
est4[sel, mort.nh.sd := vr.mort.nh.sd]
(vr.lst <-
    unique(as.character(est4$iso3[sel]))) # 118 countries VR based



#' impute missing VR in the series
#'
sel2 <- is.na(est4$mort.nh) & est4$iso3 %in% vr.lst
(vr.imp <- unique(as.character(est4$iso3[sel2])))

sel <- est4$iso3 %in% vr.imp
est4[sel, nna := sum(!is.na(mort.nh)), by = iso3] # nna non-missing values in each series
table(est4$nna)


#' only one not missing -> LOCF
#'
est4[nna == 1, .(iso3, year, mort.nh)]
est4[nna == 1, mort.nh := imputeTS::na_locf(mort.nh, na_remaining = "keep"), by = iso3]
est4[nna == 1, mort.nh.sd := imputeTS::na_locf(mort.nh.sd, na_remaining = "keep"), by = iso3]

est4[sel, nonzero := sum(mort.nh > 0 & !is.na(mort.nh)), by = iso3]
table(est4$nonzero)

#' more than one not missing, all non-missing are zeros -> LOCF
#'
est4[nna > 1 & nonzero == 0, .(iso3, year, mort.nh)]
est4[nna > 1 & nonzero == 0, mort.nh := imputeTS::na_locf(mort.nh, na_remaining = "keep"), by = iso3]
est4[nna > 1 &
       nonzero == 0, mort.nh.sd := imputeTS::na_locf(mort.nh.sd, na_remaining = 'keep'), by = iso3]


#' exactly 2 are not missing, at least one is not a zero -> simple interpolation
#'
est4[nna == 2 &
       nonzero > 0, mort.nh := imputeTS::na_interpolation(mort.nh), by = iso3]
est4[nna == 2 &
       nonzero > 0, mort.nh.sd := imputeTS::na_interpolation(mort.nh.sd), by = iso3]

nokalman <-
  c(
    'SYR',
    'KWT',
    'MKD',
    'MNE',
    'LKA',
    'TJK',
    'BLR',
    'BGR',
    'EGY',
    'EST',
    'HND',
    'ZAF',
    'AZE',
    'CYP',
    'JPN',
    'THA'
  ) # use simple interpolation to avoid predicting implausible trends
est4[sel, mm := mean(mort.nh, na.rm = T), by = iso3]

#' very low average mort.nh (<1) or ISO3 in nokalman list -> simple interpolation
#'
nokalman2 <-
  unique(as.character(est4$iso3[est4$mm < 1 &
                                  !is.na(est4$mm)])) # simple interpolation in low mortalites
est4[sel &
       iso3 %in% union(nokalman, nokalman2), mort.nh := imputeTS::na_interpolation(mort.nh), by =
       iso3]
est4[sel &
       iso3 %in% union(nokalman, nokalman2), mort.nh.sd := imputeTS::na_interpolation(mort.nh.sd), by =
       iso3]

#' use Kalman filter in the remaining series
#'
est4[sel & nna > 2, mort.nh := imputeTS::na_kalman(mort.nh), by = iso3]
est4[sel &
       nna > 2, mort.nh.sd := imputeTS::na_interpolation(mort.nh.sd), by = iso3]
summary(est4$mort.nh[sel & est4$nna > 2])

est4[sel, nna := sum(!is.na(mort.nh.sd)), by = iso3]
table(est4$nna[sel])

est4[sel &
       nna > 1, mort.nh.sd := imputeTS::na_interpolation(mort.nh.sd), by = iso3]
sel2 <- sel & is.na(est4$mort.nh.sd)
table(sel2) # should be all F

(est4['GBR', .(iso3, year, mort.nh, mort.nh.sd, vr.mort.nh, vr.mort.nh.sd)])


#' fix TJK long flat trend, assuming cst CFR
#'
sel <- est4$iso3 == 'TJK' & est4$year %in% 2006:2016
est4$mort.nh[sel] <- est4$inc[sel] * 0.07
est4$mort.nh.sd[sel] <-
  est4$mort.nh.sd[sel] * est4$mort.nh[sel] / 12.816
est4['TJK', .(iso3, year, inc, mort.nh, cfr = mort.nh / inc, mort.nh.sd)]


#+ fig.width=16, fig.height=12
for (i in 1:length(vr.imp)) {
  pdf(here(paste('output/checks/imputations/imp_vr_', vr.imp[i], '.pdf', sep = '')))
  plotNA.imputations(est4[vr.imp[i], vr.mort.nh], est4[vr.imp[i], mort.nh])
  dev.off()
}

#' clean-up
#'
est4[, nna := NULL]
est4[, nonzero := NULL]
est4[, mm := NULL]
est4$source.mort[!is.na(est4$mort.nh) &
                   est4$iso3 %in% vr.imp] <- 'VR'

est4[iso3=='RUS', source.mort := 'VR']



#' IHME
#'
#' include vr.quality==4 and poor vr.coverage or high garbage
#'
sel <- est4$vr.keep == F & !is.na(est4$vr.keep)
table(sel)
lst <- unique(as.character(est4$iso3[sel]))


#' CHN: add manually (communications CCDC June 2018)
#'
# vrchn <- data.table(
#   year = 2004:2019,
#   mort.nh = c(
#     5.4,
#     5.49,
#     5.09,
#     4.61,
#     4.15,
#     4.05,
#     3.87,
#     3.39,
#     3.04,
#     3.05,
#     2.79,
#     2.81,
#     2.79,
#     2.62,
#     NA,
#     NA
#   ),
#   mort.nh.raw = c(
#     5.4,
#     5.49,
#     4.18,
#     3.79,
#     3.41,
#     3.31,
#     3.16,
#     3.01,
#     2.99,
#     2.54,
#     2.32,
#     2.34,
#     2.32,
#     2.18,
#     NA,
#     NA
#   )
# )
# 
# vrchn[, adj := 1 - mort.nh.raw / mort.nh]
# vrchn[, mort.nh.sd := (mort.nh / 4) * (1 / (1 - adj) - 1)]
# (vrchn)

ochn <- old['CHN']
chn <- est4['CHN']
chn[, mort.nh := c(ochn$mort.nh, NA)]
chn[, mort.nh.sd := c(ochn$mort.nh.sd, NA)]
chn$mort.nh[19] <- NA
chn$mort.nh.sd[19] <- NA
chn$mort.nh <- imputeTS::na_kalman(chn$mort.nh)
chn$mort.nh.sd <- imputeTS::na_kalman(chn$mort.nh.sd)

(chn[, .(year, mort.nh, mort.nh.sd)])
(est4['CHN', .(year, mort.nh, mort.nh.sd)])
est4['CHN', mort.nh := chn$mort.nh]
est4['CHN', mort.nh.sd := chn$mort.nh.sd]



#' +ZAF: adjusts for HIV miscoding
#'
#' +IND: no VR available at WHO
#'
#' +IDN: more VR than at WHO, better analysis
#'
#' +IRQ: borderline VR, only one datapoint
#'
#' +PAK: consistent with e.mort.nh at the tail end, but more plausible lead
#'
#' +UKR: inconsistent envelopes
#' 
prevlist <-
  c('ZAF', 'IDN', 'IND', 'IRQ', 'PAK', 'KHM', 'KIR', 'GEO','UKR') # dropped VNN and added KHM following feedback in Aug

est4[, ihme.mort.nh := ihme.mort.nh.num * m / pop]
est4[, ihme.mort.nh.sd := ihme.mort.nh.num.sd * m / pop]


sel <- est4$iso3 %in% union(lst, prevlist)
table(sel)
sum(!is.na(est4$mort.nh[sel]))

est4$mort.nh[sel] <- est4$ihme.mort.nh[sel]
est4$mort.nh.sd[sel] <- est4$ihme.mort.nh.sd[sel]
est4$source.mort[sel] <- 'IHME'


#' impute missing values
#'
ihme.lst <-
  unique(as.character(est4$iso3[!is.na(est4$ihme.mort.nh)]))
sel2 <- sel & is.na(est4$mort.nh) & est4$iso3 %in% ihme.lst
ihme.imp <- unique(as.character(est4$iso3[sel2]))

est4[iso3 %in% ihme.imp, mort.nh := imputeTS::na_kalman(mort.nh), by = iso3]
est4[iso3 %in% ihme.imp, mort.nh.sd := imputeTS::na_interpolation(mort.nh.sd), by =
       iso3]
est4['PAK', .(iso3, year, mort.nh, mort.nh.sd, ihme.mort.nh, ihme.mort.nh.sd)]


#' plot imputations
#'
#+ fig.width=16, fig.height=12
for (i in 1:length(ihme.imp)) {
  pdf(here(paste('output/checks/imputations/imp_ihme_', ihme.imp[i], '.pdf', sep = '')))
  plotNA.imputations(est4[ihme.imp[i], ihme.mort.nh], est4[ihme.imp[i], mort.nh])
  dev.off()
}




#' # Indirect estimates
#'
sum(is.na(est4$mort.nh))
sum(is.na(est4$mort.nh.sd))

sel <-
  is.na(est4$mort.nh) |
  est4$iso3 %in% c('MAR')  # IHME gives an implausibly high CFR in MAR
table(sel)
sum(is.na(est4$e.mort.nh[sel]))

est4[sel, mort.nh := e.mort.nh]
est4[sel, mort.nh.sd := e.mort.nh.sd]
est4[sel, source.mort := 'Indirect']
est4[iso3 == 'CHN', source.mort := 'VR']
table(est4$source.mort)

sum(is.na(est4$mort.nh)) == 0
sum(is.na(est4$mort.nh.sd)) == 0 # should both be TRUE
sum(is.na(est4$source.mort)) == 0 # should be TRUE

est4[, test.isbinom(mort.nh / m)]
est4[, test.isbinom(mort.nh.sd / m)]


sel <- est4$iso3 == 'KHM'
sel2 <- old$iso3 == 'KHM'

est4$mort.nh[sel] <- imputeTS::na_kalman(c(old$mort.nh[sel2], NA))
est4$mort.nh.sd[sel] <- imputeTS::na_kalman(c(old$mort.nh.sd[sel2], NA))


sel <- est4$iso3 == 'VNM'
sel2 <- old$iso3 == 'VNM'

est4$mort.nh[sel] <- imputeTS::na_kalman(c(old$mort.nh[sel2], NA))
est4$mort.nh.sd[sel] <- imputeTS::na_kalman(c(old$mort.nh.sd[sel2], NA))




#' add bounds and counts
#'
sel <- est4$mort.nh > 0 & est4$mort.nh.sd > 0
table(sel)

out <- vlohi(est4$mort.nh[sel] / m, est4$mort.nh.sd[sel] / m)

est4$mort.nh.lo[sel] <- out[1,] * m
est4$mort.nh.hi[sel] <- out[2,] * m

sel <- est4$mort.nh.sd == 0 & !is.na(est4$mort.nh.sd)
table(sel)
est4$mort.nh.lo[sel] <- est4$mort.nh[sel]
est4$mort.nh.hi[sel] <- est4$mort.nh[sel]

sel <-
  (est4$mort.nh.lo > est4$mort.nh) |
  (est4$mort.nh.hi < est4$mort.nh) &
  (!is.na(est4$mort.nh) &
     !is.na(est4$mort.nh.lo) & !is.na(est4$mort.nh.hi))
table(sel)
est4[sel, .(iso3, year, inc, mort.nh, mort.nh.sd, mort.nh.lo, mort.nh.hi)]
est4[sel, mort.nh.hi := mort.nh + 1.96 * mort.nh.sd]
est4[mort.nh.lo == 0, .(iso3, year, mort.nh, mort.nh.sd, mort.nh.lo, mort.nh.hi)]
est4[mort.nh > 0, test.bounds(mort.nh, mort.nh.lo, mort.nh.hi)]

est4[is.na(mort.nh.lo), mort.nh.lo := 0]
est4[is.na(mort.nh.hi), mort.nh.hi := mort.nh + 1.96 * mort.nh.sd]


est4 <- within(est4, {
  mort.nh.num <- mort.nh * pop / m
  mort.nh.lo.num <- mort.nh.lo * pop / m
  mort.nh.hi.num <- mort.nh.hi * pop / m
  
  inc.num <- inc * pop / m
  inc.lo.num <- inc.lo * pop / m
  inc.hi.num <- inc.hi * pop / m
  inc.nh.num <- inc.nh * pop / m
  inc.nh.lo.num <- inc.nh.lo * pop / m
  inc.nh.hi.num <- inc.nh.hi * pop / m
  inc.h.num <- inc.h * pop / m
  inc.h.lo.num <- inc.h.lo * pop / m
  inc.h.hi.num <- inc.h.hi * pop / m
  
})


#' checks
#'
est4[, test.bounds(mort.nh, mort.nh.lo, mort.nh.hi)]

est4[, .(sums(mort.nh.num),
         sums(inc.num),
         sums(inc.nh.num),
         sums(mort.nh.num) / sums(inc.nh.num)), by = year]
old[, .(sums(mort.nh.num),
        sums(inc.num),
        sums(inc.nh.num),
        sums(mort.nh.num) / sums(inc.nh.num)), by = year]

wr <- unique(as.character(est4$g.whoregion))

#+ fig.width=16, fig.height=12
for (i in wr) {
  p <-
    qplot(
      year,
      mort.nh,
      data = subset(est4, g.whoregion == i),
      geom = 'line',
      colour = I('blue')
    ) +
    geom_ribbon(
      aes(year, ymin = mort.nh.lo, ymax = mort.nh.hi),
      fill = I('blue'),
      alpha = I(0.4)
    ) +
    geom_line(
      aes(year, mort.nh),
      data = subset(old, g.whoregion == i),
      colour = I('red'),
      linetype = I(2)
    ) +
    facet_wrap( ~ iso3, scales = 'free_y') + xlab('') + ylab('Rate per 100,000/year')
  suppressWarnings(print(p))
  
  suppressWarnings(ggsave(here(
    paste('output/checks/mort.nh_', i, '_compare.pdf', sep = '')),
    width = 14,
    height = 8
  ))
}

est2 <- copy(est4)
save(est2, file = here('Rdata/est2.Rdata'))
fwrite(est2, file = here(paste0('csv/est_s06_', Sys.Date(), '.csv')))
