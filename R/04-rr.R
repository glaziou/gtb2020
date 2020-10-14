#' ---
#' title: RR/MDR estimates
#' author: Philippe Glaziou
#' date: 2020-07-03
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
#'
#' # Preamble
#' (Last updated: `r Sys.Date()`)
#'
#' RR incidence estimates, Global TB Report 2019
#'
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(imputeTS))
suppressPackageStartupMessages(library(propagate))
suppressPackageStartupMessages(library(metafor))
suppressPackageStartupMessages(library(whomap))
library(here)

rm(list = ls())


load(here('Rdata/est.Rdata'))
load(here('Rdata/tb.Rdata'))
load(here('Rdata/cty.Rdata'))
load(here('Rdata/rel.Rdata'))  # dataset on relapses

#' most recent DRS data
#'
load(here('Rdata/dre.Rdata'))

#' last year's DR estimates
#'
load(here('Rdata/drh.Rdata'))

source(here('R/fun.R'))

vlohi <- Vectorize(lohi, c('ev', 'sd'))
yr <- 2019



#' add mdr regions
#'
dre <- merge(dre, est[year == yr, .(iso3, g.mdr)], by = 'iso3', all.x =
               T)
dre['ANT', g.mdr := 'AMR']
dre['SCG', g.mdr := 'EUR']

#' prop RR, surveillance >= 2017
#'
sel1 <- dre$source.new == 'Surveillance' & dre$year.new >= 2017
table(sel1)
dre[sel1, rr.new.Num := rr.new]
dre[sel1, rr.new.Den := r.rlt.new]
dre[sel1, prop.rr.new := rr.new / r.rlt.new, by = iso3]
dre[sel1, test.isbinom(prop.rr.new)]

sel2 <-
  dre$source.ret == 'Surveillance' &
  dre$year.ret >= 2017 & dre$r.rlt.ret > 0 & !is.na(dre$rr.ret)
table(sel2)
dre[sel2, rr.ret.Num := rr.ret]
dre[sel2, rr.ret.Den := r.rlt.ret]
dre[sel2, prop.rr.ret := rr.ret / r.rlt.ret, by = iso3]
dre[sel2, test.isbinom(prop.rr.ret)]



#' prop RR, surveillance < 2017
#'
sel3 <- dre$source.new == 'Surveillance' & dre$year.new < 2017
table(sel3)
dre[sel3, rr.new.Num := as.integer(rowSums(cbind(dr.r.nh.new, mdr.new, xpert.dr.r.new), na.rm =
                                             T))]
dre[sel3, rr.new.Den := as.integer(rowSums(cbind(dst.rlt.new, xpert.new), na.rm =
                                             T))]
dre[sel3, prop.rr.new := rr.new.Num / rr.new.Den]
dre[sel3 & rr.new.Den == 0, prop.rr.new := NA]
dre[is.nan(prop.rr.new), prop.rr.new := NA]
sel3 <- sel3 & !is.na(dre$prop.rr.new)
dre[sel3, test.isbinom(prop.rr.new)]


sel4 <- dre$source.ret == 'Surveillance' & dre$year.ret < 2017
table(sel4)
dre[sel4, rr.ret.Num := as.integer(rowSums(cbind(dr.r.nh.ret, mdr.ret, xpert.dr.r.ret), na.rm =
                                             T))]
dre[sel4, rr.ret.Den := as.integer(rowSums(cbind(dst.rlt.ret, xpert.ret), na.rm =
                                             T))]
dre[sel4, prop.rr.ret := rr.ret.Num / rr.ret.Den]
dre[sel4 & rr.ret.Den == 0, prop.rr.ret := NA]
dre[is.nan(prop.rr.ret), prop.rr.ret := NA]
sel4 <- sel4 & !is.na(dre$prop.rr.ret)
dre[sel4, test.isbinom(prop.rr.ret)]






#' SDs prop RR, surveillance
#' assume binomial errors
#'
out <- dre[sel1 | sel3, {
  tmp = cii(rr.new.Den, rr.new.Num)
  
  list(
    prop.rr.new.sd = tmp$se,
    prop.rr.new.lo = tmp$lower95ci,
    prop.rr.new.hi = tmp$upper95ci
  )
}, by = .(iso3)]
dre <- merge(dre, out, by = 'iso3', all.x = T)
dre[sel1 | sel3, test.isbinom(prop.rr.new.sd)]

out <- dre[sel2 | sel4, {
  tmp = cii(rr.ret.Den, rr.ret.Num)
  
  list(
    prop.rr.ret.sd = tmp$se,
    prop.rr.ret.lo = tmp$lower95ci,
    prop.rr.ret.hi = tmp$upper95ci
  )
}, by = .(iso3)]
dre <- merge(dre, out, by = 'iso3', all.x = T)
dre[sel2 | sel4, test.isbinom(prop.rr.ret.sd)]


(dre[prop.rr.new == 0, .(iso3,
                         prop.rr.new,
                         prop.rr.new.lo,
                         prop.rr.new.hi,
                         prop.rr.new.sd)])

dre[prop.rr.new.sd == 0, prop.rr.new.sd := (prop.rr.new.hi - prop.rr.new.lo) /
      3.92]
dre[prop.rr.ret.sd == 0, prop.rr.ret.sd := (prop.rr.ret.hi - prop.rr.ret.lo) /
      3.92]





#' prop MDR, surveillance
#'
dre[sel1, prop.mdr.new := mdr.new / dst.rlt.new]
dre[sel2, prop.mdr.ret := mdr.ret / dst.rlt.ret]

dre[prop.mdr.new > prop.rr.new, prop.mdr.new := prop.rr.new]
dre[prop.mdr.ret > prop.rr.ret, prop.mdr.ret := prop.rr.ret]





#' prop RR, survey >=2018
#'
sel5 <- dre$source.new != 'Surveillance' & dre$year.new >= 2018
table(sel5)
dre[sel5, prop.rr.new := rr.new.pcnt / 100]
dre[sel5, prop.rr.new.lo := rr.new.pcnt.lo / 100]
dre[sel5, prop.rr.new.hi := rr.new.pcnt.hi / 100]
dre[sel5, prop.rr.new.sd := (prop.rr.new.hi - prop.rr.new.lo) / 3.92]

sel6 <- dre$source.ret != 'Surveillance' & dre$year.ret >= 2018
table(sel6)
dre[sel6, prop.rr.ret := rr.ret.pcnt / 100]
dre[sel6, prop.rr.ret.lo := rr.ret.pcnt.lo / 100]
dre[sel6, prop.rr.ret.hi := rr.ret.pcnt.hi / 100]
dre[sel6, prop.rr.ret.sd := (prop.rr.ret.hi - prop.rr.ret.lo) / 3.92]



#' prop RR, survey <2018, phenotypic DST
#'
sel7 <- dre$source.new != 'Surveillance' & dre$year.new < 2018
table(sel7)
dre[sel7 &
      !is.na(rr.new.pcnt), prop.rr.new := rr.new.pcnt / 100, by = 'iso3']
dre[sel7 &
      !is.na(rr.new.pcnt), prop.rr.new.sd := (rr.new.pcnt.hi - rr.new.pcnt.lo) /
      3.92 / 100]
dre[sel7 &
      !is.na(rr.new.pcnt), prop.rr.new.lo := rr.new.pcnt.lo / 100]
dre[sel7 &
      !is.na(rr.new.pcnt), prop.rr.new.hi := rr.new.pcnt.hi / 100]

dre[sel7 &
      !is.na(mdr.new.pcnt), prop.mdr.new := mdr.new.pcnt / 100, by = 'iso3']
dre[sel7 &
      !is.na(mdr.new.pcnt), prop.mdr.new.sd := (mdr.new.pcnt.hi - mdr.new.pcnt.lo) /
      3.92 / 100]
dre[sel7 &
      !is.na(mdr.new.pcnt), prop.mdr.new.lo := mdr.new.pcnt.lo / 100]
dre[sel7 &
      !is.na(mdr.new.pcnt), prop.mdr.new.hi := mdr.new.pcnt.hi / 100]


dre[sel7 &
      is.na(rr.new.pcnt), prop.rr.new := rowSums(cbind(dr.r.nh.new.pcnt, mdr.new.pcnt) /
                                                   100, na.rm = T), by = 'iso3']
dre[sel7 &
      is.na(rr.new.pcnt), prop.rr.new.sd := sqrt(((
        dr.r.nh.new.pcnt.hi - dr.r.nh.new.pcnt.lo
      ) / 100 / 3.92) ^ 2 +
        ((mdr.new.pcnt.hi - mdr.new.pcnt.lo) /
           100 / 3.92) ^ 2)]

dre[sel7 &
      is.na(prop.rr.new.sd), prop.rr.new.sd := (mdr.new.pcnt.hi - mdr.new.pcnt.lo) /
      100 / 3.92]
dre[sel7 &
      !is.na(prop.rr.new) &
      is.na(prop.rr.new.sd), .(
        iso3,
        source.new,
        year.new,
        dr.r.nh.new.pcnt,
        mdr.new.pcnt,
        prop.rr.new,
        prop.rr.new.sd,
        mdr.new.pcnt.hi,
        mdr.new.pcnt.lo,
        dr.r.nh.new.pcnt.hi,
        dr.r.nh.new.pcnt.lo
      )]
dre[sel7 &
      !is.na(prop.rr.new) & is.na(prop.rr.new.sd), prop.rr.new := NA]



#' bug fix (Matteo, July 2018)
#'
sel <-
  sel7 &
  (is.na(dre$prop.rr.new.lo) |
     is.na(dre$prop.rr.new.hi)) &
  !is.na(dre$xpert.dr.r.new.pcnt.lo) &
  !is.na(dre$xpert.dr.r.new.pcnt.hi)
table(sel)
dre[sel, prop.rr.new := xpert.dr.r.new.pcnt / 100]
dre[sel, prop.rr.new.lo := xpert.dr.r.new.pcnt.lo / 100]
dre[sel, prop.rr.new.hi := xpert.dr.r.new.pcnt.hi / 100]
dre[sel, prop.rr.new.sd := (prop.rr.new.hi - prop.rr.new.lo) / 3.92]



sel8 <- dre$source.ret != 'Surveillance' & dre$year.ret < 2018
table(sel8)

dre[sel8 &
      !is.na(rr.ret.pcnt), prop.rr.ret := rr.ret.pcnt / 100, by = 'iso3']
dre[sel8 &
      !is.na(rr.ret.pcnt), prop.rr.ret.sd := (rr.ret.pcnt.hi - rr.ret.pcnt.lo) /
      3.92 / 100]
dre[sel8 &
      !is.na(rr.ret.pcnt), prop.rr.ret.lo := rr.ret.pcnt.lo / 100]
dre[sel8 &
      !is.na(rr.ret.pcnt), prop.rr.ret.hi := rr.ret.pcnt.hi / 100]

dre[sel8 &
      !is.na(mdr.ret.pcnt), prop.mdr.ret := mdr.ret.pcnt / 100, by = 'iso3']
dre[sel8 &
      !is.na(mdr.ret.pcnt), prop.mdr.ret.sd := (mdr.ret.pcnt.hi - mdr.ret.pcnt.lo) /
      3.92 / 100]
dre[sel8 &
      !is.na(mdr.ret.pcnt), prop.mdr.ret.lo := mdr.ret.pcnt.lo / 100]
dre[sel8 &
      !is.na(mdr.ret.pcnt), prop.mdr.ret.hi := mdr.ret.pcnt.hi / 100]


dre[sel8 &
      is.na(rr.ret.pcnt), prop.rr.ret := rowSums(cbind(dr.r.nh.ret.pcnt, mdr.ret.pcnt) /
                                                   100, na.rm = T), by = 'iso3']
dre[sel8 &
      is.na(rr.ret.pcnt), prop.rr.ret.sd := sqrt(((
        dr.r.nh.ret.pcnt.hi - dr.r.nh.ret.pcnt.lo
      ) / 100 / 3.92) ^ 2 +
        ((mdr.ret.pcnt.hi -
            mdr.ret.pcnt.lo) / 100 / 3.92) ^ 2)]

dre[sel8 &
      is.na(prop.rr.ret.sd), prop.rr.ret.sd := (mdr.ret.pcnt.hi - mdr.ret.pcnt.lo) /
      100 / 3.92]
dre[sel8 &
      !is.na(prop.rr.ret) &
      is.na(prop.rr.ret.sd), .(
        iso3,
        source.ret,
        year.ret,
        dr.r.nh.ret.pcnt,
        mdr.ret.pcnt,
        prop.rr.ret,
        prop.rr.ret.sd,
        mdr.ret.pcnt.hi,
        mdr.ret.pcnt.lo,
        dr.r.nh.ret.pcnt.hi,
        dr.r.nh.ret.pcnt.lo
      )]
dre[sel8 &
      !is.na(prop.rr.ret) & is.na(prop.rr.ret.sd), prop.rr.ret := NA]






#' bug fix (Matteo, July 2018)
#'
sel <-
  sel8 &
  (is.na(dre$prop.rr.ret.lo) |
     is.na(dre$prop.rr.ret.hi)) &
  !is.na(dre$xpert.dr.r.ret.pcnt.lo) &
  !is.na(dre$xpert.dr.r.ret.pcnt.hi)
table(sel)
dre[sel, prop.rr.ret := xpert.dr.r.ret.pcnt / 100]
dre[sel, prop.rr.ret.lo := xpert.dr.r.ret.pcnt.lo / 100]
dre[sel, prop.rr.ret.hi := xpert.dr.r.ret.pcnt.hi / 100]
dre[sel, prop.rr.ret.sd := (prop.rr.ret.hi - prop.rr.ret.lo) / 3.92]


dre[sel7 & !is.na(prop.rr.new), test.ispos(prop.rr.new.sd)]
dre[sel8 & !is.na(prop.rr.ret), test.ispos(prop.rr.ret.sd)]




#' prop RR, survey <2018, xpert
#'
sel9 <- dre$source.new != 'Surveillance' &
  dre$year.new %in% 2013:2017
table(sel9)
table(sel9 & is.na(dre$prop.rr.new))
sel9 <- sel9 & is.na(dre$prop.rr.new)

dre[sel9, prop.rr.new := xpert.dr.r.new.pcnt / 100]
dre[sel9, prop.rr.new.lo := xpert.dr.r.new.pcnt.lo / 100]
dre[sel9, prop.rr.new.hi := xpert.dr.r.new.pcnt.hi / 100]
dre[sel9, prop.rr.new.sd := (prop.rr.new.hi - prop.rr.new.lo) / 3.92]

sel10 <- sel8 & is.na(dre$prop.rr.ret)
table(sel10)
dre[sel10, prop.rr.ret := xpert.dr.r.ret.pcnt / 100]
dre[sel10, prop.rr.ret.lo := xpert.dr.r.ret.pcnt.lo / 100]
dre[sel10, prop.rr.ret.hi := xpert.dr.r.ret.pcnt.hi / 100]
dre[sel10, prop.rr.ret.sd := (prop.rr.ret.hi - prop.rr.ret.lo) / 3.92]

dre[sel9 & !is.na(prop.rr.new), test.ispos(prop.rr.new.sd)]
dre[sel10 & !is.na(prop.rr.ret), test.ispos(prop.rr.ret.sd)]



#' prop RR, missing bounds
#'
sel <-
  !is.na(dre$prop.rr.new) &
  (is.na(dre$prop.rr.new.lo) | is.na(dre$prop.rr.new.hi))
table(sel)
out <- vlohi(dre$prop.rr.new[sel], dre$prop.rr.new.sd[sel])
dre[sel, prop.rr.new.lo := out[1,]]
dre[sel, prop.rr.new.hi := out[2,]]

dre[!is.na(prop.rr.new), test.bounds(prop.rr.new, prop.rr.new.lo, prop.rr.new.hi)]


sel <-
  !is.na(dre$prop.rr.ret) &
  (is.na(dre$prop.rr.ret.lo) | is.na(dre$prop.rr.ret.hi))
table(sel)
out <- vlohi(dre$prop.rr.ret[sel], dre$prop.rr.ret.sd[sel])
dre[sel, prop.rr.ret.lo := out[1,]]
dre[sel, prop.rr.ret.hi := out[2,]]

dre[!is.na(prop.rr.ret), test.bounds(prop.rr.ret, prop.rr.ret.lo, prop.rr.ret.hi)]




#' prop MDR, surveys
#'
seln <- dre$source.new != 'Surveillance' & dre$surv.quality.new == 'A'
dre[seln, prop.mdr.new := mdr.new.pcnt / 100]

selr <-
  dre$source.ret != 'Surveillance' &
  dre$surv.quality.ret %in% c('A', 'R', 'S')
dre[selr, prop.mdr.ret := mdr.ret.pcnt / 100]






#' # Imputations
#'
#' pool props by region
#'
suppressWarnings(fit.new <-
                   rma(
                     yi = prop.rr.new,
                     sei = prop.rr.new.sd,
                     data = dre,
                     mods =  ~ g.mdr - 1
                   ))
suppressWarnings(fit.ret <-
                   rma(
                     yi = prop.rr.ret,
                     sei = prop.rr.ret.sd,
                     data = dre,
                     mods =  ~ g.mdr - 1
                   ))

imp.new <- fit.new$b[, 1]
imp.ret <- fit.ret$b[, 1]

imp.new.sd <- fit.new$se
imp.ret.sd <- fit.ret$se
imp.new.ui <- vlohi(unlist(imp.new), imp.new.sd)
imp.ret.ui <- vlohi(unlist(imp.ret), imp.ret.sd)

reg <- names(table(dre$g.mdr))

#' check the order of regions in reg is the same as in imp*
#'
all.equal(reg, gsub('g.mdr', '', dimnames(fit.new$beta)[[1]]))
all.equal(reg, gsub('g.mdr', '', dimnames(fit.ret$beta)[[1]]))

table(!is.na(dre$prop.rr.new))
table(!is.na(dre$prop.rr.ret))

tmp <- copy(dre)

#' not elegant, should set up a proper merge
#'
for (i in 1:length(reg)) {
  sel <- dre$g.mdr == reg[i] & is.na(dre$prop.rr.new)
  dre[sel, prop.rr.new := imp.new[i]]
  dre[sel, prop.rr.new.sd := imp.new.sd[i]]
  dre[sel, prop.rr.new.lo := imp.new.ui[1, i]]
  dre[sel, prop.rr.new.hi := imp.new.ui[2, i]]
  
  sel <- dre$g.mdr == reg[i] & is.na(dre$prop.rr.ret)
  dre[sel, prop.rr.ret := imp.ret[i]]
  dre[sel, prop.rr.ret.sd := imp.ret.sd[i]]
  dre[sel, prop.rr.ret.lo := imp.ret.ui[1, i]]
  dre[sel, prop.rr.ret.hi := imp.ret.ui[2, i]]
}


dre[, lapply(.SD, test.ispos), .SDcols = c('prop.rr.new',
                                           'prop.rr.new.sd',
                                           'prop.rr.ret',
                                           'prop.rr.ret.sd')]



#' # Ratio of MDR/RR
#'
table(is.na(dre$prop.mdr.new))
table(is.na(dre$prop.mdr.ret))

tmp <- copy(dre)
dre[, mdr.rr.new := pmin(prop.mdr.new / prop.rr.new, 1)]
dre[, mdr.rr.ret := pmin(prop.mdr.ret / prop.rr.ret, 1)]




#' # MDR checks
#'
#' mdr.rr.new and mdr.rr.ret must not exceed 1
#'
dre[!is.na(mdr.rr.new), test.isbinom(mdr.rr.new)]


#' impute missing
#'
ratio.mdr.new <- dre[, median(mdr.rr.new, na.rm = T), by = g.mdr]
ratio.mdr.ret <- dre[, median(mdr.rr.ret, na.rm = T), by = g.mdr]

table(is.na(dre$mdr.rr.new))
table(is.na(dre$mdr.rr.ret))


#' impute where missing using the regional median
#'
for (i in reg) {
  sel <- is.na(dre$mdr.rr.new) & dre$g.mdr == i
  sel2 <- ratio.mdr.new$g.mdr == i
  dre[sel, mdr.rr.new := ratio.mdr.new$V1[sel2]]
  
  sel <- is.na(dre$mdr.rr.ret) & dre$g.mdr == i
  sel2 <- ratio.mdr.ret$g.mdr == i
  dre[sel, mdr.rr.ret := ratio.mdr.ret$V1[sel2]]
}

dre[, test.isbinom(mdr.rr.new)]
dre[, test.isbinom(mdr.rr.ret)]





#' # RR incidence
#'
#' $I_\text{rr} = I \left[ (1-f) p_n ((1-r) + r \rho) + fp_r \right]$
#'
#' $I$ = TB incidence
#'
#' $f$ = proportion retreatment for failure or default
#'
#' $p_n$ = prop.rr.new
#'
#' $r$ = proportion relapses
#'
#' $\rho$ = risk ratio of rr in relapse vs new
#'
#' $p_r$ = prop.rr.ret
#'
incdr <-
  function(inc,
           inc.sd,
           f,
           f.sd,
           pn,
           pn.sd,
           r,
           r.sd,
           rho,
           rho.sd,
           pr,
           pr.sd) {
    require(propagate)
    
    DT <- cbind(
      inc = c(inc, inc.sd),
      f = c(f, f.sd),
      pn = c(pn, pn.sd),
      r = c(r, r.sd),
      rho = c(rho, rho.sd),
      pr = c(pr, pr.sd)
    )
    EXPR <- expression(inc * ((1 - f) * pn * ((1 - r) + r * rho) + f * pr))
    
    out <-
      propagate(
        expr = EXPR,
        data = DT,
        type = 'stat',
        do.sim = F,
        second.order = T
      )
  }


#' prep data, f, r
#'
rr <-
  dre[, .(
    iso3,
    g.whoregion,
    g.mdr,
    year.new,
    source.new,
    all.areas.covered.new,
    surv.quality.new,
    year.ret,
    source.ret,
    all.areas.covered.ret,
    surv.quality.ret,
    prop.rr.new,
    prop.rr.new.sd,
    prop.rr.new.lo,
    prop.rr.new.hi,
    prop.rr.ret,
    prop.rr.ret.sd,
    prop.rr.ret.lo,
    prop.rr.ret.hi,
    mdr.rr.new,
    mdr.rr.ret
  )]


rrtb <- tb[year > yr - 5, .(iso3, year, c.newinc, c.ret, ret.nrel)]
setkey(rrtb, iso3, year)

rrtb[, ret.rel := c.ret - ret.nrel]
rrtb[!is.na(ret.rel), test.ispos(ret.rel)]

rrtb[, f := ret.nrel / c.newinc]
rrtb[c.newinc == 0, f := 1]
rrtb[c.newinc < c.ret, f := 1]
rrtb[!is.na(f), test.isbinom(f)]

rrtb[f < 1, r := pmin(ret.rel / c.newinc / (1 - f), 1)]
rrtb[f == 1, r := 0]
rrtb[!is.na(r), test.isbinom(r)]

out <- rrtb[, .(
  r = mean(r, na.rm = T),
  r.sd = sd(r, na.rm = T),
  f = mean(f, na.rm = T),
  f.sd = sd(f, na.rm = T)
), by = iso3]

out[!is.na(r), test.ispos(r)]
out[!is.na(r.sd), test.ispos(r.sd)]
out[!is.na(f), test.ispos(f)]
out[!is.na(f.sd), test.ispos(f.sd)]
out[is.na(r.sd) & !is.na(r), r.sd := r * .2]
out[is.na(f.sd) & !is.na(f), f.sd := f * .2]

out[is.nan(r), r := NA]
out[is.nan(r.sd), r.sd := NA]
out[is.nan(f), f := NA]
out[is.nan(f.sd), f.sd := NA]
out[is.infinite(f), f := NA]
out[is.infinite(f.sd), f.sd := NA]


#' imputation - pool r and f by mdr region
out <- merge(out, est[year == yr, .(iso3, g.mdr)], by = 'iso3')
fit.r <- rma(
  yi = r,
  sei = r.sd,
  data = out[r.sd > 0],
  mods =  ~ g.mdr - 1
)
fit.f <- rma(
  yi = f,
  sei = f.sd,
  data = out[f.sd > 0],
  mods =  ~ g.mdr - 1
)

imp.r <- fit.r$b[, 1]
imp.f <- fit.f$b[, 1]
imp.r.sd <- fit.r$se
imp.f.sd <- fit.f$se

#' check the order of regions in reg is the same as in imp*
all.equal(reg, gsub('g.mdr', '', dimnames(fit.r$beta)[[1]]))
all.equal(reg, gsub('g.mdr', '', dimnames(fit.f$beta)[[1]]))

out[, g.mdr := NULL]
rr <- merge(rr, out, by = 'iso3', all.x = T)

table(!is.na(rr$r))
table(!is.na(rr$f))

for (i in 1:length(reg)) {
  sel <- rr$g.mdr == reg[i] & (is.na(rr$r) | rr$r == 0)
  rr[sel, r := imp.r[i]]
  rr[sel, r.sd := imp.r.sd[i]]
  
  sel <- rr$g.mdr == reg[i] & (is.na(rr$f) | rr$f == 0)
  rr[sel, f := imp.f[i]]
  rr[sel, f.sd := imp.f.sd[i]]
}



#' $\rho$ and its SD
#' (see rel.R)
#' 
# rho <- 3.2798
# rho.sd <- .3247
rho <- 5.5 # reverting to last year's after USAID discussions
rho.sd <- (6.8 - 4.4) / 3.92
rho.cv <- rho.sd / rho

rel2 <- rel[level=='national',.(rho.rr=last(ratio.rel.new), 
                               rho.sd.rr=last(ratio.rel.new.sd)), by=iso3]
rel3 <- rel2[rho.rr<=rho & rho.rr>1]
rr2 <- merge(rr, rel3, by='iso3', all.x=TRUE)
dim(rr); dim(rr2)
rr <- copy(rr2)

rr[!is.na(rho.rr), rho.rr := ifelse(prop.rr.ret / prop.rr.new > rho.rr, rho.rr, prop.rr.ret /
                        prop.rr.new)]

rr[is.na(rho.rr), rho.rr := ifelse(prop.rr.ret / prop.rr.new > rho, rho, prop.rr.ret / prop.rr.new)]
rr[is.na(rho.sd.rr), rho.sd.rr := rho.rr * rho.cv]
rr[is.na(rho.rr), rho.rr := rho]
rr[is.na(rho.sd.rr), rho.sd.rr := rho.sd]


rr[, rho.mdr := ifelse(
  prop.rr.ret * mdr.rr.ret / (prop.rr.new * mdr.rr.new) > rho.rr,
  rho.rr,
  prop.rr.ret * mdr.rr.ret / (prop.rr.new * mdr.rr.new)
)]
rr[, rho.sd.mdr := rho.cv * rho.mdr]
rr[is.na(rho.mdr), rho.mdr := rho]
rr[is.na(rho.sd.mdr), rho.sd.mdr := rho.sd]


#' inc.rr
#'
rr <- merge(rr, est[year == yr, .(iso3, inc, inc.sd, pop)], all.x = T)
out <-
  rr[, {
    tmp = incdr(
      inc,
      inc.sd,
      f,
      f.sd,
      pn = prop.rr.new,
      pn.sd = prop.rr.new.sd,
      r,
      r.sd,
      rho = rho.rr,
      rho.sd = rho.sd.rr,
      pr = prop.rr.ret,
      pr.sd = prop.rr.ret.sd
    )$prop
    
    list(inc.rr = tmp[2],
         inc.rr.sd = tmp[4])
  },
  by = .(iso3)]

rr <- merge(rr, out, by = 'iso3', all.x = T)




#' # MDR incidence
#'
#' $I_\text{mdr} = I \left[ (1-f) q_n ((1-r) + r \rho) + fq_r \right]$
#'
#' $q_n$ = prop.rr.new * mdr.rr.new
#'
#' $q_r$ = prop.rr.ret * mdr.rr.ret
#'

#' inc.mdr
#'
out2 <-
  rr[, {
    tmp = incdr(
      inc,
      inc.sd,
      f,
      f.sd,
      pn = prop.rr.new * mdr.rr.new,
      pn.sd = prop.rr.new.sd * mdr.rr.new,
      r,
      r.sd,
      rho = rho.mdr,
      rho.sd = rho.sd.mdr,
      pr = prop.rr.ret * mdr.rr.ret,
      pr.sd = prop.rr.ret.sd * mdr.rr.ret
    )$prop
    
    list(inc.mdr = tmp[2],
         inc.mdr.sd = tmp[4])
  },
  by = .(iso3)]

rr <- merge(rr, out2, by = 'iso3', all.x = T)

sel <- rr$prop.rr.new > 0 & rr$inc > 0

rr[sel, lapply(.SD, test.ispos), .SDcols = c('inc.rr', 'inc.rr.sd', 'inc.mdr', 'inc.mdr.sd')]

rr[inc.rr < inc.mdr] # should be empty


#' bounds
#'
m <- 1e5
sel <- rr$inc.rr.sd > rr$inc.rr
table(sel)

(rr[sel, .(iso3, inc.rr, inc.rr.sd, inc.rr / m * (1 - inc.rr /
                                                              m))])

rr[sel, inc.rr.sd := inc.rr * .25]

sel <- rr$inc.rr > 0 & !is.na(rr$inc.rr)
table(sel)
out <- vlohi(rr$inc.rr[sel] / m, rr$inc.rr.sd[sel] / m)
rr[sel, inc.rr.lo := out[1,] * m]
rr[sel, inc.rr.hi := out[2,] * m]

rr[inc.rr == 0, inc.rr.lo := 0]
rr[inc.rr == 0, inc.rr.hi := 1.96 * inc.rr.sd]

rr[!is.na(inc.rr), test.bounds(inc.rr, inc.rr.lo, inc.rr.hi)]


sel <- rr$inc.mdr.sd > rr$inc.mdr
table(sel)

(rr[sel, .(iso3, inc.mdr, inc.mdr.sd, inc.mdr / m * (1 - inc.mdr /
                                                                 m))])

rr[sel, inc.mdr.sd := inc.mdr * .25]

sel <- rr$inc.mdr > 0 & !is.na(rr$inc.mdr)
table(sel)
out <- vlohi(rr$inc.mdr[sel] / m, rr$inc.mdr.sd[sel] / m)
rr[sel, inc.mdr.lo := out[1,] * m]
rr[sel, inc.mdr.hi := out[2,] * m]

rr[inc.mdr == 0, inc.mdr.lo := 0]
rr[inc.mdr == 0, inc.mdr.hi := 1.96 * inc.mdr.sd]

rr[!is.na(inc.mdr), test.bounds(inc.mdr, inc.mdr.lo, inc.mdr.hi)]




#' check key params visually
#'

#+ fig.width=12, fig.height=8
whomap(rr[, .(iso3, var = cut(f, b = c(0, .02, 0.05, 0.1, 0.2, 0.3, 1)))]) +
  scale_fill_brewer('failure/default\nrate')
ggsave(here('output/checks/rr_failure_map.pdf'), width=10, height=8)

whomap(rr[, .(iso3, var = cut(r, b = c(0, .02, 0.05, 0.1, 0.2, 0.3, 1)))]) +
  scale_fill_brewer('relapse rate')
ggsave(here('output/checks/rr_relapse_map.pdf'), width=10, height=8)



#' check output visually
#'
#'
#+ fig.width=12, fig.height=8
whomap(rr[, .(iso3, var = cut(mdr.rr.new, b = c(0, .5, 0.6, 0.7, 0.9, 1)))])
whomap(rr[, .(iso3, var = cut(mdr.rr.ret, b = c(0, .5, 0.6, 0.7, 0.9, 1)))])
whomap(rr[, .(iso3, var = cut(inc.rr, b = c(0, 0.2, 0.5, 1, 2, 5, 10, Inf)))]) +
  scale_fill_brewer('RR incidence\nper 100 000/year')
ggsave(here('output/checks/rr_incidence_map.pdf'), width=10, height=8)




rr[, inc.rr.num := inc.rr * pop / 1e5]
rr[, inc.rr.lo.num := inc.rr.lo * pop / 1e5]
rr[, inc.rr.hi.num := inc.rr.hi * pop / 1e5]

rr[, inc.mdr.num := inc.mdr * pop / 1e5]
rr[, inc.mdr.lo.num := inc.mdr.lo * pop / 1e5]
rr[, inc.mdr.hi.num := inc.mdr.hi * pop / 1e5]


rr <- rr[iso3 %ni% c('ANT', 'SCG')]
rr <- merge(rr, est[year == yr, .(iso3, g.hbmdr)], by = 'iso3', all.x =
              T)

#' save
#'
save(rr, file = here('Rdata/rr.Rdata'))
fwrite(rr, file = here(paste0('csv/rr_', Sys.Date(), '.csv')))

load(here('../gtb2019/Rdata/rr.Rdata'))
orr <- copy(rr)
save(orr, file = here('Rdata/orr.Rdata'))

load(here('Rdata/rr.Rdata'))


#' # RR aggregates and Table 3.9 RR incidence
M <- 1e5
m <- 1e3

#' incidence
#'
global.rr <-
  data.table(rr[, addXY(inc.rr / M, r.sd = inc.rr.sd / M, weights = pop)])
global.rr[, region := 'Global']

regional.rr <-
  data.table(rr[, addXY(inc.rr / M, r.sd = inc.rr.sd / M, weights = pop), by =
                  g.whoregion])
setnames(regional.rr, 'g.whoregion', 'region')

hb.rr <-
  data.table(rr[g.hbmdr == TRUE, addXY(inc.rr / M, r.sd = inc.rr.sd / M, weights =
                                         pop)])
hb.rr[, region := 'MDR/RR HBCs']


global.mdr <-
  data.table(rr[, addXY(inc.mdr / M, r.sd = inc.mdr.sd / M, weights = pop)])
global.mdr[, region := 'Global']

regional.mdr <-
  data.table(rr[, addXY(inc.mdr / M, r.sd = inc.mdr.sd / M, weights = pop), by =
                  g.whoregion])
setnames(regional.mdr, 'g.whoregion', 'region')

hb.mdr <-
  data.table(rr[g.hbmdr == TRUE, addXY(inc.mdr / M, r.sd = inc.mdr.sd / M, weights =
                                         pop)])
hb.mdr[, region := 'MDR/RR HBCs']

A <- rr[g.hbmdr == TRUE, .(
  region = iso3,
  inc.rr.num = inc.rr.num / m,
  inc.rr.num.lo = inc.rr.lo.num / m,
  inc.rr.num.hi = inc.rr.hi.num / m,
  inc.rr,
  inc.rr.lo,
  inc.rr.hi
)]
B <- hb.rr[, .(
  region,
  inc.rr.num = r.num / m,
  inc.rr.num.lo = r.lo.num / m,
  inc.rr.num.hi = r.hi.num / m,
  inc.rr = r * M,
  inc.rr.lo = r.lo * M,
  inc.rr.hi = r.hi * M
)]

C <- regional.rr[, .(
  region,
  inc.rr.num = r.num / m,
  inc.rr.num.lo = r.lo.num / m,
  inc.rr.num.hi = r.hi.num / m,
  inc.rr = r * M,
  inc.rr.lo = r.lo * M,
  inc.rr.hi = r.hi * M
)]
D <- global.rr[, .(
  region,
  inc.rr.num = r.num / m,
  inc.rr.num.lo = r.lo.num / m,
  inc.rr.num.hi = r.hi.num / m,
  inc.rr = r * M,
  inc.rr.lo = r.lo * M,
  inc.rr.hi = r.hi * M
)]

tab1 <- rbind(A, B, C, D)

mdrA <- rr[g.hbmdr == T, inc.mdr / inc.rr]
mdrB <- rr[g.hbmdr == T, sums(inc.mdr.num) / sums(inc.rr.num)]
mdrC <- rr[, sums(inc.mdr.num) / sums(inc.rr.num), by = g.whoregion]
mdrD <- rr[, sums(inc.mdr.num) / sums(inc.rr.num)]

tab1[, pctMDR := 100 * c(mdrA, mdrB, mdrC$V1, mdrD)]

(tab1)


#' prop RR
#'
rr[, new := inc * pop * pmax(0, 1 - f - r) / 1e5]
rr[, ret := inc * pop * pmin(f + r, 1) / 1e5]

global.prop <- rr[, .(
  prop.rr.new = weighted.mean(prop.rr.new, w = new),
  prop.rr.new.sd = sqrt(sum(prop.rr.new.sd ^ 2 * new /
                              sum(new))),
  prop.rr.ret = weighted.mean(prop.rr.ret, w = ret),
  prop.rr.ret.sd = sqrt(sum(prop.rr.ret.sd ^ 2 * new /
                              sum(ret)))
)]
global.prop[, region := 'Global']

hb.prop <-
  rr[g.hbmdr == TRUE, .(
    prop.rr.new = weighted.mean(prop.rr.new, w = new),
    prop.rr.new.sd = sqrt(sum(prop.rr.new.sd ^
                                2 * new / sum(new))),
    prop.rr.ret = weighted.mean(prop.rr.ret, w =
                                  ret),
    prop.rr.ret.sd = sqrt(sum(prop.rr.ret.sd ^
                                2 * new / sum(ret)))
  )]
hb.prop[, region := 'MDR/RR HBCs']

regional.prop <-
  rr[, .(
    prop.rr.new = weighted.mean(prop.rr.new, w = new),
    prop.rr.new.sd = sqrt(sum(prop.rr.new.sd ^ 2 * new /
                                sum(new))),
    prop.rr.ret = weighted.mean(prop.rr.ret, w = ret),
    prop.rr.ret.sd = sqrt(sum(prop.rr.ret.sd ^ 2 * new /
                                sum(ret)))
  ),
  by = g.whoregion]
setnames(regional.prop, 'g.whoregion', 'region')


tab4 <- rbind(hb.prop, regional.prop, global.prop)


#' bounds
#'
out1 <- vlohi(tab4$prop.rr.new, tab4$prop.rr.new.sd)
out2 <- vlohi(tab4$prop.rr.ret, tab4$prop.rr.ret.sd)

tab4[, prop.rr.new.lo := out1[1,]]
tab4[, prop.rr.new.hi := out1[2,]]

tab4[, prop.rr.ret.lo := out2[1,]]
tab4[, prop.rr.ret.hi := out2[2,]]

tab4[, prop.rr.new.sd := NULL]
tab4[, prop.rr.ret.sd := NULL]
cty.prop <-
  rr[g.hbmdr == TRUE, .(
    region = iso3,
    prop.rr.new,
    prop.rr.new.lo,
    prop.rr.new.hi,
    prop.rr.ret,
    prop.rr.ret.lo,
    prop.rr.ret.hi
  )]

tab5 <- rbind(cty.prop, tab4)

tab5[, test.bounds(prop.rr.new, prop.rr.new.lo, prop.rr.new.hi)]
tab5[, test.bounds(prop.rr.ret, prop.rr.ret.lo, prop.rr.ret.hi)]
tab6 <- merge(tab5, tab1, by = 'region', sort = FALSE)
tab6[, 2:7 := lapply(.SD, function(x)
  ftb(x * 100)), .SDcols = 2:7]
tab6[, 8:14 := lapply(.SD, function(x)
  ftb(x)), .SDcols = 8:14]

fwrite(tab6, file = here('report/tabs/tab4_10_RRincidence.csv'))

tab7 <- merge(tab5, tab1, by = 'region', sort = FALSE)
rrg <- tab7[31:38]
rrg[, 8:10 := lapply(.SD, function(x)
  x * m), .SDcols = 8:10]



#' # Global results
#'
#' global RR incidence = `r rr[, as.integer(sums(inc.rr.num))]`
#'
#' global MDR incidence = `r rr[, as.integer(sums(inc.mdr.num))]`
#'

#' save
save(rrg, file = here('Rdata/rrg.Rdata'))
fwrite(rrg, file = here(paste0('csv/rrg_', Sys.Date(), '.csv')))


rrg[, 8:10 := lapply(.SD, function(x)
  as.integer(x)), .SDcols = 8:10]
(rrg)






#' # Comparisons with last year, grouped by WHO region
#'
wr <- unique(rr$g.whoregion)
# drh <- merge(drh, cty[, .(iso3, g.whoregion)], by = 'iso3')

#' prop rr in new
#'
#' blue: new estimates
#'
#' orange: last year's estimates
#'
for (i in wr) {
  p <-
    qplot(prop.rr.new,
          iso3,
          data = rr[g.whoregion == i],
          colour = I('blue')) +
    geom_segment(
      aes(
        x = prop.rr.new.lo,
        xend = prop.rr.new.hi,
        y = iso3,
        yend = iso3
      ),
      colour = I('blue'),
      size = I(1)
    ) +
    geom_point(aes(e.rr.prop.new, iso3),
               data = drh[g.whoregion == i &
                            year == yr - 1],
               colour = I('orange')) +
    geom_segment(
      aes(
        x = e.rr.prop.new.lo,
        xend = e.rr.prop.new.hi,
        y = iso3,
        yend = iso3
      ),
      data = drh[g.whoregion == i &
                   year == yr - 1],
      colour = I('orange'),
      alpha = I(.4),
      size = I(3)
    ) +
    xlab('Proportion RR in new cases') + ylab('')
  suppressWarnings(print(p))
  suppressWarnings(ggsave(
    filename = here(paste0('output/checks/rr_new', i, '_compare.pdf')),
    height = 8,
    width = 6
  ))
}

#' prop rr in ret
#'
for (i in wr) {
  p <-
    qplot(prop.rr.ret,
          iso3,
          data = rr[g.whoregion == i],
          colour = I('blue')) +
    geom_segment(
      aes(
        x = prop.rr.ret.lo,
        xend = prop.rr.ret.hi,
        y = iso3,
        yend = iso3
      ),
      colour = I('blue'),
      size = I(1)
    ) +
    geom_point(aes(e.rr.prop.ret, iso3),
               data = drh[g.whoregion == i &
                            year == yr - 1],
               colour = I('orange')) +
    geom_segment(
      aes(
        x = e.rr.prop.ret.lo,
        xend = e.rr.prop.ret.hi,
        y = iso3,
        yend = iso3
      ),
      data = drh[g.whoregion == i &
                   year == yr - 1],
      colour = I('orange'),
      alpha = I(.4),
      size = I(3)
    ) +
    xlab('Proportion RR in retreatment cases') + ylab('')
  suppressWarnings(print(p))
  suppressWarnings(ggsave(
    filename = here(paste0('output/checks/rr_ret', i, '_compare.pdf')),
    height = 8,
    width = 6
  ))
}


#' rr inc
#'
for (i in wr) {
  p <-
    qplot(inc.rr, iso3, data = rr[g.whoregion == i], colour = I('blue')) +
    geom_segment(
      aes(
        x = inc.rr.lo,
        xend = inc.rr.hi,
        y = iso3,
        yend = iso3
      ),
      colour = I('blue'),
      size = I(1)
    ) +
    geom_point(aes(e.inc.rr.100k, iso3),
               data = drh[g.whoregion == i &
                            year == yr - 1],
               colour = I('orange')) +
    geom_segment(
      aes(
        x = e.inc.rr.100k.lo,
        xend = e.inc.rr.100k.hi,
        y = iso3,
        yend = iso3
      ),
      data = drh[g.whoregion == i &
                   year == yr - 1],
      colour = I('orange'),
      alpha = I(.4),
      size = I(3)
    ) +
    xlab('RR incidence') + ylab('')
  suppressWarnings(print(p))
  suppressWarnings(ggsave(
    filename = here(paste0('output/checks/incrr', i, '_compare.pdf')),
    height = 8,
    width = 6
  ))
}



# Primary versus acquired DR
#
# this does not work very well, SDs too large
# prop primary RR ------------------------------------
# primarydr <- function(inc, inc.sd, pn, pn.sd, inc.rr, inc.rr.sd){
#   require(propagate)
#
#   DT <- cbind(inc=c(inc, inc.sd),
#               pn=c(pn, pn.sd),
#               inc.rr=c(inc.rr, inc.rr.sd))
#   EXPR <- expression(inc * pn / inc.rr)
#
#   out <- propagate(expr=EXPR, data=DT, type='stat', do.sim=F, second.order=T)
# }
#
# out <- rr[inc>60 & inc.rr>2,{tmp = primarydr(inc=inc, inc.sd=inc.sd, pn=prop.rr.new,
#                                              pn.sd=0, inc.rr=inc.rr,
#                                              inc.rr.sd=inc.rr.sd)$prop;
# list(primary.rr = tmp[2], primary.rr.sd = tmp[4])}, by=.(iso3)]
#
# # rr <- merge(rr, out2, by='iso3', all.x=T)
#
#
#
#
# rrn <- function(iso3='CHN', year=2017, conf=FALSE){
#   # RR cases among pulmonary cases
#   # conf = restricted to bacteriologically confirmed
#
#   seltb <- tb$iso3==iso3 & tb$year==year
#   selrr <- rr$iso3==iso3
#
#   rr.new <- rr$prop.rr.new[selrr]
#   rr.new.sd <- rr$prop.rr.new.sd[selrr]
#   rr.ret <- rr$prop.rr.ret[selrr]
#   rr.ret.sd <- rr$prop.rr.ret.sd[selrr]
#
#   if (conf){
#     cp.new <- tb$new.labconf[seltb]
#     cp.ret <- tb$ret.rel.labconf[seltb]
#   } else {
#     cp.new <- sums(tb$new.labconf[seltb], tb$new.clindx[seltb])
#     cp.ret <- sums(tb$ret.rel.labconf[seltb], tb$rel.ret.clindx[seltb])
#   }
#
#   tot <- sums(cp.new, cp.ret)
#   best <- rr.new*cp.new + rr.ret * cp.ret
#   best.sd <- sqrt((rr.new.sd^2 * cp.new^2 ) + (rr.ret.sd^2 * cp.ret^2 ))
#   lo <- best - 1.96 * best.sd
#   hi <- best + 1.96 * best.sd
#   return(c(tot, best, lo, hi))
# }
# (rrn())
# (rrn(conf=T))
#
