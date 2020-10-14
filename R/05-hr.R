#' ---
#' title: HR estimates
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
#' (Last updated: `r Sys.Date()`)
#'
#' RR incidence estimates, Global TB Report 2019
#'
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(propagate))
suppressPackageStartupMessages(library(metafor))
suppressPackageStartupMessages(library(whomap))
library(here)


load(here('Rdata/est.Rdata'))
load(here('Rdata/tb.Rdata'))
load(here('Rdata/cty.Rdata'))
load(here('Rdata/rr.Rdata'))
load(here('Rdata/rrg.Rdata'))


#' most recent HR data
#'
load(here('Rdata/hr.Rdata'))

vlohi <- Vectorize(lohi, c('ev', 'sd'))
yr <- 2019


hr <-
  merge(hr,
        est[year == 2019, .(iso3, g.mdr)],
        by = 'iso3',
        all.x = T,
        all.y = F)
hr <- hr[!is.na(g.mdr)]


#' Anna's formulas
#'
#' Any INH resistance in new (regardless of rifampicin status)
#' Data collected in years ≥2018:  dst_rlt_hr_new / dst_rlt_new
#' Data collected in years <2018: (dr_h_nr_new + mdr_new) / dst_rlt_new
#'
#' INH-resistant, RIF-susceptible TB in new
#' Data collected in years ≥2018: (dst_rlt_hr_new - mdr_new) / dst_rlt_new
#' Data collected in years <2018: dr_h_nr_new / dst_rlt_new
#'
#' Any INH resistance in retreatment (regardless of rifampicin status)
#' Data collected in years ≥2018:  dst_rlt_hr_ret / dst_rlt_ret
#' Data collected in years <2018: (dr_h_nr_ret + mdr_ret) / dst_rlt_ret
#'
#' INH-resistant, RIF-susceptible TB in retreatment
#' Data collected in years ≥2018: (dst_rlt_hr_ret - mdr_ret) / dst_rlt_ret
#' Data collected in years <2018: dr_h_nr_ret / dst_rlt_ret


#' Any INH resistance
#'
hr[year.new >= 2018, prop.hr.new := dst.rlt.hr.new / dst.rlt.new]
hr[year.new < 2018, prop.hr.new := (dr.h.nr.new + mdr.new) / dst.rlt.new]
hr[!is.na(prop.hr.new), test.isbinom(prop.hr.new)]
hr[dst.rlt.ret > 0, prop.hr.new.sd := sqrt(prop.hr.new *  (1 - prop.hr.new)  / dst.rlt.new)]

hr[year.ret >= 2018, prop.hr.ret := dst.rlt.hr.ret / dst.rlt.ret]
hr[year.ret < 2018, prop.hr.ret := (dr.h.nr.ret + mdr.ret) / dst.rlt.ret]

# exclude BLZ
hr[prop.hr.ret > 1, prop.hr.ret := NA]
hr[!is.na(prop.hr.ret), test.isbinom(prop.hr.ret)]
hr[dst.rlt.ret > 0, prop.hr.ret.sd := sqrt(prop.hr.ret *  (1 - prop.hr.ret)  / dst.rlt.ret)]


# HR RS
#
# hr[year.new >= 2018, prop.hr.rs.new := (dst.rlt.hr.new - mdr.new) / dst.rlt.new]
# hr[year.new < 2018, prop.hr.rs.new := dr.h.nr.new / dst.rlt.new]
# hr[!is.na(prop.hr.rs.new), test.isbinom(prop.hr.rs.new)]
# hr[dst.rlt.ret > 0, prop.hr.rs.new.sd := sqrt(prop.hr.rs.new *  (1 - prop.hr.rs.new)  / dst.rlt.new)]
#
# hr[year.ret >= 2018, prop.hr.rs.ret := (dst.rlt.hr.ret - mdr.ret) / dst.rlt.ret]
# hr[year.ret < 2018, prop.hr.rs.ret := dr.h.nr.ret / dst.rlt.ret]
# hr[!is.na(prop.hr.rs.ret), test.isbinom(prop.hr.rs.ret)]
# hr[dst.rlt.ret > 0, prop.hr.rs.ret.sd := sqrt(prop.hr.rs.ret *  (1 - prop.hr.rs.ret)  / dst.rlt.ret)]


#' missing values
#'
(hr[, sum(is.na(prop.hr.new))])
# (hr[, sum(is.na(prop.hr.rs.new))])
(hr[, sum(is.na(prop.hr.ret))])
# (hr[, sum(is.na(prop.hr.rs.ret))])

(hr[, sum(is.na(prop.hr.new.sd))])
# (hr[, sum(is.na(prop.hr.rs.new.sd))])
(hr[, sum(is.na(prop.hr.ret.sd))])
# (hr[, sum(is.na(prop.hr.rs.ret.sd))])

hr[is.na(prop.hr.new.sd), prop.hr.new := NA]
hr[is.na(prop.hr.ret.sd), prop.hr.ret := NA]
# hr[is.na(prop.hr.rs.new.sd), prop.hr.rs.new := NA]
# hr[is.na(prop.hr.rs.ret.sd), prop.hr.rs.ret := NA]


#' imputations, regional pooling
#'
suppressWarnings(fit.new <-
                   rma(
                     yi = prop.hr.new,
                     sei = prop.hr.new.sd,
                     data = hr,
                     mods = ~ g.mdr - 1
                   ))
suppressWarnings(fit.new.global <-
                   rma(
                     yi = prop.hr.new,
                     sei = prop.hr.new.sd,
                     data = hr,
                   ))
suppressWarnings(fit.ret <-
                   rma(
                     yi = prop.hr.ret,
                     sei = prop.hr.ret.sd,
                     data = hr,
                     mods = ~ g.mdr - 1
                   ))
imp.new <- unlist(fit.new$b[, 1])
imp.new[1] <- unlist(fit.new.global$b[, 1]) # too little data in AFR
imp.ret <- unlist(fit.ret$b[, 1])

imp.new.sd <- fit.new$se
imp.new.sd[1] <- fit.new.global$se
imp.ret.sd <- fit.ret$se

reg <- names(table(hr$g.mdr))

#' check the order of regions in reg is the same as in imp*
#'
all.equal(reg, gsub('g.mdr', '', dimnames(fit.new$beta)[[1]]))
all.equal(reg, gsub('g.mdr', '', dimnames(fit.ret$beta)[[1]]))

table(!is.na(hr$prop.hr.new))
table(!is.na(hr$prop.hr.ret))

tmp <- copy(hr)

for (i in 1:length(reg)) {
  sel <- hr$g.mdr == reg[i] & is.na(hr$prop.hr.new)
  hr[sel, prop.hr.new := imp.new[i]]
  hr[sel, prop.hr.new.sd := imp.new.sd[i]]
  
  sel <- hr$g.mdr == reg[i] & is.na(hr$prop.hr.ret)
  hr[sel, prop.hr.ret := imp.ret[i]]
  hr[sel, prop.hr.ret.sd := imp.ret.sd[i]]
}

hr[, lapply(.SD, test.ispos), .SDcols = c('prop.hr.new',
                                          'prop.hr.new.sd',
                                          'prop.hr.ret',
                                          'prop.hr.ret.sd')]


#' HR incidence
#'
#' $I_\text{rr} = I \left[ (1-f) p_n ((1-r) + r \rho) + fp_r \right]$
#'
#' $I$ = TB incidence
#'
#' $f$ = proportion retreatment for failure or default
#'
#' $p_n$ = prop.hr.new
#'
#' $r$ = proportion relapses
#'
#' $\rho$ = risk ratio of hr in relapse vs new
#'
#' $p_r$ = prop.hr.ret
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
    EXPR <-
      expression(inc * ((1 - f) * pn * ((1 - r) + r * rho) + f * pr))
    
    out <-
      propagate(
        expr = EXPR,
        data = DT,
        type = 'stat',
        do.sim = F,
        second.order = T
      )
  }




#' inc.hr
#'
hr <- merge(hr, est[year == yr, .(iso3, inc, inc.sd, pop)], all.x = T)
hr <- merge(hr, rr[, .(iso3, f, f.sd, r, r.sd, inc.mdr, inc.mdr.sd)])
out <-
  hr[, {
    tmp = incdr(
      inc,
      inc.sd,
      f,
      f.sd,
      pn = prop.hr.new,
      pn.sd = prop.hr.new.sd,
      r,
      r.sd,
      rho = 1.76,
      rho.sd = .44,
      pr = prop.hr.ret,
      pr.sd = prop.hr.ret.sd
    )$prop
    
    list(inc.hr = tmp[2],
         inc.hr.sd = tmp[4])
  },
  by = .(iso3)]

hr <- merge(hr, out, by = 'iso3', all.x = T)


rr[inc.rr < inc.mdr] # should be empty


#' bounds
#'
m <- 1e5
sel <- hr$inc.hr.sd > hr$inc.hr
table(sel)

(hr[sel, .(iso3, inc.hr, inc.hr.sd, inc.hr / m * (1 - inc.hr /
                                                    m))])

hr[sel, inc.hr.sd := inc.hr * .25]

sel <- hr$inc.hr > 0 & !is.na(hr$inc.hr)
table(sel)
out <- vlohi(hr$inc.hr[sel] / m, hr$inc.hr.sd[sel] / m)
hr[sel, inc.hr.lo := out[1, ] * m]
hr[sel, inc.hr.hi := out[2, ] * m]

hr[inc.hr == 0, inc.hr.lo := 0]
hr[inc.hr == 0, inc.hr.hi := 1.96 * inc.hr.sd]

hr[!is.na(inc.hr), test.bounds(inc.hr, inc.hr.lo, inc.hr.hi)]


whomap(hr[, .(iso3, var = cut(inc.hr, b = c(0, 0.5, 1, 2, 5, 10, 20, Inf)))]) +
  scale_fill_brewer('HR incidence\nper 100 000/year')
ggsave(here('output/checks/hr_incidence_map.pdf'),
       width = 10,
       height = 8)


#' # HR aggregates and Table 3.9 HR incidence
M <- 1e5
m <- 1e3



#' global HR incidence
#'
global.hr <-
  data.table(hr[, addXY(inc.hr / M, r.sd = inc.hr.sd / M, weights = pop)])
#' prop RR
#'
hr[, new := inc * pop * pmax(0, 1 - f - r) / 1e5]
hr[, ret := inc * pop * pmin(f + r, 1) / 1e5]

global.prop <- hr[, .(
  prop.hr.new = weighted.mean(prop.hr.new, w = new),
  prop.hr.new.sd = sqrt(sum((
    prop.hr.new.sd * new /
      sum(new)
  ) ^ 2)),
  prop.hr.ret = weighted.mean(prop.hr.ret, w = ret),
  prop.hr.ret.sd = sqrt(sum((
    prop.hr.ret.sd * new /
      sum(ret)
  ) ^ 2))
)]

fwrite(global.prop, file=here('report/tabs/globalHR.csv'))


#' Table 3.9, H vs R
#'
p <- sum(hr$pop)
global.inc <-
  est[year == 2019, addXY(inc / M, r.sd = inc.sd / M, weights = pop)]
setnames(
  global.inc,
  c(
    "r",
    "r.lo",
    "r.hi",
    "r.sd",
    "r.num",
    "r.lo.num",
    "r.hi.num",
    "pop"
  ),
  c(
    'inc',
    'inc.lo',
    'inc.hi',
    'inc.sd',
    'inc.num',
    'inc.lo.num',
    'inc.hi.num',
    'pop'
  )
)

g <- rrg[8,]

tab <-
  data.table(
    INH = c('HR', 'HS', 'total'),
    RR = c(
      g$inc.rr.num * g$pctMDR / 100,
      g$inc.rr.num * (1 - g$pctMDR / 100),
      g$inc.rr.num
    ),
    RR.sd = NA,
    RR.lo = c(NA, NA, g$inc.rr.num.lo),
    RR.hi = c(NA, NA, g$inc.rr.num.hi),
    RS = NA,
    RS.sd = integer(3L),
    total = c(global.hr$r.num, NA, global.inc$inc.num),
    total.sd = NA,
    total.lo = c(global.hr$r.lo.num, NA, global.inc$inc.lo.num),
    total.hi = c(global.hr$r.hi.num, NA, global.inc$inc.hi.num)
  )
tab$total[2] <- tab$total[3] - tab$total[1]
tab[, total.sd := (total.hi - total.lo) / 3.92]
tab$total.sd[2] <- sqrt(tab$total.sd[3] ^ 2 - tab$total.sd[1] ^ 2)
tab[, RR.sd := (RR.hi - RR.lo) / 3.92]
tab$RR.sd[1] <- tab$RR.sd[3] * g$pctMDR / 100
tab$RR.sd[2] <- tab$RR.sd[3] * (1 - g$pctMDR / 100)
tab[, RS := total - RR]
tab$RS.sd <- sqrt(tab$total.sd ^ 2 - tab$RR.sd[1] ^ 2)

tab[is.na(total.lo), total.lo := as.integer(total - 1.96 * total.sd)]
tab[is.na(total.hi), total.hi := as.integer(total + 1.96 * total.sd)]

tab[is.na(RR.lo), RR.lo := as.integer(RR - 1.96 * RR.sd)]
tab[is.na(RR.hi), RR.hi := as.integer(RR + 1.96 * RR.sd)]

tab[, RS.lo := as.integer(RS - 1.96 * RS.sd)]
tab[, RS.hi := as.integer(RS + 1.96 * RS.sd)]
(ghr <-
    tab[, .(INH, RR, RR.lo, RR.hi, RS, RS.lo, RS.hi, total, total.lo, total.hi)])

save(ghr, file = here('Rdata/ghr.Rdata'))


ghr[, 2:10 := lapply(.SD, ftb), .SDcols = 2:10]
ghr[3, INH := 'Total']
setnames(ghr, c(' ', 'RR', ' ', ' ', 'RS', ' ', ' ', 'Total', ' ', ' '))

fwrite(ghr, file = here('report/tabs/tab4_9_globalDR.csv'))



