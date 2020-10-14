#' title: RR risk in relapse
#' author: Philippe Glaziou
#' date: 17/04/2020 
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

#' Last updated: 2020-04-17
#'
#'
#'
#' # Load libraries and data
#'
library(metafor)
load('Rdata/old.Rdata')

rel <- fread('input/dr/new_vs_relapse_3Apr2020.csv')
setkey(rel, iso3)
rel <- merge(rel, old[year==2018,.(iso3, g.mdr)], all.x=T, all.y=F)

out1 <- rel[, {
  tmp = cii(dst_rlt_new, mdr_rr_new)
  
  list(
    prop.rr.new = tmp$prob,
    prop.rr.new.sd = tmp$se,
    prop.rr.new.lo = tmp$lower95ci,
    prop.rr.new.hi = tmp$upper95ci
  )
}, by = .(iso3, level, year)]
rel <- merge(rel, out1, by = c('iso3','level','year'), all.x = T)


out2 <- rel[, {
  tmp = cii(dst_rlt_rel, mdr_rr_rel)
  
  list(
    prop.rr.rel = tmp$prob,
    prop.rr.rel.sd = tmp$se,
    prop.rr.rel.lo = tmp$lower95ci,
    prop.rr.rel.hi = tmp$upper95ci
  )
}, by = .(iso3, level, year)]
rel <- merge(rel, out2, by = c('iso3','level','year'), all.x = T)

out3 <- rel[, {
  tmp = divXY(prop.rr.rel, prop.rr.new, prop.rr.rel.sd^2, prop.rr.new.sd^2)
  
  list(
    ratio.rel.new = tmp[[1]],
    ratio.rel.new.sd = sqrt(tmp[[2]])
  )
}, by = .(iso3, level, year)]
rel <- merge(rel, out3, by = c('iso3','level','year'), all.x = T)

out4 <- rel[, {
  tmp = cii(dst_rlt_ret, mdr_rr_ret)
  
  list(
    prop.rr.ret = tmp$prob,
    prop.rr.ret.sd = tmp$se,
    prop.rr.ret.lo = tmp$lower95ci,
    prop.rr.ret.hi = tmp$upper95ci
  )
}, by = .(iso3, level, year)]
rel <- merge(rel, out4, by = c('iso3','level','year'), all.x = T)
rel[region=='SEAR', region := 'SEA']



suppressWarnings(fit.relnew1 <-
                   rma(
                     yi = ratio.rel.new,
                     sei = ratio.rel.new.sd,
                     data = rel[!is.nan(ratio.rel.new) & ratio.rel.new > 0],
                     mods =  ~ region - 1
                   ))

(summary(fit.relnew1)) 


suppressWarnings(fit.relnew2 <-
                   rma(
                     yi = ratio.rel.new,
                     sei = ratio.rel.new.sd,
                     data = rel[!is.nan(ratio.rel.new) & ratio.rel.new > 0],
                     mods =  ~ g.mdr - 1
                   ))

(summary(fit.relnew2)) 


#' lots of heterogeneity even after accounting for regions
#' let us try something else
#' 
rel[, dst.retnr := dst_rlt_ret - dst_rlt_rel]
rel[, mdr.rr.retnr := mdr_rr_ret - mdr_rr_rel]
rel[dst.retnr > 0, prop.rr.retnr := mdr.rr.retnr / dst.retnr]
rel[, prop.rel := dst_rlt_rel / dst_rlt_ret]
(rel[, .(prop.rr.new, prop.rr.ret, prop.rr.rel, prop.rr.retnr, prop.rel)][])

fit2 <- lm(ratio.rel.new ~ I(prop.rr.ret/prop.rr.new) + region, data=rel)
(summary(fit2))
fit3 <- lm(ratio.rel.new ~ I(prop.rr.ret/prop.rr.new), data=rel)
(summary(fit3))


suppressWarnings(fit.relnew0 <-
                   rma(
                     yi = ratio.rel.new,
                     sei = ratio.rel.new.sd,
                     data = rel[!is.nan(ratio.rel.new) & ratio.rel.new > 0],
                   ))

(summary(fit.relnew0)) 

save(rel, file='Rdata/rel.Rdata')
