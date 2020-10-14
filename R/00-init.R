#' ---
#' title: init GTB2019
#' author: Philippe Glaziou
#' date: 2020-05-30
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---
library(data.table)

setwd('../gtb2020')

load('../gtb2019/Rdata/est.Rdata')
load('../gtb2019/Rdata/cty.Rdata')
load('../gtb2019/Rdata/pop.Rdata')
load('../gtb2019/Rdata/vr.Rdata')
old <- copy(est)
vrcov <- vr[,.(iso3, year, vr.coverage, codqual)]


save(old, file = 'Rdata/old.Rdata')
save(cty, file = 'Rdata/cty.Rdata')
save(pop, file = 'Rdata/pop.Rdata')
save(vrcov, file = 'Rdata/vrcov.Rdata')
