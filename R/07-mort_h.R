#' ---
#' title: TB mortality HIV-positive
#' author: Philippe Glaziou
#' date: 09/07/2020
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
#' TB mortality HIV-positive
#'
#' # Load libraries and data
#'
library(data.table)
library(here)

yr <- 2019
m <- 1e5
source(here('R/fun.R'))

load(here('Rdata/est2.Rdata'))
load(here('Rdata/cty.Rdata'))
load(here('Rdata/old.Rdata'))
load(here('Rdata/pop.Rdata'))
load(here('Rdata/grpmbr.Rdata'))
load(here('Rdata/unaids.Rdata'))


#' vectorize lohi
#'
vlohi <- Vectorize(lohi, c('ev', 'sd'))


#' # Indirect mortality HIV-positive
#'
est2$art.coverage <-
  est2$hiv.art.p / (est2$c.notified * est2$tbhiv) # ART in treated TB

sel <- est2$year >= 2015  # change in vars
est2$art.coverage[sel] <-
  est2$newrel.art[sel] / (est2$c.notified[sel] * est2$tbhiv[sel])
sel <- is.infinite(est2$art.coverage)
table(sel)
est2$art.coverage[sel] <- NA
est2[!is.na(art.coverage), test.ispos(art.coverage)]

sel <- is.na(est2$art.coverage) & est2$g.income %in% c('HIC', 'UMC')
table(sel)
est2$art.coverage[sel] <- 1

sel <- is.na(est2$art.coverage)
table(sel)
(unique(as.character(est2$iso3[sel])))
est2$art.coverage[sel] <- 0
summary(est2$art.coverage)
sel <- est2$art.coverage > 1 & !is.na(est2$art.coverage)
table(sel)
est2$art.coverage[sel] <- 1

out <-
  est2[, {
    tmp = inc2mort(inc,
                   inc.sd,
                   imp.newinc / (1 - U),
                   tbhiv,
                   tbhiv.sd,
                   art.coverage,
                   noHIV = FALSE)$prop
    
    list(e.mort.h = tmp[2],
         e.mort.h.sd = tmp[4])
  },
  by = .(iso3, year)] # takes a while
out[, test.isbinom(e.mort.h / m)]
out[, test.isbinom(e.mort.h.sd / m)]

est3 <- merge(est2, out, by = c('iso3', 'year'), all.x = T)
setnames(est3, c('e.mort.h', 'e.mort.h.sd'), c('mort.h', 'mort.h.sd'))


#' # Inconsistencies with HIV mortality
#'
sel <-
  est3$mort.h > est3$mort.hiv &
  !is.na(est3$mort.hiv) & !is.na(est3$mort.h)
table(sel)
lst <- unique(est3$iso3[sel])

exclude <- c('LSO', 'LAO')
sel2 <- sel & est3$iso3 %ni% exclude
est3$mort.h[sel2] <-
  est3$mort.hiv[sel2] * .5 # arbitrary high proportion

sum(is.na(est3$mort.h))





#' bounds
#'
sel <- est3$mort.h > 0 & est3$mort.h.sd > 0
table(sel)
out <- with(est3[sel], vlohi(mort.h / m, mort.h.sd / m))

est3$mort.h.lo[sel] <- out[1, ] * m
est3$mort.h.hi[sel] <- out[2, ] * m

sel <- est3$mort.h.sd == 0 & !is.na(est3$mort.h.sd)
table(sel)

sel <-
  (est3$mort.h.lo > est3$mort.h) |
  (est3$mort.h.hi < est3$mort.h) &
  (!is.na(est3$mort.h) &
     !is.na(est3$mort.h.lo) & !is.na(est3$mort.h.hi))
table(sel)
# est3[sel, .(iso3, year, inc, mort.h, mort.h.sd, mort.h.lo, mort.h.hi)]
# est3[sel, mort.h.hi := mort.h + 1.96 * mort.h.sd]

sel <- est3$mort.h == 0 & est3$mort.h.sd > 0
table(sel)
est3[sel, mort.h.lo := 0]
est3[sel, mort.h.hi := 1.96 * mort.h.sd]

sel <- est3$mort.h.sd > est3$mort.h
table(sel)
est3[sel, mort.h.sd := mort.h * 0.25]

est3[, test.bounds(mort.h, mort.h.lo, mort.h.hi)]

est3 <- within(est3, {
  mort.h.num <- mort.h * pop / m
  mort.h.lo.num <- mort.h.lo * pop / m
  mort.h.hi.num <- mort.h.hi * pop / m
})


est3[, .(sum(mort.h.num),
         sum(inc.h.num, na.rm = T),
         sum(mort.h.num) / sum(inc.h.num, na.rm = T)), by = year]
old[, .(sum(mort.h.num),
        sum(inc.h.num, na.rm = T),
        sum(mort.h.num) / sum(inc.h.num, na.rm = T)), by = year]





#' plots
#'
wr <- unique(as.character(est3$g.whoregion))
comp <-
  merge(est3[, . (iso3, year, g.whoregion, mort.h, mort.h.lo, mort.h.hi)],
        old[year > 1999, .(
          iso3,
          year,
          omort.h = mort.h,
          omort.h.lo = mort.h.lo,
          omort.h.hi = mort.h.hi
        )],
        by = c('iso3', 'year'), all.x = T)
dim(comp)

#+ fig.width=16, fig.height=12
for (i in wr) {
  p <-
    qplot(
      year,
      mort.h,
      data = subset(comp, g.whoregion == i),
      geom = 'line',
      colour = I('blue')
    ) +
    geom_line(aes(year, omort.h),
              colour = I('red'),
              linetype = I(2)) +
    geom_ribbon(
      aes(year, ymin = mort.h.lo, ymax = mort.h.hi),
      fill = I('blue'),
      alpha = I(0.3)
    ) +
    facet_wrap(~ iso3, scales = 'free_y') + xlab('') + ylab('TB deaths')
  suppressWarnings(print(p))
  suppressWarnings(ggsave(here(
    paste('output/checks/mort.h_', i, '_compare.pdf', sep = '')),
    width = 14,
    height = 8
  ))
}




#' # Total TB mortality
#'
est3$mort <- rowSums(cbind(est3$mort.nh, est3$mort.h), na.rm = T)
est3[, list(sum(mort * pop / m, na.rm = TRUE)), by = year] # looks ok

est3$mort.sd <-
  sqrt(rowSums(cbind(est3$mort.nh.sd ^ 2, est3$mort.h.sd ^ 2), na.rm = T))

summary(est3$mort)
summary(est3$mort.sd)
est3[, test.ispos(mort)]
est3[, test.ispos(mort.sd)]

sel <- est3$mort.sd > est3$mort & est3$mort > 0
table(sel)


sel <- est3$mort > 0 & est3$mort.sd > 0
out <- with(est3[sel], vlohi(mort / m, mort.sd / m))
est3$mort.hi[sel] <- out[2, ] * m
est3$mort.lo[sel] <- out[1, ] * m

sel <- est3$mort.sd == 0 & !is.na(est3$mort.sd)
table(sel)
est3$mort.lo[sel] <- est3$mort[sel]
est3$mort.hi[sel] <- est3$mort[sel]

# sel <- (est3$mort.lo > est3$mort) | (est3$mort.hi < est3$mort) & (!is.na(est3$mort) & !is.na(est3$mort.lo) & !is.na(est3$mort.hi))
# table(sel)
# est3[sel,.(iso3,year,inc,mort,mort.sd,mort.lo,mort.hi)]
# est3[sel,mort.hi:=mort + 1.96 * mort.sd]

sel <-
  est3$mort.sd > 0 &
  !is.na(est3$mort.sd) & est3$mort == 0 & !is.na(est3$mort)
table(sel)
est3$mort.lo[sel] <- 0
est3$mort.hi[sel] <- est3$mort.sd[sel] * 1.96


est3[, test.bounds(mort, mort.lo, mort.hi)]



#' absolute numbers
#'
est3 <- within(est3, {
  mort.num <- mort * pop / m
  mort.lo.num <- mort.lo * pop / m
  mort.hi.num <- mort.hi * pop / m
})



#' Save
#'
save(est3, file = here('Rdata/est3.Rdata'))
fwrite(est3, file = here(paste0('csv/est_s07_', Sys.Date(), '.csv')))
