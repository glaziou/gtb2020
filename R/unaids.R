#' ---
#' title: UNAIDS estimates
#' author: Philippe Glaziou
#' date: 2020-06-10
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

library(data.table)
library(bit64)
library(readxl)

load('Rdata/cty.Rdata')

#' vectorize lohi
#'
source('R/fun.R')
vlohi <- Vectorize(lohi, c('ev', 'sd'))
m <- 1e5

dic <- fread('input/UNAIDS/dic.csv')

# typz <- c(rep('text', 5),
#           rep('numeric', 66))
# rangez <- 4:(4+169)
# l <- 4 + 174 * 1:19
# for (i in l) rangez <- c(rangez, i:(i+169)) 
# rm(l)
# 
# suppressWarnings(unaids <-
#                    as.data.table(
#                      read_excel(
#                        'input/UNAIDS/HIV2020Estimates_ForTB_9June2020.xlsx',
#                        sheet = 1,
#                        skip = 1,
#                        col_types = typz
#                      )
#                    ))

# unaids <- fread('input/UNAIDS/HIV2020Estimates_ForTB_9June2020.csv')
unaids <- fread('input/UNAIDS/HIV2020Estimates_ForTB_17June2020.csv')
# unaids[, UN := NULL]

setnames(unaids, dic$varnames)
# unaids <- unaids[rangez, ]

setkey(unaids, iso3, year)
unaids[, country := NULL]

unaids[, pop014 := pop04 + pop514]
unaids[, pop := pop014 + pop15 + pop65]
unaids[, hiv.num := h014 + h15]

unaids <- within(unaids, {
  hiv15 <- h15 / pop15
  hiv15.lo <- h15.lo / pop15
  hiv15.hi <- h15.hi / pop15
  hiv014 <- h014 / pop
  hiv014.lo <- h014.lo / pop
  hiv014.hi <- h014.hi / pop
  inchiv <- (inch014 + inch15) / pop
  inchiv.num <- inch014 + inch15
  d <- rowSums(cbind(d15, d014), na.rm = T)
})

unaids <- within(unaids, {
  hiv <- hiv.num / pop
  hiv.sd <-
    sqrt(((h014.hi - h014.lo) / 3.92) ^ 2 + ((h15.hi - h15.lo) /
                                               3.92) ^ 2) / pop
  inchiv.sd <-
    sqrt(((inch014.hi - inch014.lo) / 3.92) ^ 2 + ((inch15.hi - inch15.lo) /
                                                     3.92) ^ 2) / pop
  mort.hiv <-
    d / pop * m
  mort.hiv.sd <-
    sqrt(((d15.hi - d15.lo) / 3.92) ^ 2 + ((d014.hi - d014.lo) /
                                             3.92) ^ 2) / pop * m
})

sel <-
  unaids$hiv > 0 &
  !is.na(unaids$hiv) &
  unaids$hiv.sd > 0 &
  unaids$hiv.sd^2 < unaids$hiv * (1 - unaids$hiv)
(table(sel))
out <- vlohi(unaids$hiv[sel], unaids$hiv.sd[sel])
unaids[sel, hiv.lo := out[1,]]
unaids[sel, hiv.hi := out[2,]]
unaids[hiv>hiv.hi, hiv.hi := hiv + 1.96 * hiv.sd]
unaids[sel, test.bounds(hiv, hiv.lo, hiv.hi)]


sel <-
  unaids$mort.hiv > 0 &
  unaids$mort.hiv.sd > 0 & !is.na(unaids$mort.hiv)
(table(sel))
out <- vlohi(unaids$mort.hiv[sel] / m, unaids$mort.hiv.sd[sel] / m)
unaids$mort.hiv.lo[sel] <- out[1,] * m
unaids$mort.hiv.hi[sel] <- out[2,] * m
unaids$mort.hiv.lo[!sel] <- 0
unaids$mort.hiv.hi[!sel] <- 0

unaids[sel, test.bounds(mort.hiv, mort.hiv.lo, mort.hiv.hi)]

unaids <- within(unaids, {
  hiv.lo.num <- hiv.lo * pop
  hiv.hi.num <- hiv.hi * pop
  mort.hiv.num <- d
  mort.hiv.lo.num <- mort.hiv.lo * pop / m
  mort.hiv.hi.num <- mort.hiv.hi * pop / m
})


#' check aggregates
#' 
(unaids[, sums(hiv.num), by=year])


#' save
#'
save(unaids, file = 'Rdata/unaids.Rdata')
fwrite(unaids, file = paste0('csv/unaids_', Sys.Date(), '.csv'))
