#' ---
#' title: VR data
#' author: Philippe Glaziou
#' date: 2020/05/30
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
#'
#' # Load libraries and data
#'

#' Load libraries and data
library(data.table)
library(imputeTS)
library(haven) # read_dta
library(readxl)
library(stringr)
library(here)

source(here('R/fun.R'))


#' load data
#'
load(here('Rdata/old.Rdata'))
load(here('Rdata/cty.Rdata'))
load(here('Rdata/pop.Rdata'))
load(here('Rdata/vrcov.Rdata'))
m <- 1e5


#' convert VR data in excel
#' adapting code from Pete (2019)
#'
#' function for reformatting
refrm <- function(indat) {
  # rename & aggregate
  indat <- indat[, .(
    Country,
    name,
    Year,
    icd,
    Cause,
    cause1,
    Sex,
    `0-4` = rowSums(cbind(
      Deaths2, Deaths3, Deaths4, Deaths5, Deaths6
    ), na.rm = T),
    `5-14` = rowSums(cbind(Deaths7, Deaths8), na.rm = T),
    `15-24` = rowSums(cbind(Deaths9, Deaths10), na.rm = T),
    `25-34` = rowSums(cbind(Deaths11, Deaths12), na.rm = T),
    `35-44` = rowSums(cbind(Deaths13, Deaths14), na.rm = T),
    `45-54` = rowSums(cbind(Deaths15, Deaths16), na.rm = T),
    `55-64` = rowSums(cbind(Deaths17, Deaths18), na.rm = T),
    `65plus` = rowSums(
      cbind(
        Deaths19,
        Deaths20,
        Deaths21,
        Deaths22,
        Deaths23,
        Deaths24,
        Deaths25
      ),na.rm = T))
    ]
  
  # Sequelae
  seq <- c("B90", "B900", "B901", "B902", "B908", "B909", "B077")
  indat[Cause %in% seq, cause1 := "tbseq"]
  indat[, Cause := NULL]
  
  ## reshape
  MM <-
    melt(indat, id = c("Country", "name", "Year", "icd", "cause1", "Sex"))
  MM$sex <- c('M', 'F')[as.numeric(MM$Sex)]
  MM[is.na(sex), sex := 'U']
  MM$sex <- factor(MM$sex)
  MM$age <- factor(MM$variable, levels = agz3, ordered = TRUE)
  MM[, age_group := gsub('-', '_', age)]
  MM[, age := NULL]
  MM
}

#' Some useful age range vectors:
agz <-
  c('04', '514', '1524', '2534', '3544', '4554', '5564', '65') #for extract
agz2 <-
  c('0_4',
    '5_14',
    '15_24',
    '25_34',
    '35_44',
    '45_54',
    '55_64',
    '65plus') #for labels
agz3 <- gsub('_', '-', agz2)
agz4 <- c(rev(rev(agz3)[-1]), "\u2265 65")
kds <- c('0_4', '5_14')
kdz <- c('04', '514')
AA <-
  data.table(
    a1 = agz,
    age_group = agz2,
    age = agz3,
    agegt = agz4
  ) #for conversion


typz <- c('text',
          'text',
          'text',
          'text',
          rep('text', 5),
          rep('numeric', 26)) #specify column types to avoid warning
xlfn <- 'input/VR/info_TB_may20.xls'



## -------- ICD>10.1
M1 <-
  as.data.table(read_excel(
    xlfn,
    sheet = 4,
    skip = 1,
    col_types = typz
  )) 
M1 <- refrm(M1)


## -------- ICD 10.1
M2 <-
  as.data.table(read_excel(
    xlfn,
    sheet = 5,
    skip = 0,
    col_types = typz
  )) 
M2 <- refrm(M2)


## -------- ICD 9
M3 <-
  as.data.table(read_excel(
    xlfn,
    sheet = 6,
    skip = 0,
    col_types = typz
  ))
# remove double counted B02
M3[, dble := sums(str_count(Cause, 'B02[0-9]')), by = .(name, Year, Sex)]
M3[Cause == 'B02' & dble > 0, drop := TRUE]
(M3[, table(drop)])
M3 <- M3[is.na(drop)]
M3[, drop := NULL]
M3[, dble := NULL]
M3 <- refrm(M3)


## -------- ICD 8
M4 <-
  as.data.table(read_excel(
    xlfn,
    sheet = 7,
    skip = 0,
    col_types = typz
  )) 
M4 <- refrm(M4)


## --- join ---
VR <- rbind(M1, M2, M3, M4)


## Differences in names to be done by hand:
(vrbad <- setdiff(VR[, unique(name)],
                  pop[, unique(country)]))

## direct renamings:
done <- c()
for (i in seq_along(vrbad)) {
  newnm <- grep(vrbad[i], pop[, unique(country)], value = TRUE)
  if (length(newnm) > 0) {
    print(newnm)
    VR[name == vrbad[i], name := newnm]
    done <- c(done, i)
  }
}
vrbad <- vrbad[-done]

## others for renaming
vrbad
(newnm <- grep('Czech', pop[, unique(country)], value = TRUE))
VR[name == grep('Czech', vrbad, value = TRUE), name := newnm]
(newnm <- grep('Serbia', pop[, unique(country)], value = TRUE)[1])
VR[name == grep('Serb', vrbad, value = TRUE), name := newnm]
(newnm <-
    grep('Macedonia', pop[, unique(country)], value = TRUE)[1])
VR[name == grep('Mace', vrbad, value = TRUE), name := newnm]
(newnm <- grep('Vincent', pop[, unique(country)], value = TRUE)[1])
VR[name == grep('Vincent', vrbad, value = TRUE), name := newnm]

## those still bad
(vrbad <- vrbad[!str_detect(vrbad, 'Cze|Serb|Mace')])

## sub-countries
VR[name %in% c(
  "French Guiana",
  "Martinique",
  "Reunion",
  "Mayotte",
  "Guadeloupe",
  "Saint Pierre and Miquelon"
), name := "France"]
VR[name %in% c("Rodrigues"), name := "Mauritius"]
VR[name %in% c("Occupied Palestinian Territory"), name := "West Bank and Gaza Strip"]

## check
(dropname <- setdiff(VR[, unique(name)], cty[, unique(country)]))
VR <- VR[name %ni% dropname]

VR[, year := as.integer(Year)]
VR[, Year := NULL]

## Add iso3, and tidy up
VR <- merge(VR,
            cty[, .(iso3, country)],
            by.x = 'name',
            by.y = "country",
            all.x = TRUE)
(VR[is.na(iso3), unique(name)])

VR[, Sex := NULL]
VR[, age_group := variable]
VR$age_group <-
  factor(gsub("-", "_", VR$age_group),
         levels = agz2,
         ordered = TRUE)
VR[, variable := NULL]

setkey(VR, iso3)
rm(M1, M2, M3, M4)

save(VR, file = 'Rdata/VR.Rdata')



#' aggregate and reshape long to wide
#'
vr <- VR[, .(value = sums(value)),
         by = .(iso3, year, cause1)]

vr <- dcast(vr, ... ~ cause1)
setkey(vr, iso3)
setnames(vr, c('iso3','year','hiv','tb','ill_def','tb_seq','total'))


#' import coverage and coding quality
#' 
vr <- merge(vr, vrcov[,.(iso3, year, codqual)], by=c('iso3','year'), all.x = TRUE)



#' # Country data additions to WHO database
#'

#' * Azerbaijan additions
#'
#' source: epi review, the Hague 2017
#'
addAZE <- vr['AZE']
addAZE$year <- 2010:2015
dim(vr)
vr2 <- rbind(vr, addAZE, use.names = TRUE)
setkey(vr2, iso3)
sel <- vr2$iso3 == 'AZE' & vr2$year %in% 2010:2015

# vr2$vr.coverage[sel] <- rep(1, 6)
vr2$total[sel] <- c(53580, 53726, 55017, 54383, 55648, 54697)
vr2$tb[sel] <- c(709, 577, 373, 378, 372, 485)
vr2$ill_def[sel] <- c(1343, 1771, 1836, 1892, 2440, 1864)
vr2['AZE']


#' RUS 2019
#' 
#' Total deaths = 1798307
#' TB deaths = 7536
#' Ill-defined = 124940
addRUS <- vr2['RUS'][year==2018]
addRUS$year <- 2019
addRUS$total <- 1798307
addRUS$tb <- 7536
addRUS$ill_def <- 124940
vr3 <- rbind(vr2, addRUS, use.names = TRUE)
setkey(vr3, iso3)

vr <- copy(vr3)
(dim(vr))



#' check that the additions do not end-up duplicating country-year entries
#' 
sum(duplicated(vr[,.(iso3,year)]))==0


#' # process raw VR data
#'
(dim(vr))
# miss.est <- unique(vr$iso3[is.na(vr$vr.coverage)])

# all in small countries, plus HKG, PRK, MAC, part of ZAF series, ZWE 1990
# we assume vrcoverage=1 except for ZAF and ZWE
# 
# vr[iso3 %in% miss.est, vrcov.nm := sum(!is.na(vr.coverage)), by = iso3]
# vr[, table(vrcov.nm)]
# vr[is.na(vr.coverage) & vrcov.nm < 2, vr.coverage := 1]
# vr[vrcov.nm >= 2, vr.coverage := na_interpolation(vr.coverage), by = iso3]
vr[, codqual.nm := sum(!is.na(codqual)), by = iso3]
vr[, table(codqual.nm)]
vr[codqual.nm>= 2, codqual := na_interpolation(codqual), by = iso3]
# (summary(vr$vr.coverage))
(vr['ZAF'])
(vr['RUS'])
#' clean up
#'
# vr[, vrcov.nm := NULL]
vr[, codqual.nm := NULL]



#' missing tb deaths
miss.tb <- unique(vr$iso3[is.na(vr$tb)])

vr[iso3 %in% miss.tb, tb.nm := sum(!is.na(tb)), by = iso3]
vr[, table(tb.nm)]
vr[tb.nm == 0]

(dim(vr))
vr <- vr[tb.nm > 0 | iso3 %ni% miss.tb]
(dim(vr))

vr[tb.nm == 1] # not yet imputable

vr[, prop.tb := tb / total]
vr[tb.nm >= 2, tb.imp := na_interpolation(round(prop.tb * total)), by =
     iso3]
(summary(vr$tb.imp))
vr[!is.na(tb) & is.na(tb.imp), tb.imp := tb]
(vr['ZAF'])
(vr['RUS'])

#' clean up
#'
vr[, tb.nm := NULL]
rm(miss.tb)



#' missing ill-defined
miss.ill <- unique(vr$iso3[is.na(vr$ill_def)])

vr[iso3 %in% miss.ill, ill.nm := sum(!is.na(ill_def)), by = iso3]
(vr[, table(ill.nm)])

vr[, garbage := ill_def / total]

vr[ill.nm >= 2, garbage.imp := na_interpolation(garbage), by = iso3]
(summary(vr$garbage.imp))
vr[!is.na(garbage) & is.na(garbage.imp), garbage.imp := garbage]
(vr['ZAF'])
(vr['RUS'])



#' clean up
#'
vr[, ill.nm := NULL]
rm(miss.ill)


#' missing tb_seq
#'
vr[, sum(is.na(tb_seq))]


#' total TB deaths
#'
vr[, totaltb := rowSums(cbind(vr$tb.imp, vr$tb_seq), na.rm = TRUE)]

#' proportion sequelae
#'
vr[, seq.prop := tb_seq / totaltb]


setkey(vr, iso3, year)





#' # GHO envelopes
#'
#'
ghe <- fread('input/GHE/ghe2016_deaths_country_allages.csv')
ghe <- subset(ghe, cause2015 == 0 & sex == 'BTSX')

ghe2 <- melt(ghe[, -c(2:6), with = F], id.vars = 1)
ghe2[, year := as.integer(gsub('\\D', '', variable))]
ghe2[, variable := gsub('\\d', '', variable)]
ghe3 <- dcast(ghe2, iso3 + year ~ variable, value.var = 'value')

ghe <- copy(ghe3)
setkey(ghe, iso3, year)
setnames(ghe,
         c('dths', 'low', 'upp'),
         c('ghe.env', 'ghe.env.lo', 'ghe.env.hi'))

save(ghe, file = 'Rdata/ghe.Rdata')
fwrite(ghe, file = paste0('csv/ghe_', Sys.Date(), '.csv'))





#' Add GHE envelopes in the vr file
#' -- following a call on 13 July 2017 with Doris Mafat
#' explaining why the totals differ with GHE
#' GHE is the WHO reference to use
vr2 <- merge(
  vr,
  ghe,
  by = c('iso3', 'year'),
  all.x = T,
  all.y = F
)
dim(vr)
dim(vr2)

vr2[, env := ghe.env]                    # default envelope
vr2[iso3 %in% c('HKG','DMA','PRI','SMR','KNA'), env := total]

vr2[, env.nm := sum(!is.na(env)), by=iso3]
vr2[, table(env.nm)]
vr2[env.nm >=2, env := na_interpolation(env), by=iso3] 
vr2[, sum(is.na(env))]


vr2[, vr.coverage := pmin(1, total / env)] # VR coverage
vr2[is.na(vr.coverage), vr.coverage := total] # non-GHE estimate where missing GHE
# 
# env needs to be consistent in 2017 (countries with 2017 data)
# lst <- vr2[year == 2017, iso3]
# vr2[iso3 %in% lst, ghe.env.nm := sum(!is.na(ghe.env)), by = iso3]
# vr2[ghe.env.nm >= 2 & year > 2014]
# vr2[ghe.env.nm >= 2, env := na_kalman(ghe.env), by = iso3]

# qplot(year, env, data = vr2[ghe.env.nm >= 2], geom = 'line') +
#   geom_point(aes(year, ghe.env), colour = I('red')) +
#   facet_wrap( ~ iso3, scales = 'free_y')



#' proportion of TB deaths out of well documented deaths
#'
vr2[, tb.prop := totaltb / (total - ill_def)]


#'
#' adjusted TB deaths
#'
vr2[, tb.adj := env * tb.prop]

vr <- copy(vr2)

#' check
#'
vr[!is.na(tb.prop), test.isbinom(tb.prop)]





#' # VR quality
#' COD quality (WHS 2018, SDG target 17.19, p57)
#' 1 = high
#' 2 = medium
#' 3 = low
#' 4 = very low
#'
#'
vr[, .(mean(vr.coverage, na.rm = T),
       min(vr.coverage, na.rm = T),
       max(vr.coverage, na.rm = T)), by = codqual][order(codqual)]
vr[, .(mean(garbage, na.rm = T),
       min(garbage, na.rm = T),
       max(garbage, na.rm = T)), by = codqual][order(codqual)]

vr[codqual %in% 1:3, keep.vr := T]
vr[codqual %in% 4, keep.vr := F]
vr[, summary(keep.vr)]

vr[is.na(keep.vr), summary(vr.coverage)]
vr[is.na(keep.vr), summary(garbage)]

vr[is.na(keep.vr) & !is.na(vr.coverage), sum(vr.coverage > .8)]
vr[is.na(keep.vr), sum(garbage < .2, na.rm = T)]
vr[is.na(keep.vr), sum(vr.coverage > .8 & garbage < .2, na.rm = T)]
vr[is.na(keep.vr), keep.vr := vr.coverage > .8 & garbage < .2]
vr[is.na(keep.vr), keep.vr := F]


vr[, summary(keep.vr)]






#' # SDs
#'
#' assume TB deaths between 0.5 and 1.5 times observed $t$ rate among
#' garbage $g$ and non covered $c$
#'
#' $t_{adj} = \frac{t}{c(1-g)}$
#'
#' $\text{SD}(t_{adj}) = \frac{t}{4} \left(\frac{1}{c(1-g)} - 1\right)$
#'
vr[, tb.adj.sd := (tb.adj / 4) * (1 / (vr.coverage * (1 - garbage)) - 1)]

vr[!is.na(tb.adj.sd), test.ispos(tb.adj.sd)]
vr[keep.vr == T & is.na(tb.adj.sd), tb.adj.sd := tb.adj * .2]
vr[keep.vr == T, summary(tb.adj.sd / tb.adj)]


vr2  <-
  merge(vr[iso3 %ni% 'VIR'], pop[, .(iso3, year, pop = e.pop.num)], by =
          c('iso3', 'year'), all.x = TRUE)
dim(vr[iso3 %ni% 'VIR'])
dim(vr2)
vr2[is.na(pop)]

#' exclude MNE and SRB 2000:2004
vr <- copy(vr2[!is.na(pop)])




#'  save vr dataset
#'
save(vr, file = 'Rdata/vr.Rdata')
fwrite(vr, file = paste0('csv/vr_', Sys.Date(), '.csv'))




#' Codecorrect from IHME
#'
load('Rdata/ihme.Rdata')   # see ihme.R


#' envelope ratios GHO/IHME
#'
ihme.all <-
  ihme[year >= 2000 &
         metric_name == 'Number' &
         cause_name == 'All causes' & measure_name == 'Deaths']
ihme.all <- ihme.all[!is.na(iso3)]
gbd <-
  merge(ihme.all[, .(iso3,
                     year,
                     ihme = val,
                     ihme.lo = lower,
                     ihme.hi = upper)],
        ghe,
        by = c('iso3', 'year'),
        all.y = T)


#' missing IHME envelopes?
#'
gbd[, sapply(.SD, function(x)
  sum(is.na(x)))]


#' WHO/IHME env ratio
#'
gbd[, env.ratio := ghe.env / ihme]


#' IHME estimates of TB deaths
#'
ihme.tb <-
  ihme[year >= 2000 &
         metric_name == 'Number' &
         cause_name == 'Tuberculosis' &
         measure_name == 'Deaths' & !is.na(iso3)]


#' save
#'
save(gbd, file = here('Rdata/gbd.Rdata'))
save(ihme.tb, file = here('Rdata/ihme.tb.Rdata'))
save(ihme.all, file = here('Rdata/ihme.all.Rdata'))

fwrite(gbd, file = paste0(here('csv/gbd_'), Sys.Date(), '.csv'))
fwrite(ihme.tb, file = paste0(here('csv/ihme.tb_'), Sys.Date(), '.csv'))
fwrite(ihme.all, file = paste0(here('csv/ihme.all_'), Sys.Date(), '.csv'))

