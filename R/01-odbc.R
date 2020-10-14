#' ---
#' title: GTB views
#' author: Philippe Glaziou
#' date: 2020-06-29
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
#' Download GTB views
#' 

library(RODBC)
library(data.table)
library(here)

rm(list = ls())


source(here('R/fun.R'))
yr <- 2019

#' # Create connection string
#' 
connection_string <-
  "driver={SQL Server}; server=ssdb231.who.int; database=TMEData; uid=TMEDATA_reader; pwd=Tmere@d3r"

#' Connect to the database
#' 
channel <- odbcDriverConnect(connection_string)

pop <-
  sqlQuery(channel,
           paste("select * from view_TME_estimates_population"),
           na.strings = "")
tb <-
  sqlQuery(channel,
           paste("select * from view_TME_master_notification"),
           na.strings = "")
cty <-
  sqlQuery(channel,
           paste("select * from view_TME_master_report_country"),
           na.strings = "")
sty <-
  sqlQuery(channel,
           paste("select * from view_TME_master_strategy"),
           na.strings = "")
tx <-
  sqlQuery(channel,
           paste("select * from view_TME_master_outcomes"),
           na.strings = "")
dic <-
  sqlQuery(channel,
           paste("select * from view_TME_data_dictionary"),
           na.strings = "")
codes <-
  sqlQuery(channel,
           paste("select * from view_TME_data_codes"),
           na.strings = "")
drs <-
  sqlQuery(channel,
           paste("select * from view_TME_master_drs"),
           na.strings = "") # everything
dre <-
  sqlQuery(channel,
           paste("select * from view_DRS_most_recent_for_estimation"),
           na.strings = "") # curated, most recent
drnew <-
  sqlQuery(channel,
           paste("select * from view_DRS_for_estimation_new"),
           na.strings = "") # curated new, all years
drret <-
  sqlQuery(channel,
           paste("select * from view_DRS_for_estimation_ret"),
           na.strings = "") # curated ret, all years
drh <-
  sqlQuery(channel,
           paste("select * from view_TME_estimates_drtb_rawvalues"),
           na.strings = "") # history of estimates >= 2015
drh2 <-
  sqlQuery(channel,
           paste("select * from view_TME_estimates_drtb"),
           na.strings = "") # history of estimates < 2015

hr <-
  sqlQuery(channel,
           paste("select * from view_DRS_most_recent_for_estimation_INH"),
           na.strings = "") # history of estimates < 2015

dr <-
  sqlQuery(channel,
           paste("select * from view_TME_master_dr_surveillance"),
           na.strings = "")

ntp <-
  sqlQuery(
    channel,
    paste(
      "select iso3, ntp_name, ntp_email from view_TME_master_data_collection"
    ),
    na.strings = ""
  )
gp <-
  sqlQuery(channel,
           paste("select * from view_TME_estimates_global_plan"),
           na.strings = "")
dots <-
  sqlQuery(
    channel,
    paste("select * from view_TME_master_notification_dots_ndots"),
    na.strings = ""
  )
agg <-
  sqlQuery(channel,
           paste("select * from view_TME_master_TBHIV_for_aggregates"),
           na.strings = "")
mdrN <-
  sqlQuery(channel,
           paste("select * from view_TME_estimates_mdr_in_notified"),
           na.strings = "")
prev <-
  sqlQuery(
    channel,
    paste("select * from survey.view_prevalence_survey_estimates"),
    na.strings = ""
  )

dr.est <-
  sqlQuery(channel,
           paste("select * from view_TME_estimates_drtb"),
           na.strings = "")
aggdr <-
  sqlQuery(channel,
           paste("select * from view_TME_aggregated_estimates_drtb"),
           na.strings = "")

grptypes <-
  sqlQuery(channel,
           paste("select * from view_country_group_types"),
           na.strings = "")
grp <-
  sqlQuery(channel,
           paste("select * from view_country_groups"),
           na.strings = "")
grpmbr <-
  sqlQuery(channel,
           paste("select * from view_country_group_membership"),
           na.strings = "")
datacoll <-
  sqlQuery(channel,
           paste("select * from view_TME_master_data_collection"),
           na.strings = "")

sdgdef <-
  sqlQuery(
    channel,
    paste("select * from external_indicators.view_indicator_definition"),
    na.strings = ""
  )
sdg <-
  sqlQuery(
    channel,
    paste("select * from external_indicators.view_indicator_data"),
    na.strings = ""
  )

ltbi <-
  sqlQuery(channel,
           paste("select * from view_TME_estimates_ltbi"),
           na.strings = "")

svy.agegr <-
  sqlQuery(channel,
           paste("select * from survey.age_group"),
           na.strings = "")
svy.casetype <-
  sqlQuery(channel,
           paste("select * from survey.case_type"),
           na.strings = "")
svy.screen <-
  sqlQuery(channel,
           paste("select * from survey.screen_group"),
           na.strings = "")
svy.sex <-
  sqlQuery(channel,
           paste("select * from survey.sex"),
           na.strings = "")
svy.patientgr <-
  sqlQuery(channel,
           paste("select * from survey.patient_group "),
           na.strings = "")
svy.areatype <-
  sqlQuery(channel,
           paste("select * from survey.area_type"),
           na.strings = "")
svy.cc <-
  sqlQuery(channel,
           paste("select * from survey.view_catastrophic_costs_survey"),
           na.strings = "")
svy.prevchar <-
  sqlQuery(channel,
           paste("select * from survey.view_prevalence_survey"),
           na.strings = "")
svy.prevcases <-
  sqlQuery(channel,
           paste("select * from survey.view_prevalence_survey_cases"),
           na.strings = "")
svy.prev <-
  sqlQuery(
    channel,
    paste("select * from survey.view_prevalence_survey_estimates"),
    na.strings = ""
  )
latest <-
  sqlQuery(
    channel,
    paste("select * from dcf.latest_notification"),
    na.strings = ""
  )

monthly <-
  sqlQuery(
    channel,
    paste("select * from dcf.latest_tbcases_monthly"),
    na.strings = ""
  )

covid <- 
  sqlQuery(
    channel,
    paste("select * from view_TME_master_covid_unhlm"),
    na.strings = ""
  )

odbcClose(channel)


tb2 <-
  merge(
    tb,
    pop[pop$year %in% 1980:yr, -c(1, 3, 5)],
    by = c("iso3", "year"),
    all.x = TRUE,
    all.y = TRUE
  )
tb <- copy(tb2)
rm (tb2)

dtnames <-
  c(
    'pop',
    'tb',
    'sty',
    'cty',
    'tx',
    'dic',
    'codes',
    'drs',
    'dre',
    'drnew',
    'drret',
    'drh',
    'drh2',
    'hr',
    'dr',
    'ntp',
    'dots',
    'agg',
    'mdrN',
    'prev',
    'dr.est',
    'aggdr',
    'grptypes',
    'grp',
    'grpmbr',
    'datacoll',
    'sdgdef',
    'sdg',
    'ltbi',
    'svy.agegr',
    'svy.casetype',
    'svy.screen',
    'svy.sex',
    'svy.patientgr',
    'svy.areatype',
    'svy.cc',
    'svy.prevchar',
    'svy.prevcases',
    'svy.prev',
    'latest',
    'monthly',
    'covid'
  )

dtlist <- mget(dtnames)
uname <- function(dt) {
  names(dt) <- gsub('_', '.', names(dt))
  dt <- setDT(dt)
  return(dt)
}

dtlist <- lapply(dtlist, uname)
list2env(dtlist, .GlobalEnv)


#' set keys
#' 
dtnames2 <-
  c(
    'pop',
    'tb',
    'sty',
    'cty',
    'tx',
    'drs',
    'dre',
    'drnew',
    'drret',
    'drh',
    'drh2',
    'dr',
    'hr',
    'ntp',
    'dots',
    'agg',
    'mdrN',
    'prev',
    'dr',
    'dr.est',
    'grpmbr',
    'datacoll',
    'sdg',
    'ltbi',
    'svy.cc',
    'svy.prevchar',
    'svy.prevcases',
    'svy.prev',
    'covid'
  )

dtlist2 <- mget(dtnames2)
dtlist2 <- lapply(dtlist2, function(x)
  setkey(x, iso3))
list2env(dtlist2, .GlobalEnv)
setkey(latest, iso2)

#' 
#' 
levels(dic$variable.name) <- gsub('_', '.', levels(dic$variable.name))

tb$pop <- as.double(tb$e.pop.num)
pop$pop <- as.double(pop$e.pop.num)

#' add notification rate per 100,000 pop
#' 
tb[, newinc := c.newinc * 100000 / e.pop.num]


dim(tb)
tb <-
  merge(
    tb,
    tx[, .(iso3, year, c.new.tsr, c.ret.tsr)],
    by = c('iso3', 'year'),
    all.x = TRUE,
    all.y = FALSE
  )
dim(tb)


#' add useful vars
tb[, ch := rowSums(cbind(newrel.m014, newrel.f014), na.rm = TRUE) / c.newinc]

tb[year >= 2013,
   conf := rowSums(cbind(new.labconf, ret.rel.labconf), na.rm = TRUE) / rowSums(cbind(new.labconf, ret.rel.labconf, new.clindx, ret.rel.clindx),
                                                                                na.rm = TRUE)] 
tb[year < 2013 & g.whoregion != 'EUR',
   conf := new.sp / rowSums(cbind(new.sp, new.sn, new.su), na.rm = TRUE)] 
tb[year %in% 2002:2012 & g.whoregion == 'EUR',
   conf := new.labconf / rowSums(cbind(new.sp, new.sn, new.su), na.rm = TRUE)]
tb[year %in% 2000:2001 & g.whoregion == 'EUR',
   conf := new.sp / rowSums(cbind(new.sp, new.sn, new.su), na.rm = TRUE)]

tb[, ep := rowSums(cbind(new.ep, ret.rel.ep), na.rm = TRUE) / c.newinc]
tb[, newrel.m := rowSums(cbind(newrel.m014, newrel.m15plus))]
tb[, newrel.f := rowSums(cbind(newrel.f014, newrel.f15plus))]
tb[newrel.f > 0, sex.ratio := newrel.m / newrel.f]

sdg[indicator.id=='EG.CFT.ACCS.ZS', ind := 'fuel']
sdg[indicator.id=='EN_NLD_SLUM', ind := 'slum']
sdg[indicator.id=='FINPROTECTION_CATA_TOT_10_POP', ind := 'hhexpend']
sdg[indicator.id=='GHED_CHE_pc_PPP_SHA2011', ind := 'healthexp']
sdg[indicator.id=='M_Est_smk_curr_std', ind := 'smoking']
sdg[indicator.id=='MDG_0000000029', ind := 'hiv']
sdg[indicator.id=='NCD_GLUC_04', ind := 'diabetes']
sdg[indicator.id=='NY.GDP.PCAP.PP.KD', ind := 'gdp']
sdg[indicator.id=='per_allsp.cov_pop_tot', ind := 'socprot']
sdg[indicator.id=='SA_0000001462', ind := 'alcohol']
sdg[indicator.id=='SI.POV.GINI', ind := 'gini']
sdg[indicator.id=='SI_POV_DAY1', ind := 'poverty']
sdg[indicator.id=='SN.ITK.DEFC.ZS', ind := 'undernutrition']
sdg[indicator.id=='UHC_INDEX_REPORTED', ind := 'uhc']

detach(package:RODBC)


#' # save data.tables
#'
dtlist <- mget(dtnames)

invisible(lapply(names(dtlist), function(u) {
  assign(u, dtlist[[u]])
  save(list = u, file = here(paste0("Rdata", "/", u, ".Rdata")))
  fwrite(dtlist[[u]], file = here(paste0("csv", "/", u, '_', Sys.Date(), '.csv')))
}))
