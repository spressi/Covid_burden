library(tidyverse)
library(haven)   #read spss files: better data format
library(foreign) #read spss files: better at finding correct encoding

path = "C:/Data/Corona Psychosozial/Data/" #@work
#path = path %>% gsub("C:/Data", "D:/Arbeit", ., fixed=T) #@home

se = function(x, na.rm = TRUE) {
  sd(x, na.rm) / sqrt(if(!na.rm) length(x) else sum(!is.na(x)))
}

myGgTheme = theme_bw() + theme(
  plot.title = element_text(hjust = 0.5),
  panel.background = element_rect(fill="white", color="white"),
  legend.background = element_rect(fill="white", color="grey"),
  legend.key=element_rect(fill='white'),
  axis.text = element_text(color="black"),
  axis.ticks.x = element_line(color="black"),
  axis.line.x = element_line(color="black"),
  axis.line.y = element_line(color="black"),
  legend.text = element_text(size=14, color="black"),
  legend.title = element_text(size=14, color="black"),
  strip.text.x = element_text(size=12, color="black"),
  strip.text.y = element_text(size=12, color="black"),
  axis.text.x = element_text(size=16, color="black"),
  axis.text.y = element_text(size=16, color="black"),
  axis.title = element_text(size=16, color="black"))

writeFiles = F

#the newest sample was gathered during the 2021 recession of the pandemic (T2)
#it unites subjects from different sub samples (SFB and Studi) that also have been interrogated before (see below)
data.t2.raw = "data_cnb_2021-10-19_12-06_noDescription.csv" %>% paste0(path, .) %>% read_csv2() %>% 
  filter(REF %>% grepl("test", ., ignore.case=T) == F) %>% #discard tests
  mutate(B002_01 = ifelse(is.na(B002_01), "XX", B002_01 %>% str_to_upper()), #treat NA as "XX"
         B002_02 = ifelse(is.na(B002_02), "XX", B002_02 %>% str_to_upper()),
         REF = ifelse(is.na(REF), paste(B002_01, B002_02, B002_03), REF), #use codeword as REF for subjects without REF (i.e., studi sample)
         REF = REF %>% gsub("https://www.google.de", "", .), #get rid of curious case with google in REF :D
         REF = REF %>% gsub("Z02-", "Z02_", .) %>% gsub("Z0W", "Z02W", .) %>% str_to_upper(), #manually fix errors in the REF
         REF = ifelse(REF == "XX XX NA", NA, REF)) %>% filter(REF %>% is.na() == F) #handle fully NA codewords / REFs

#data.t2.raw = data.t2.raw %>% filter(FINISHED==1) #don't filter, use as much data as possible

matching = data.t2.raw %>% select(CASE, REF, STARTED, QUESTNNR) %>% rename(datestamp2 = STARTED) #smaller data frame just for matching
  

# SFB (unique serial) ------------------------------------------------------
# this data set is comprised of two independent waves (FP2 & FP3) prior to the pandemic (T0)
# and one data set during the 2nd 2020 peak of the pandemic in fall (T1.2; not used for analysis)
path.old.sfb = "sfb_z02_cities_sharepoint-master/Online Nachbefragung/Analysis/input/" %>% paste0(path, .)
path.old.sfb.timestamps = "sfb_z02_cities_sharepoint-master/timestamps/" %>% paste0(path, .)

### T0 (pre-pandemic)
#FP2
data.old.sfb.t0.fp2 = "fp2fp3/Z02-2_WUE_MS_HH_COMPLETE.sav" %>% paste0(path.old.sfb, .) %>% 
  foreign::read.spss(to.data.frame=T) %>% tibble() %>% 
  filter(Site=="WÜ") %>% #only subjects from Würzburg (others have not been contacted for T2)
  rename(REF = Vp, sex = AFB_0, age = AFB_1) %>% #AFB_0 was really asked as sex instead of gender (gender was asked in T2)
  mutate(REF = REF %>% str_remove_all("\\s") %>% str_to_upper(),
         sex = ifelse(sex=="weiblich", "female", "male")) #convert to English (no non-binary entries)
data.old.sfb.t0.fp2 = data.old.sfb.t0.fp2 %>% 
  mutate(across(starts_with("LE") & ends_with("ges"), function(x) {x %>% as.integer() %>% ifelse(is.na(.), 0, .)}),
         across(starts_with("LE") & ends_with("ges"), function(x) {ifelse(x == 3, 1, 0)}), #only count adverse events (Scharfenort et al., 2016; https://doi.org/10.1038/tp.2016.126))
         CTQ_problem = CTQ_emotMissbr >= 13 | CTQ_körperlMissh >= 10 | CTQ_sexMissbr >= 8 | CTQ_emotVernachl >= 15 | CTQ_körperlVernachl >= 10,
         LE_sum = rowSums(across(starts_with("LE") & ends_with("ges")))) %>% #note: LE_*_ges have already been coded to be numeric by research assistants
  select(REF, sex, age, ASI_sum:PSWQ_sum, SPAI_m:LE_sum); attr(data.old.sfb.t0.fp2, "variable.labels") = NA #reduce to essential columns & erase SPSS labels

data.old.sfb.t0.fp2.wue.timestamps = "Z02_WUEEP_limesurvey_datestamps.sav" %>% paste0(path.old.sfb.timestamps, .) %>% haven::read_sav() %>% rename(REF = token) #%>% mutate(across(everything(), as.vector))
#data.old.sfb.t0.fp2.hh.timestamps = "datestamps_fp2_HH.csv" %>% paste0(path.old.sfb.timestamps, .) %>% read_delim() %>% rename(REF = VP_Code)
data.old.sfb.t0.fp2 = data.old.sfb.t0.fp2.wue.timestamps %>% #bind_rows(data.old.sfb.t0.fp2.hh.timestamps) %>% 
  right_join(data.old.sfb.t0.fp2, by="REF"); rm(data.old.sfb.t0.fp2.wue.timestamps)#; rm(data.old.sfb.t0.fp2.hh.timestamps)
#data.old.sfb.t0.fp2 %>% summarise(Missing_Timestamps = datestamp %>% is.na() %>% mean())
data.old.sfb.t0.fp2 %>% summarise(min = min(datestamp, na.rm=T), mean = mean(datestamp, na.rm=T), md = median(datestamp, na.rm=T), max = max(datestamp, na.rm=T))

#data.old.sfb.t0.fp2 %>% filter(REF %in% REF[REF %>% duplicated()]) %>% View() #check duplicates
data.old.sfb.t0.fp2 = data.old.sfb.t0.fp2 %>% filter(REF %>% duplicated() == F) #always take first instance of duplicates (seems like no difference)


#FP3
data.old.sfb.t0.fp3 = "fp2fp3/Daten_Sammlung_Z02_FP3_precorona_WU_HH_MS_sums.sav" %>% paste0(path.old.sfb, .) %>% 
  haven::read_sav() %>% rename(REF = VP, Site = Standort, sex = AFB_0, age = AFB_1) %>% 
  filter(sex %>% is.na() == F) %>% #kick out subjects that have empty data
  filter(Site==1) %>% #only subjects from Würzburg (others have not been contacted for T2)
  mutate(CTQ_problem = CTQ_emotMissbr >= 13 | CTQ_körperlMissh >= 10 | CTQ_sexMissbr >= 8 | CTQ_emotVernachl >= 15 | CTQ_körperlVernachl >= 10)

# Problems with LE questionnaire (all solved now)
# data.old.sfb.t0.fp3 %>% mutate(LE_sum = ifelse(LE_sum %>% is.na(), 0, LE_sum), #note: LE_sum is the sum of LE_*_gesN (not LE_*_ges as in FP2). LE_*_gesN is a numeric recode of LE_*_ges with "+" -> 1, "x" -> 2, "-" -> 3
#                                across(starts_with("LE") & ends_with("gesN"), function(x) {x %>% as.integer() %>% ifelse(is.na(.), 0, .)}),
#                                LE_sum2 = rowSums(across(starts_with("LE") & ends_with("gesN")))) %>% 
#   summarise(LE_check = mean(LE_sum == LE_sum2)) #check works
#data.old.sfb.t0.fp3 %>% select(starts_with("LE") & contains("ges")) %>% .[, order(names(.))] %>% View("ges vs. gesN") #classification is suboptimal => redo
data.old.sfb.t0.fp3 %>% pivot_longer(cols=starts_with("LE") & ends_with("ges")) %>% pull(value) %>% unique() %>% sort() #TODO interpret this :D
#data.old.sfb.t0.fp3 %>% pivot_longer(cols=starts_with("LE") & ends_with("ges")) %>% pull(value) %>% table() #get counts

data.old.sfb.t0.fp3 = data.old.sfb.t0.fp3 %>% 
  mutate(REF = REF %>% gsub("Z0W", "Z02W", .) %>% str_to_upper(),
         sex = ifelse(sex==0, "female", "male"),
         across(starts_with("LE") & ends_with("ges"), 
                function(x) {
                  x %>% tolower() %>%
                    #treat special cases in the data first
                    {case_when(. %in% c("", "0", "0x") ~ 0, #exact matches (to exclude e.g. "10x")
                              grepl("+", ., fixed=T) & grepl("-", ., fixed=T) ~ 2, #set ambivalence (+ AND -) equal to neutral (x)
                              grepl("*", ., fixed=T) ~ 1, #assume * is + with capslock on (German keyboard)
                              
                              grepl("+", ., fixed=T) ~ 1,
                              grepl("x", ., fixed=T) ~ 2,
                              grepl("-", ., fixed=T) ~ 3,
                              T ~ 2)} #treat anything else as a 2
                }, .names="{.col}N"),
         across(starts_with("LE") & ends_with("ges"), function(x) {ifelse(x == 3, 1, 0)}), #only count adverse events (Scharfenort et al., 2016; https://doi.org/10.1038/tp.2016.126))
         LE_sum = rowSums(across(starts_with("LE") & ends_with("gesN"))),
         GSE_sum = rowSums(across(starts_with("GSE")))) %>% #select(starts_with("LE") & contains("ges")) %>% .[, order(names(.))] %>% View("ges vs. gesN corrected")
  select(REF, sex, age, ASI_sum:BDI_sum, GSE_sum, CTQ_problem) #reduce to essential columns
data.old.sfb.t0.fp3 = "Datestamp_Z02_FP3_WÜ_.sav" %>% paste0(path.old.sfb.timestamps, .) %>% 
  haven::read_sav() %>% rename(REF = VP_Code) %>% na.omit() %>% 
  right_join(data.old.sfb.t0.fp3, by="REF")
#data.old.sfb.t0.fp3 %>% summarise(Missing_Timestamps = datestamp %>% is.na() %>% mean())
#data.old.sfb.t0.fp3 %>% filter(datestamp %>% is.na())
data.old.sfb.t0.fp3 %>% summarise(min = min(datestamp, na.rm=T), mean = mean(datestamp, na.rm=T), md = median(datestamp, na.rm=T), max = max(datestamp, na.rm=T))

#data.old.sfb.t0.fp3 %>% filter(REF %in% REF[REF %>% duplicated()]) #check duplicates (none)


#T0: FP2 + FP3
data.sfb.t0 = data.old.sfb.t0.fp2 %>% bind_rows(data.old.sfb.t0.fp3) %>% 
  rename(asi3_phys = ASI3_PH, stai_trait = STAI_sum) %>% 
  rename_with(function(x) {x %>% tolower() %>% gsub("_sum", "", .) %>% gsub("_mean", "", .) %>% gsub("_m$", "", .) %>% gsub("soz", "soc", .)}, -REF) %>% 
  mutate(strain = bind_cols(scale(asi3), scale(stai_trait), scale(pswq)) %>% rowMeans(),
         cerq_maladapt = cerq_selfblame + cerq_rumination + cerq_catastrophizing + cerq_otherblame,
         cerq_adapt = cerq_positiverefocusing + cerq_refocusonplanning + cerq_positivereappraisal + cerq_puttingintoperspective) %>% 
  select(REF:stai_trait, strain, everything())
data.sfb.t0 %>% mutate(ctqcheck = data.sfb.t0 %>% select(contains("ctq_")) %>% rowSums() %>% {. - ctq_bagatellisierung}) %>% mutate(diff = ctq-ctqcheck) %>% pull(diff) %>% table() #ctq is not the sum of subscores?
rm("data.old.sfb.t0.fp2", "data.old.sfb.t0.fp3")
#data.sfb.t0 %>% write_rds("data.sfb.t0.rds" %>% paste0(path, .))


### T1.2 (upturn in fall 2020, not used)
data.sfb.t1 = "survey_z02_nachbefragung/data_sfbz2_2020-10-30_12-24.csv" %>% paste0(path.old.sfb, .) %>% 
  read_tsv(locale = locale(encoding = "UTF-16LE")) %>% 
  filter(REF != "") %>%   # Remove cases with empty participant id 
  filter(!REF %in% c(
    "3Z02WU168","3Z02WU267","3Z02WU307","3Z02WU494","3Z02WU496", "3Z02WU503", "3Z02WU504", "3Z02WU505","3Z02WU512", "3Z02WU513", "3Z02WU462", "3Z02WU461", "3Z02WU063", "3Z02WU102", "3Z02WU223", "3Z02WU500", "3Z02WU495"   # These participants where accidentally contacted but invalid, since they completed to peivous part of the longitudinal study (FP2, FP3) already during the corona pandemic.
  )) %>%
  filter(!CASE %in% c(
    224, 1854, 1550, 1551, 1552, 1554, 1555, 1557, 2422, 2423, 2424, 2425, 1845, 2110, 1054, 543, 633, 636, 979, 745   # Remove these cases with duplicate participant id. Participants where asked after the survey which is their correct bank number. We compared this number with the bank number entered into the survey to identify the correct case for the participant.
  )) %>%
  filter(!REF %in% c(
    "*3Z02HH999*"  # Test participant HH
  )) %>%
  filter(LASTPAGE >= 55) %>%   # Remove cases that did not reach last page
  distinct(REF, .keep_all = TRUE) %>%   # Keep only first case if several cases with same participant id exist
  rename("Site" = "B053_01") %>%
  filter(Site %in% c("Ham", "Mue", "Wue")) %>%   # remove cases with invalid location
  filter(B029 == 1) %>%   # Only participants with correct validation answer
  mutate(REF = str_replace(REF, "^(\\d)$", "Z02_HH_00\\1")) %>%   # Adjust Hamburg participant codes to fp2fp3 nomenclature 
  mutate(REF = str_replace(REF, "^(\\d\\d)$", "Z02_HH_0\\1")) %>%
  mutate(REF = str_replace(REF, "^(\\d\\d\\d)$", "Z02_HH_\\1")) %>% 
  
  filter(REF %>% grepl("WU", .)) #only subjects from Würzburg (others have not been contacted for T2)
  #TODO reduce to essential columns (canceled because data set not used)

#data.sfb.t1 %>% summarise(min = min(STARTED, na.rm=T), mean = mean(STARTED, na.rm=T), md = median(STARTED, na.rm=T), max = max(STARTED, na.rm=T))

#data.sfb.t1 %>% filter(REF %in% REF[REF %>% duplicated()]) #check duplicates (none)
#data.sfb.t1 %>% write_rds("data.sfb.t1.rds" %>% paste0(path, .))


### Matching to T2
matching.t0 = data.sfb.t0 %>% select(REF:datestamp) %>% rename(datestamp0 = datestamp)
matching.t1 = data.sfb.t1 %>% select(REF, STARTED) %>% rename(datestamp1 = STARTED)
#matching.t0 %>% full_join(matching.t1, by="REF") %>% filter(datestamp0 %>% is.na()) #one person in T2 who was not in T0?
matching.sfb = matching.t0 %>% filter(datestamp0 %>% is.na() == F) %>% left_join(matching.t1, by="REF"); rm(matching.t0); rm(matching.t1)

#matching
matching.sfb = matching.sfb %>% right_join(matching %>% filter(QUESTNNR=="sfb") %>% select(REF, datestamp2), by="REF")
#matching.sfb %>% filter(datestamp0 %>% is.na()) #one serial not found?
#data.sfb.t0 %>% filter(REF=="3Z02WU461") #no entry but invited to T2?
matching.sfb = matching.sfb %>% filter(datestamp0 %>% is.na() == F)
matching.sfb %>% summarise(n0 = sum(datestamp0 %>% is.na() == F), n1 = sum(datestamp1 %>% is.na() == F), n2 = sum(datestamp2 %>% is.na() == F))
#TODO may contain massive NAs in questionnaires

matching.sfb %>% pivot_longer(datestamp0:datestamp2, names_to="time", values_to="datestamp") %>% mutate(time = time %>% gsub("datestamp", "", .)) %>% 
  group_by(time) %>% summarise(min = min(datestamp, na.rm=T), mean = mean(datestamp, na.rm=T), md = median(datestamp, na.rm=T), max = max(datestamp, na.rm=T))


serials.matched = matching.sfb %>% filter(datestamp0 %>% is.na() == F) %>% pull(REF)
#serials.matched %>% duplicated() %>% sum() #check duplicate matches (none)
#data.t2.raw %>% filter(REF %in% serials.matched) %>% View()


# Selective Attrition -----------------------------------------------------
data.sfb.t0 %>% mutate(stayed = REF %in% serials.matched) %>% 
  group_by(stayed) %>% summarise(n = n(),
                                 across(c("asi3", "pswq", "stai_trait", "spai", "lsas", "gse", "cerq_maladapt", "cerq_adapt", "cerq_acceptance", "ctq", "lte", "le"),
                                        .fn=function(x) mean(x, na.rm=T)))
                                        #.fn=mean, na.rm=T))
                                        #.fn=list(m = "mean", sd = "sd"), na.rm=T))
df = data.sfb.t0 %>% mutate(stayed = REF %in% serials.matched)
attrition = tibble(var = c("asi3", "pswq", "stai_trait", "spai", "lsas", "gse", "cerq_maladapt", "cerq_adapt", "cerq_acceptance", "ctq", "lte", "le")) %>% 
  rowwise() %>% 
  mutate(test = t.test(df %>% filter(stayed==T) %>% pull(var),
                       df %>% filter(stayed==F) %>% pull(var), 
                       var.equal=T) %>% list(),
         p = test$p.value,
         sig = case_when(p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", T ~ ""),
         apa = test %>% apa::t_apa(print=F, es_ci=T)) %>% ungroup()
#attrition
attrition %>% filter(p == min(p))


# Studi / "Base" ----------------------------------------------------------
# this data set is comprised of one sample during the first peak of the pandemic in April 2020 (T1)

#find duplicates of code in the new data set
codes = matching %>% filter(QUESTNNR=="base") %>% pull(REF) #only for base subsample ("studi")
duplicates = codes[codes %>% duplicated()] %>% unique()
#data.t2.raw %>% filter(REF %in% duplicates) %>% arrange(REF) %>% View()

data.old.studi.raw = "Studierendenbefragung/data_corona-2020_2020-05-01_11-59_noDescription.csv" %>% paste0(path, .) %>% read_csv2() %>%
  #filter((B002_01 %>% is.na() | B002_02 %>% is.na() | B002_03 %>% is.na()) == F) %>% 
  filter(B002_03 %>% is.na() == F) %>% 
  mutate(B002_01 = ifelse(is.na(B002_01), "XX", B002_01 %>% str_to_upper()),
         B002_02 = ifelse(is.na(B002_02), "XX", B002_02 %>% str_to_upper())) %>% 
  separate(B002_03, into=c("DD", "MM", "YYYY"), remove=F) %>% 
  mutate(SERIAL = paste(B002_01, B002_02, paste(YYYY, MM, DD, sep="-")) %>% str_to_upper(),
         STARTED = STARTED %>% as.Date(format="%d.%m.%Y"), LASTDATA = LASTDATA %>% as.Date(format="%d.%m.%Y")) %>% 
  rename(datestamp = STARTED)
#data.old.studi.raw %>% pull(datestamp) %>% summary()
data.old.studi.raw %>% summarise(min = min(datestamp, na.rm=T), mean = mean(datestamp, na.rm=T), md = median(datestamp, na.rm=T), max = max(datestamp, na.rm=T))
data.old.studi.raw %>% pull(YYYY) %>% as.integer() %>% summary() #incorrect entries

codes.old.studi = data.old.studi.raw %>% pull(SERIAL) %>% sort()
duplicates.old.studi = codes.old.studi[codes.old.studi %>% duplicated()] %>% unique()

codes.unique = codes %>% setdiff(duplicates) %>% setdiff(duplicates.old.studi) %>% sort()
#{codes %>% setdiff(duplicates) %>% length()} - {codes %>% setdiff(duplicates) %>% setdiff(duplicates.old.studi) %>% length()}
codes.matched = codes.unique[codes.unique %in% codes.old.studi]
codes.unmatched = codes.unique[codes.unique %in% codes.old.studi == F]

#data.t2.raw %>% filter(QUESTNNR=="base") %>% filter(REF %in% codes.matched == F) %>% summarise(N = n()) #failed matches before manual matching

# Aid Manual Matching
matching.studi = matching %>% filter(QUESTNNR=="base") %>% separate(REF, into=c("B002_01", "B002_02", "YYYY", "MM", "DD"), remove=F)
matching.manual = tibble(concordance = 0, new = codes.unmatched, old = "", comment = "")
for (c in codes.unmatched) { #find best matches
  matching.studi.vp = matching.studi %>% filter(REF==c)
  match.vp = data.old.studi.raw %>% select(CASE:YYYY) %>%
    mutate(B002_01.match = B002_01==matching.studi.vp$B002_01, B002_02.match = B002_02==matching.studi.vp$B002_02,
           DD.match = DD==matching.studi.vp$DD, MM.match = MM==matching.studi.vp$MM, YYYY.match = YYYY==matching.studi.vp$YYYY,
           match = DD.match + MM.match + YYYY.match + B002_01.match + B002_02.match) %>%
    filter(DD.match | MM.match | YYYY.match | B002_01.match | B002_02.match) %>% arrange(desc(match))
  
  match.vp.max = match.vp %>% pull(match) %>% max()
  matching.manual[matching.manual$new==c, "concordance"] = match.vp.max

  #all codes with maximum matching (match.vp.max)
  match.vp.potential = match.vp %>% filter(match==match.vp.max)
  if (nrow(match.vp.potential) == 1) { #only one match => assume correct
    matching.manual[matching.manual$new==c, "old"] = match.vp.potential$SERIAL
    matching.manual[matching.manual$new==c, "comment"] = "auto match"
  }
  
  if (writeFiles) match.vp.potential %>% write_excel_csv2(file=paste0(path, "matching/", match.vp.max, "_", c, ".csv"))
}
if (writeFiles) matching.manual %>% write_excel_csv2(file=paste0(path, "matching/matching.csv"))


# Apply Manual Matching
matching.manual.all = "matching/matching.xlsx" %>% paste0(path, .) %>% readxl::read_xlsx(sheet="matching")
#matching.manual.all %>% filter(comment!="auto match", decision %in% c("+", "?")) %>% arrange(desc(concordance), decision) %>% View("examples of manual matching")

matching.manual = matching.manual.all %>% filter(decision %in% c("+", "?", "??", "???"))
#matching.manual %>% filter(comment %>% grepl("Marina", .) | comment %>% grepl("Iris", .) | comment %>% grepl("Tippfehler", .) | (concordance == 4 & comment %>% grepl("auto match", .))) %>% arrange(concordance) %>% View("example manual matches")
matching.duplicates = matching.manual %>% filter(comment %>% grepl("doppelt", .)) %>% rename(valid = comment2) %>% 
  separate(comment, into=c("keyword", "num1", "num2")) %>% mutate(invalid = ifelse(valid==num1, num2, num1)) %>% 
  select(new, old, valid, invalid)

data.old.studi.raw = data.old.studi.raw %>% filter(CASE %in% matching.duplicates$invalid == F) #kick old duplicates that have been matched manually

data.t2.raw$REF[data.t2.raw$REF %in% matching.manual$new] = matching.manual$old #switch new unmatched codes to their old manual match

#redo matching procedure (but with manipulated data)
codes = data.t2.raw %>% filter(QUESTNNR=="base") %>% pull(REF)
duplicates = codes[codes %>% duplicated()] %>% unique()
#data.t2.raw %>% filter(REF %in% duplicates) %>% arrange(REF) %>% View()

codes.unique = codes %>% setdiff(duplicates) %>% setdiff(duplicates.old.studi) %>% sort()
#{codes %>% setdiff(duplicates) %>% length()} - {codes %>% setdiff(duplicates) %>% setdiff(duplicates.old.studi) %>% length()}
codes.matched = codes.unique[codes.unique %in% codes.old.studi]
codes.unmatched = codes.unique[codes.unique %in% codes.old.studi == F]

#data.t2.raw %>% filter(QUESTNNR=="base") %>% filter(REF %in% codes.matched == F) %>% summarise(N = n()) #failed matches after manual matching


# Usable data set ---------------------------------------------------------
data.t2.raw.matched = data.t2.raw %>% filter(REF %in% codes.matched | REF %in% serials.matched)
data.t2.raw.matched %>% group_by(QUESTNNR) %>% summarise(N = n())

#data.t2.raw.matched %>% write_rds("data_t2_raw_matched.rds" %>% paste0(path, .)); data.t2.raw.matched %>% write_excel_csv2("data_t2_raw_matched.csv" %>% paste0(path, .), na="")

#age conversion
age.sfb = data.sfb.t0 %>% select(REF, datestamp, age) %>% rename(datestamp0 = datestamp, age0 = age) %>% 
  right_join(data.t2.raw.matched %>% filter(QUESTNNR=="sfb") %>% select(REF, STARTED) %>% rename(datestamp2 = STARTED)) %>% 
  mutate(age2 = age0 + {difftime(datestamp2, datestamp0, units="days") / 365.25} %>% round() %>% as.integer())
age.sfb %>% pull(age0) %>% summary()
age.sfb %>% pull(age2) %>% summary()
