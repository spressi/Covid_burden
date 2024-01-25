library(tidyverse)
library(psych); library(ltm) #Cronbach's alpha
library(nFactors) #factor analysis
#library(Matching); library(rgenoud)

missing.max = .2
writeFiles = F

missingItems.vars = function(vars) return(rowSums(is.na(across(vars))) / ncol(across(vars)))
sumscore.vars = function(vars, missing.max=1) {
  ifelse(missingItems.vars(vars) > missing.max, NA,
         rowMeans(across(vars), na.rm=T) * ncol(across(vars)))
}

missingItems = function(identifier) return(rowSums(is.na(across(starts_with(identifier)))) / ncol(across(starts_with(identifier))))
sumscore = function(identifier, missing.max=1) {
  ifelse(missingItems(identifier) > missing.max, NA,
         rowMeans(across(starts_with(identifier)), na.rm=T) * ncol(across(starts_with(identifier))))
}



# Preprocessing -----------------------------------------------------------
#data.t2.raw.matched = read_rds("data_t2_raw_matched.rds" %>% paste0(path, .))

data.t2 = data.t2.raw %>% 
  mutate( ## sum scores of standardized questionnaires
    #order of columns is different from temporal order of questionnaires! => re-establish temporal order here
    
    #CD RISC 10 (CD-RISC)
    across(starts_with("B039"), ~ . - 1), #lowest = 0
    cdr = sumscore("B039", missing.max),
    cdr_flex = sumscore.vars(c("B039_01", "B039_05"), missing.max),
    cdr_effic = sumscore.vars(c("B039_02", "B039_04", "B039_09"), missing.max),
    cdr_regulate = B039_10,
    cdr_optim = sumscore.vars(c("B039_03", "B039_06", "B039_08"), missing.max),
    cdr_focus = B039_07,
    
    #Body Questionnaire: Somatic Symptom Scale (Body)
    across(starts_with("B040"), ~ . - 1), #lowest = 0
    body = sumscore("B040", missing.max),
    
    #Fear of Covid-19 Scale (FoC-19)
    foc = sumscore("FC01", missing.max),
    
    #Patient Health Questionnaire (PHQ-2) [Depression Screening]
    across(starts_with("PH01"), ~ . - 1), #lowest = 0
    phq2 = sumscore("PH01", missing.max),
    
    #State-Trait Anxiety Inventory (STAI): State
    across(c(ST03_01, ST03_15, ST03_16), ~ 5 - .), #STAI state reverse items
    stai_state = sumscore("ST03", missing.max),
    
    #Intolerance of Uncertainty (IUS)
    ius_p = sumscore.vars(c("IU02_07", "IU02_08", "IU02_10", "IU02_11", "IU02_18", "IU02_19", "IU02_21"), missing.max),
    ius_i = sumscore.vars(c("IU02_09", "IU02_12", "IU02_15", "IU02_20", "IU02_25"), missing.max),
    
    #Penn State Worry Questionnaire (PSWQ)
    across(c(PS01_01, PS01_03, PS01_08, PS01_10, PS01_11), ~ 6 - .), #reverse items
    pswq = sumscore("PS01", missing.max),
    
    #Anxiety Sensitivity Index (ASI-3)
    across(starts_with("AS01"), ~ . - 1), #lowest value 0 instead of 1
    asi3_phys = sumscore.vars(c("AS01_03", "AS01_04", "AS01_07", "AS01_08", "AS01_12", "AS01_15"), missing.max),
    asi3_cog = sumscore.vars(c("AS01_02", "AS01_05", "AS01_10", "AS01_14", "AS01_16", "AS01_18"), missing.max),
    asi3_soc = sumscore.vars(c("AS01_01", "AS01_06", "AS01_09", "AS01_11", "AS01_13", "AS01_17"), missing.max),
    asi3 = asi3_phys + asi3_cog + asi3_soc,
    
    #Illness Attitude Scales (IAS)
    across(starts_with("IA01"), ~ . - 1), #lowest = 0
    ias = sumscore("IA01", missing.max),
    
    #State-Trait Anxiety Inventory (STAI): Trait
    across(c(ST02_01, ST02_03, ST02_06, ST02_07, ST02_10, ST02_13, ST02_14, ST02_16, ST02_19), ~ 5 - .), #STAI trait reverse items
    stai_trait = sumscore("ST02", missing.max),
    
    #Allgemeine Depressionsskala Kurzform [General Depression Scale, short] (ADS-K)
    across(starts_with("AD01"), ~ . - 1), #lowest = 0
    across(c(AD01_09, AD01_12), ~ 3 - .), #reverse items
    adsk = sumscore("AD01", missing.max),
    
    #Loneliness and Isolation during Social Distancing (LISD)
    across(c(SD01_09, SD04_03, SD04_07, SD04_08, SD04_10), ~ 6 - .), #reverse items
    lisd_state_lonely  = sumscore.vars(c("SD01_01", "SD01_03", "SD01_04", "SD01_05", "SD01_06", "SD01_07", "SD01_10", "SD01_11", "SD01_12"), missing.max)/9, #LISD uses mean item scores instead of sums
    lisd_state_support = sumscore.vars(c("SD01_02", "SD01_08", "SD01_09"), missing.max)/3, #LISD uses mean item scores instead of sums
    lisd_trait_lonely  = sumscore.vars(c("SD04_02", "SD04_09", "SD04_12", "SD04_13"), missing.max)/4, #LISD uses mean item scores instead of sums
    lisd_trait_belong  = sumscore.vars(c("SD04_01", "SD04_04", "SD04_05", "SD04_10", "SD04_11"), missing.max)/5, #LISD uses mean item scores instead of sums
    lisd_trait_close   = sumscore.vars(c("SD04_03", "SD04_06", "SD04_07", "SD04_08"), missing.max)/4, #LISD uses mean item scores instead of sums
    
  ) %>% 
  
  mutate( ##summarise non-standardized questions
    #BASIC SOCIODEMOPGRAPHICS
    #sociodemographics
    gender = case_when(B013 == 1 ~ "male",
                       B013 == 2 ~ "female",
                       B013 == 3 ~ "diverse",
                       T ~ "NA") %>% ifelse(. == "NA", NA, .) %>% as_factor(),
    
    age = {difftime(STARTED, B002_03, units="days") / 365.25} %>% round() %>% as.integer(), #only works for QUESTNNR=="base" (inject SFB age later)
    
    vaccination = case_when(B037 == 1 ~ "full",
                            B037 == 2 ~ "once",
                            B037 == 3 ~ "planned",
                            B037 == 4 ~ "no interest",
                            B037 == 5 ~ "not possible",
                            T ~ "NA") %>% ifelse(. == "NA", NA, .) %>% as_factor(), 
    vaccination_full = vaccination=="full",
    risk_group = B011_14, 
    city = case_when(B034 == 1 ~ "rural",
                     B034 == 2 ~ "city <= 100K",
                     B034 == 3 ~ "city > 100K",
                     T ~ "NA") %>% ifelse(. == "NA", NA, .) %>% as_factor(), 
    validation = B024 == 1,
    
    citizenship = case_when(B004 == 36 ~ "german", B004 != 36 ~ "other", T ~ "NA") %>% ifelse(. == "NA", NA, .) %>% as_factor(),
    
    #living situation
    living_with_partner = factor(B017_01, levels = c ("2", "1"), labels = c ("yes", "no")),
    living_with_pet = factor(B017_06, levels = c ("2", "1"), labels = c ("yes", "no")),
    living_with_child = case_when(B017_02 == 2|B017_07== 2 ~ "yes",
                                  B017_02 == 1|B017_07== 1 ~ "no",
                                  T ~ "NA") %>% ifelse(. == "NA", NA, .) %>% as_factor(),
    living_with_other = case_when(B017_08 == 2|B017_03== 2|B017_04 == 2|B017_05== 2 ~ "yes",
                                  B017_08 == 1|B017_03== 1|B017_04 == 1|B017_05== 1 ~ "no",
                                  T ~ "NA") %>% ifelse(. == "NA", NA, .) %>% as_factor(),
    living_alone = living_with_partner == "no" & living_with_pet == "no" & living_with_child == "no" & living_with_other == "no",
    
    #occupation
    occupation = case_when(B019_06 ==2 |B019_09 ==2 ~ "studies/apprenticeship",
                           B019_07 ==2 ~ "unemployed",
                           B019_06 ==1 & B019_09 ==1 & B019_01 ==2 ~ "employed",    #TODO "sonstiges" is missing here!!
                           T ~ "NA") %>% ifelse(. == "NA", NA, .) %>% as_factor(),
    unemployed = occupation=="unemployed"
    
  ) %>% 
  
  #inject SFB age
  left_join(age.sfb %>% select(REF, age2), by="REF") %>% mutate(age = ifelse(QUESTNNR=="base", age, age2)) %>% select(-age2) %>% 
  
  #(RECENT) DAILY STRUCTURE
  #positive_dayroutine
  
  mutate_at(   # polarity reversal of some daily structure items 
    paste0("B032_", sprintf("%02d", c(3,6,8))),
    ~ recode(., `5` = 1, `4` = 2, `3` = 3, `2` = 4, `1` = 5)
  ) %>%
  mutate(daystructure = rowMeans(.[paste0("B032_", sprintf("%02d", 1:8))]))%>%   # higher value indicates positive day routine
  
  #COVID EXPOSURE
  #COVID infection others
  rowwise()%>% mutate(infection_private = ifelse(B025_01_4 ==2, 4,ifelse(B025_01_3 ==2, 3,ifelse(B025_01_2 ==2, 2, ifelse(B025_01_1 ==2, 1, 0)))))%>%
  rowwise()%>% mutate(infection_job = ifelse(B025_02_4 ==2, 4,ifelse(B025_02_3 ==2, 3,ifelse(B025_02_2 ==2, 2, ifelse(B025_02_1 ==2, 1, 0)))))%>%
  rowwise()%>% mutate(infection_other = ifelse(B025_03_4 ==2, 4,ifelse(B025_03_3 ==2, 3,ifelse(B025_03_2 ==2, 2, ifelse(B025_03_1 ==2, 1, 0)))))%>%
  rowwise()%>% mutate(infection_acquaintances = max(c(infection_private,infection_job,infection_other)))%>% # summary three variables above
  #COVID infection self
  rowwise()%>% mutate(infection_self = ifelse(B025_04_3 ==2, 3,ifelse(B025_04_2 ==2, 2,ifelse(B025_04_1 ==2, 1, 0))))%>%
  ungroup()%>%
  
  #COVID risk
  mutate(
    covid_risk = rowMeans(.[paste0("B011_", sprintf("%02d",c(1,2,14,4,5)))]),
    #COVID AFFECTEDNESS
    covid_behavior_change = rowMeans(.[paste0("B011_", sprintf("%02d",6:8))]),
    covid_opinion_change = rowMeans(.[paste0("B011_", sprintf("%02d",12:13))]), # better use items separately?
    
    covid_efficiency_change = rowMeans(.[paste0("B028_", sprintf("%02d",1:2))]), # zusätzlich Kategorien zusammenfassen:1-2,3, 4-5?
    covid_psychological_change = 6 - rowMeans(.[paste0("B028_", sprintf("%02d",c(3,8)))]), #revert scale to match greater ~ more psychological strain
    covid_psychological_change_emo = 6 - B028_03,
    covid_psychological_change_sleep = 6 - B028_08,
    covid_situation_change = rowMeans(.[paste0("B028_", sprintf("%02d",4:7))]),# better use items separately?
    
    covid_media_update = case_when(B009 ==1 ~ "never",
                                   B009 ==2 ~ "every few days",
                                   B009 ==3 ~ "once a day", 
                                   B009 >3 ~ "several times a day", 
                                   T ~ "NA") %>% ifelse(. == "NA", NA, .) %>% as_factor(),
    
    covid_worry_livelihood = rowMeans(.[paste0("B011_", sprintf("%02d",9:10))]),
    #safety behaviour
    #summarize categories
    contacts_real_01 = case_when(B012_01  == 1 ~ "1",#"no"
                            B012_01  >1 & B012_01 <5 ~ "2", #"partly"
                            B012_01  == 5 ~ "3", #"absolutely"
                            T ~ "NA") %>% ifelse(. == "NA", NA, .) %>% as.integer(),
    contacts_real_02 = case_when(B012_02  == 1 ~ "1",
                            B012_02  >1 & B012_02 <5 ~ "2", 
                            B012_02  == 5 ~ "3", 
                            T ~ "NA") %>% ifelse(. == "NA", NA, .) %>% as.integer(),
    contacts_real_03 = case_when(B012_03  == 1 ~ "1",
                            B012_03  >1 & B012_03 <5 ~ "2", 
                            B012_03  == 5 ~ "3", 
                            T ~ "NA") %>% ifelse(. == "NA", NA, .) %>% as.integer(),
    contacts_real_04 = case_when(B012_04  == 1 ~ "1",
                            B012_04 >1 & B012_04 <5 ~ "2", 
                            B012_04  == 5 ~ "3", 
                            T ~ "NA") %>% ifelse(. == "NA", NA, .) %>% as.integer(),
    contacts_real_05 = case_when(B012_05  == 1 ~ "1",
                            B012_05  >1 & B012_05 <5 ~ "2", 
                            B012_05  == 5 ~ "3", 
                            T ~ "NA") %>% ifelse(. == "NA", NA, .) %>% as.integer(),
    contacts_real_06 = case_when(B012_06  == 1 ~ "1",
                            B012_06  >1 & B012_06 <5 ~ "2", 
                            B012_06  == 5 ~ "3", 
                            T ~ "NA") %>% ifelse(. == "NA", NA, .) %>% as.integer()
  )%>%
  mutate(
    
    # contacts_avoidance = rowMeans(.[paste0("B012", sprintf("%02d",c(1,2,5)))]), # original response scale
    #  contacts_safety = rowMeans(.[paste0("B012", sprintf("%02d",c(3,4,6)))]), # original response scale
    contacts_real_avoidance = rowMeans(.[paste0("contacts_real_", sprintf("%02d",c(1,2,5)))]),
    contacts_real_safety = rowMeans(.[paste0("contacts_real_", sprintf("%02d",c(3,4,6)))]),
    safety_behaviour = rowMeans(.[paste0("contacts_real_", sprintf("%02d",1:6))])# summary of two variables above
  )%>%
  
  #SOCIAL FACTORS
  #need for social contacts
  mutate_at(   # polarity reversal of some contacts_need items 
    paste0("B021_", sprintf("%02d", c(5,7))),
    ~ recode(., `5` = 1, `4` = 2, `3` = 3, `2` = 4, `1` = 5)
  ) %>%
  mutate( #summarize categories
    contacts_need_01 = case_when(B021_01 <4 ~ "1", # "1"=low
                                 B021_01 >3 ~ "2", # "2"=high
                                 T ~ "NA") %>% ifelse(. == "NA", NA, .) %>% as.integer(),
    contacts_need_02 = case_when(B021_02 <4 ~ "1",
                                 B021_02 >3 ~ "2", 
                                 T ~ "NA") %>% ifelse(. == "NA", NA, .) %>% as.integer(),
    contacts_need_03 = case_when(B021_05 <4 ~ "1",
                                 B021_05 >3 ~ "2", 
                                 T ~ "NA") %>% ifelse(. == "NA", NA, .) %>% as.integer(),
    contacts_need_04 = case_when(B021_06 <4 ~ "1",
                                 B021_06 >3 ~ "2", 
                                 T ~ "NA") %>% ifelse(. == "NA", NA, .) %>% as.integer(),
    contacts_need_05 = case_when(B021_07 <4 ~ "1",
                                 B021_07 >3 ~ "2", 
                                 T ~ "NA") %>% ifelse(. == "NA", NA, .) %>% as.integer()
  )%>%
  mutate(
    contacts_need = rowMeans(.[paste0("contacts_need_", sprintf("%02d",1:5))]), 
    
    #outgroup opinion
    #summarize categories
    opinion_outgroup_01 = case_when(B021_03 <3 ~ "2",# positive
                                    B021_03 >2 ~ "1", # negative
                                    T ~ "NA") %>% ifelse(. == "NA", NA, .) %>% as.integer(),
    opinion_outgroup_02 = case_when(B021_04 <3 ~ "2",
                                    B021_04 >2 ~ "1", 
                                    T ~ "NA") %>% ifelse(. == "NA", NA, .) %>% as.integer()
  )%>%
  mutate(
    opinion_outgroup = rowMeans(.[paste0("opinion_outgroup_", sprintf("%02d",1:2))])# higher values indicate more positive opinion
    
  )%>% 
  
  filter(cdr %>% is.na() == F) %>% 
  #mutate(strain = bind_cols(scale(asi3), scale(stai_trait), scale(pswq)) %>% rowMeans()) %>% #do this later for all data to not lose information on baseline differences
  select(CASE:QUESTNNR, gender, age, vaccination:validation,
         #strain, 
         asi3, stai_trait, pswq, #criteria
         cdr:lisd_trait_close, #standardized questionnaires (predictors)
         citizenship:opinion_outgroup, #unstandardized questions (predictors)
         
         B002_01:B002_03, #individual parts of REF (if for QUESTNNR=="base")
         STARTED, LASTDATA, FINISHED, LASTPAGE:DEG_TIME, TIME_SUM, TIME001:TIME022)

data.t2 %>% pull(validation) %>% summary() #TODO include NAs? (i.e., subjects who did not complete the questionnaire)?
#data.t2 = data.t2 %>% filter(validation | validation %>% is.na()) %>% select(-validation) %>% 
data.t2 = data.t2 %>% filter(validation) %>% select(-validation) %>% 
  mutate(strain = bind_cols(scale(asi3), scale(stai_trait), scale(pswq)) %>% rowMeans()) %>% 
  select(CASE:pswq, strain, everything())
#data.t2 %>% write_rds("data_t2.rds" %>% paste0(path, .))

#glance at missing values (doesn't work anymore to make the data more compact; need to call "missingItems" explicitly in mutate of data.t2.raw)
#data.t2 %>% select(CASE, contains(".missing")) %>% pivot_longer(-CASE, names_to="Questionnaire", values_to="missing") %>% filter(missing != 0, missing != 1)
#data.t2 %>% select(contains(".missing")) %>% summarise(., missing.mean = rowMeans(.)) %>% filter(missing.mean != 0, missing.mean != 1)


data.old.studi = data.old.studi.raw %>% 
  filter(SERIAL %in% codes.matched) %>% 
  mutate(
    #Anxiety Sensitivity Index (ASI-3)
    across(starts_with("AS01"), ~ . - 1), #lowest value 0 instead of 1
    asi3_phys = sumscore.vars(c("AS01_03", "AS01_04", "AS01_07", "AS01_08", "AS01_12", "AS01_15"), missing.max),
    asi3_cog = sumscore.vars(c("AS01_02", "AS01_05", "AS01_10", "AS01_14", "AS01_16", "AS01_18"), missing.max),
    asi3_soc = sumscore.vars(c("AS01_01", "AS01_06", "AS01_09", "AS01_11", "AS01_13", "AS01_17"), missing.max),
    asi3 = asi3_phys + asi3_cog + asi3_soc,
    
    #Penn State Worry Questionnaire (PSWQ)
    across(c(PS01_01, PS01_03, PS01_08, PS01_10, PS01_11), ~ 6 - .), #reverse items
    pswq = sumscore("PS01", missing.max),
    
    #State-Trait Anxiety Inventory (STAI): Trait
    across(c(ST02_01, ST02_03, ST02_06, ST02_07, ST02_10, ST02_13, ST02_14, ST02_16, ST02_19), ~ 5 - .), #STAI trait reverse items
    stai_trait = sumscore("ST02", missing.max),
  ) %>% select(CASE:REF, datestamp, asi3, contains("asi"), pswq, stai_trait) %>% na.omit()

codes.matched = codes.matched[codes.matched %in% data.old.studi$SERIAL]

# Descriptives ------------------------------------------------------------
#data.t2 = read_rds("data_final.rds" %>% paste0(path, .))
data.t2 %>% group_by(gender, QUESTNNR) %>% summarise(N = n()) %>% arrange(gender, desc(QUESTNNR))
data.t2 %>% select(asi3, stai_trait, pswq) %>% summarise(across(.fn=list(m = mean, sd = sd)))
data.t2 %>% select(asi3, stai_trait, pswq) %>% reframe(across(.fns=range))

# Reliability -------------------------------------------------------------
#ASI
data.t2.raw %>% select(starts_with("AS01")) %>% #no reverse items
  #psych::alpha() #comprehensive analysis
  ltm::cronbach.alpha(na.rm=T) #minimal output

#PSWQ
data.t2.raw %>% select(starts_with("PS01")) %>% 
  mutate(across(c(PS01_01, PS01_03, PS01_08, PS01_10, PS01_11), ~ 6 - .)) %>% #reverse items
  #psych::alpha() #comprehensive analysis
  ltm::cronbach.alpha(na.rm=T) #minimal output

#STAI-T
data.t2.raw %>% select(starts_with("ST02")) %>% 
  mutate(across(c(ST02_01, ST02_03, ST02_06, ST02_07, ST02_10, ST02_13, ST02_14, ST02_16, ST02_19), ~ 5 - .)) %>% #reverse items
  #psych::alpha() #comprehensive analysis
  ltm::cronbach.alpha(na.rm=T) #minimal output


# Exploratory Factor Analysis ---------------------------------------------
#png("0 Factor Analysis.png" %>% paste0("C:/Users/mar84qk/Dropbox/Arbeit/Corona Psychosocial/Ergebnisse/", .))
data.t2.raw %>% select(starts_with("AS01") | starts_with("PS01") | starts_with("ST02")) %>% 
  na.omit() %>% 
  mutate(across(c(PS01_01, PS01_03, PS01_08, PS01_10, PS01_11), ~ 6 - .), #reverse score PSWQ
         across(c(ST02_01, ST02_03, ST02_06, ST02_07, ST02_10, ST02_13, ST02_14, ST02_16, ST02_19), ~ 5 - .)) %>% #reverse score STAI-T
  #factanal(., factors = ncol(.)-1, rotation="varimax") #doesn't work
  #psych::nfactors()
  cor() %>% eigenComputes() %>% nScree() %>% plotnScree() #nFactors package, see https://www.statmethods.net/advstats/factor.html#:~:text=Determining%20the%20Number%20of%20Factors%20to%20Extract
#dev.off()

# For Prevention Center ---------------------------------------------------
prevention = data.t2 %>% select(CASE, REF, QUESTNNR, gender, age, STARTED) %>% 
  left_join(data.old.studi.raw %>% select(SERIAL, B026_05, datestamp) %>% 
              rename(REF = SERIAL, Semester_old = B026_05, T1 = datestamp), by="REF") %>% 
  mutate(extraSemester = {difftime(STARTED, T1 %>% as.Date(format="%d.%m.%Y"), units="days") / 365 * 2} %>% floor() %>% as.integer(), #since T1 was at very beginning of semester, floor instead of round
         Semester = Semester_old + extraSemester) %>% select(CASE:age, Semester) %>% rename(Age = age)
prevention = "Präventionszentrum/Präventionszentrum+Demogr+LISD.tsv" %>% paste0(path, .) %>% read_tsv() %>% 
  select(CASE:Semester) %>% select(-age) %>% 
  left_join(prevention %>% select(CASE, gender:Age), by="CASE")
#prevention %>% write_excel_csv2(paste0(path, "Präventionszentrum/Präventionszentrum+Demogr+LISD_final.csv"), na="")
