library(tidyverse)

dichotomize <- function(data, column_var) {
  data %>% 
    mutate(bin = 1) %>%
    pivot_wider(
      names_from = {{ column_var }}, # https://tidyeval.tidyverse.org/dplyr.html
      values_from = bin,
      values_fill = 0,
      names_prefix = str_c(quo_name(enquo(column_var)), "_")
    )
}

#data.t2 = read_rds("data_t2.rds" %>% paste0(path, .))
#data.sfb.t0 = read_rds("data.sfb.t0.rds" %>% paste0(path, .))
data.sfb.t0 %>% select(asi3, stai_trait, pswq) %>% summarise(across(.fn=list(m = mean, sd = sd)))
data.sfb.t0 %>% select(asi3, stai_trait, pswq) %>% reframe(across(.fns=range))

#rescale strain based on baseline (T0) => changes in strain across time won't be masked
data.sfb.baselines = data.sfb.t0 %>% filter(REF %in% serials.matched) %>% summarise(
  asi3.m0 = mean(asi3, na.rm=T), asi3.sd0 = sd(asi3, na.rm=T),
  pswq.m0 = mean(pswq, na.rm=T), pswq.sd0 = sd(pswq, na.rm=T),
  stai_trait.m0 = mean(stai_trait, na.rm=T), stai_trait.sd0 = sd(stai_trait, na.rm=T))

#data.t2 %>% filter(data.t2 %>% select(cdr:lisd_trait_close) %>% rowSums() %>% is.na()) %>% View() #subjects that didn't finish the standardized questionnaires
data.t2.sfb = data.t2 %>% filter(REF %in% serials.matched) %>% mutate(time = 2) %>% 
  select(REF, time, everything()) %>% select(-CASE, -SERIAL, -QUESTNNR) %>% rename(datestamp = STARTED)
data.sfb.criteria = data.t2.sfb %>% select(REF, time, datestamp, asi3, asi3_phys:asi3_soc, pswq, stai_trait) %>% 
  #inner_join(data.sfb.t0 %>% select(REF, asi3, asi3_phys:asi3_soc, pswq, stai_trait) %>% filter(REF %in% serials.matched) %>% mutate(time = 0), by="REF") %>% 
  bind_rows(data.sfb.t0 %>% select(REF, datestamp, asi3, asi3_phys:asi3_soc, pswq, stai_trait) %>% filter(REF %in% serials.matched) %>% mutate(time = 0)) %>% 
  
  #rescale strain to T0 baseline
  mutate(strain = bind_cols((asi3-data.sfb.baselines$asi3.m0)/data.sfb.baselines$asi3.sd0, 
                            (pswq-data.sfb.baselines$pswq.m0)/data.sfb.baselines$pswq.sd0, 
                            (stai_trait-data.sfb.baselines$stai_trait.m0)/data.sfb.baselines$stai_trait.sd0) %>% rowMeans())

data.sfb.predictors = data.t2.sfb %>% select(REF, !any_of(data.sfb.criteria %>% names())) %>% select(!B002_01:TIME022) %>% 
  select(-contains("contacts_real_0"), -contains("contacts_need_0"), -contains("opinion_outgroup_")) %>% 
  inner_join(data.sfb.criteria %>% pivot_wider(id_cols=REF, names_from=time, values_from=datestamp) %>% mutate(days = difftime(`2`, `0`, units="days") %>% as.double()) %>% select(REF, days), ., by="REF")
#data.sfb.criteria %>% write_rds("data.sfb.criteria.rds" %>% paste0(path, .)); data.sfb.predictors %>% write_rds("data.sfb.predictors.rds" %>% paste0(path, .))


# Criteria Change & Interrelation ------------------------------------------
#data.sfb.criteria = read_rds("data.sfb.criteria.rds" %>% paste0(path, .)); data.sfb.predictors = read_rds("data.sfb.predictors.rds" %>% paste0(path, .))
data.sfb.wide = data.sfb.criteria %>% filter(time == 0) %>% select(-time) %>%
  inner_join(data.sfb.criteria %>% filter(time == 2) %>% select(-time),
             by = "REF", suffix = c("_t0", "_t2")) %>%
  #inner_join(data.sfb.predictors, by = "REF") %>%
  mutate(strain_t2minust0 = strain_t2 - strain_t0,
         asi3_t2minust0 = asi3_t2 - asi3_t0,
         pswq_t2minust0 = pswq_t2 - pswq_t0,
         stai_trait_t2minust0 = stai_trait_t2 - stai_trait_t0) %>% select(REF, contains("t2minust0"), everything())

{
  cat("strain_change = ", data.sfb.wide$strain_t2minust0 %>% mean(na.rm=T) %>% signif(3), ", ", sep=""); data.sfb.wide$strain_t2minust0 %>% t.test() %>% apa::t_apa(es_ci=T); cat("\n")
  cat("\n")
  cat("ASI_change = ", data.sfb.wide$asi3_t2minust0 %>% mean(na.rm=T) %>% signif(3), ", ", sep=""); data.sfb.wide$asi3_t2minust0 %>% t.test() %>% apa::t_apa(es_ci=T); cat("\n")
  cat("PSWQ_change = ", data.sfb.wide$pswq_t2minust0 %>% mean(na.rm=T) %>% signif(3), ", ", sep=""); data.sfb.wide$pswq_t2minust0 %>% t.test() %>% apa::t_apa(es_ci=T); cat("\n")
  cat("STAI_change = ", data.sfb.wide$stai_trait_t2minust0 %>% mean(na.rm=T) %>% signif(3), ", ", sep=""); data.sfb.wide$stai_trait_t2minust0 %>% t.test() %>% apa::t_apa(es_ci=T); cat("\n")
  cat("\n")
  cat("ASI  & PSWQ:   "); data.sfb.wide %>% with(cor.test(asi3_t2minust0, pswq_t2minust0)) %>% apa::cor_apa(r_ci=T); cat("\n")
  cat("ASI  & STAI-T: "); data.sfb.wide %>% with(cor.test(asi3_t2minust0, stai_trait_t2minust0)) %>% apa::cor_apa(r_ci=T); cat("\n")
  cat("PSWQ & STAI-T: "); data.sfb.wide %>% with(cor.test(pswq_t2minust0, stai_trait_t2minust0)) %>% apa::cor_apa(r_ci=T); cat("\n")
}



# Lack of Variance --------------------------------------------------------
data.sfb.predictors %>% pull(gender) %>% summary()
data.sfb.predictors %>% pull(age) %>% summary()

data.sfb.predictors %>% pull(vaccination) %>% summary() #no differentiation between non fully vaccinated
data.sfb.predictors %>% pull(vaccination_full) %>% mean() #only <9% not fully vaccinated => self-selection in sample? (could compare to vaccinations in this age group by timestamp)
data.sfb.predictors %>% pull(infection_self) %>% {. > 0} %>% mean() #only 6% of sample was infected at least once
data.sfb.predictors %>% pull(risk_group) %>% table(useNA="ifany") #possible as continuous predictor?
#data.sfb.predictors %>% pull(covid_risk) %>% summary()
data.sfb.predictors %>% pull(covid_media_update) %>% table(useNA="ifany") #possible as continuous predictor?

data.sfb.predictors %>% pull(living_alone) %>% mean() #16% living alone
data.sfb.predictors %>% pull(unemployed) %>% mean(na.rm=T) #<2% unemployed
data.sfb.predictors %>% pull(occupation) %>% summary() #{summary(.) / length(.)}

data.sfb.predictors %>% pull(city) %>% summary()


# Prepare for Elastic Net -------------------------------------------------
#data.sfb.predictors %>% with(cor.test(covid_psychological_change_emo, covid_psychological_change_sleep)) %>% apa::cor_apa()
#emotional change and sleep are moderately correlated => treat separately

data.sfb.predictors.elastic = data.sfb.predictors %>% 
  select(-vaccination, #use vaccination_full instead
         -unemployed,  #use occupation instead
         -contains("living_with") #use living_alone instead
  ) %>%  filter(gender != "diverse", #too little variance => dichotomize
                occupation %in% c("employed", "studies/apprenticeship")) %>% #too little variance => employed vs. studies/apprenticeship 
  mutate(female = if_else(gender == "female", 1, 0),
         employed = if_else(occupation == "employed", 1, 0),
         vaccination_full = if_else(vaccination_full, 1, 0),
         living_alone = if_else(living_alone, 1, 0),
         city100K = if_else(city == "city > 100K", 1, 0),
         german = if_else(citizenship == "german", 1, 0)) %>% 
  select(-gender, -occupation, -city, -citizenship, -days) %>% 
  select(-covid_psychological_change) %>% 
  #select(-contains("covid_psychological_change_")) %>% 
  dichotomize(covid_media_update)
  
#data.sfb.predictors.elastic %>% write_rds("data.sfb.data.sfb.predictors.elastic.rds" %>% paste0(path, .))

#network graphs of predictors as overview
data.sfb.predictors.elastic %>% select(-REF, -contains("cdr_"), -contains("contacts_real_")) %>% corrr::correlate() %>% corrr::network_plot(colors=c("red", "green"))
data.sfb.predictors.elastic %>% select(-REF, -contains("cdr_"), -contains("contacts_real_")) %>% 
  select(cdr:lisd_trait_close, 
         contains("covid_psychological_change"), covid_situation_change, covid_efficiency_change, 
         daystructure, covid_worry_livelihood, contacts_need,
         covid_risk, covid_behavior_change, safety_behaviour) %>% corrr::correlate() %>% corrr::network_plot(colors=c("red", "green"))

# Elastic Net -------------------------------------------------------------
#data.sfb.predictors.elastic = read_rds("data.sfb.data.sfb.predictors.elastic.rds" %>% paste0(path, .))
#data.sfb.criteria = read_rds("data.sfb.criteria.rds" %>% paste0(path, .))

data <- data.sfb.criteria %>% 
  filter(time == 0) %>%
  select(-time) %>%
  inner_join(
    data.sfb.criteria %>% 
      filter(time == 2) %>%
      select(-time),
    by = "REF",
    suffix = c("_t0", "_t2")
  ) %>%
  inner_join(
    data.sfb.predictors.elastic,
    by = "REF"
  ) %>% 
  mutate(
    strain_t2minust0 = strain_t2 - strain_t0,
    asi3_t2minust0 = asi3_t2 - asi3_t0,
    pswq_t2minust0 = pswq_t2 - pswq_t0,
    stai_trait_t2minust0 = stai_trait_t2 - stai_trait_t0,
  ) %>%
  filter(
    !is.na(strain_t2minust0),
    !is.na(asi3_t2minust0), 
    !is.na(pswq_t2minust0), 
    !is.na(stai_trait_t2minust0)
  ) %>%
  rename_all(~ gsub(" ", "_", .))

set.seed(1)
data_split <- rsample::initial_split(data)
data_train <- rsample::training(data_split)
data_test <- rsample::testing(data_split)

get_bics <- function(target) {
  # https://www.jihongzhang.org/post/2019-02-19-lasso-regression-with-glmnet/
  X <- data_train %>%
    select(-REF, -ends_with("t0"), -ends_with("t2")) %>%
    as.matrix()
  y <- data_train %>%
    select(target) %>%
    as.matrix()
  
  fit <- glmnet::glmnet(
    X, 
    y, 
    alpha = 0.5,
    nlambda = 2000
  )
  
  #plot(fit, xvar = "lambda", label = TRUE)
  
  order <- fit$beta %>%
    as.matrix() %>%
    as_tibble(rownames = "predictor") %>%
    mutate_if(is.numeric, ~ if_else(. == 0, 0, 1)) %>%
    mutate(lambda_resistance = rowSums(select_if(., is.numeric))) %>%
    arrange(desc(lambda_resistance)) %>%
    select(predictor, lambda_resistance)
  
  get_model_values <- function(data.sfb.predictors.elastic) {
    c(
      target,
      "~",
      data.sfb.predictors.elastic %>%
        str_flatten(collapse = "+")
    ) %>%
      str_flatten() %>%
      as.formula() %>%
      lm(data = data) %>%
      broom::glance() %>%
      select(BIC, adj.r.squared, p.value)
  }
  
  order %>% 
    mutate(
      row = row_number(),
      data.sfb.predictors.elastic = map(row, ~ order$predictor[1:.]),
      model_values = map(data.sfb.predictors.elastic, get_model_values),
      bic = map_dbl(model_values, ~ .$BIC),
      r_sq_adj = map_dbl(model_values, ~ .$adj.r.squared),
      p = map_dbl(model_values, ~ .$p.value),
      bic_diff_to_min = bic - min(bic)
    ) %>%
    select(-row, -data.sfb.predictors.elastic, -model_values) %>% 
    rowwise() %>%
    mutate(
      correlation_with_target_r = cor(
        data_train %>% pull(target), 
        data_train %>% pull(predictor),
        use="pairwise.complete.obs"
      ),
      correlation_with_target_p = cor.test(
        data_train %>% pull(target), 
        data_train %>% pull(predictor)
      )$p.value
    ) %>%
    ungroup()
}

results <- tribble(
  ~target,
  "strain_t2minust0",
  "asi3_t2minust0",
  "pswq_t2minust0",
  "stai_trait_t2minust0"
) %>%
  mutate(
    bics = map(target, get_bics)
  )

View(results)
#results %>% write_rds("elastic 1 (covid_psy combined).rds" %>% paste0(path, .))
#results %>% write_rds("elastic 2 (covid_psy sep).rds" %>% paste0(path, .))
#for (i in results$target %>% seq_along()) { results[[i, "bics"]][[1]] %>% clipr::write_clip(); readline(paste0(results[i, 1], " copied to clipboard (press enter to proceed)")) }
#TODO try https://www.r-bloggers.com/2019/08/creating-excel-workbooks-with-multiple-sheets-in-r/ ?


# Twin Matching (quasi-longitudinal) --------------------------------------
matching.vars = results[[1, "bics"]][[1]]$predictor[1:which(results[[1, "bics"]][[1]]$bic_diff_to_min==0)] %>% 
  c("gender", "age", .)

# data.sfb.predictors.overview = data.sfb.predictors %>% select(matching.vars) %>% pivot_longer(cols=-gender, names_to="predictor") %>% 
#   group_by(predictor, gender) %>% summarise(M = mean(value, na.rm=T), SD = sd(value, na.rm=T), `SD/2` = SD/2, `SD/3` = SD/3, `SD/4` = SD/4)
# for (var in data.sfb.predictors %>% select(matching.vars) %>% names() %>% setdiff("gender") %>% sort()) {
#   print(data.sfb.predictors %>% ggplot(aes_string(x=var)) + 
#           geom_histogram(color="black", fill="grey") + 
#           #geom_errorbarh(data=data.sfb.predictors.overview %>% filter(predictor==var), aes(xmin=M-SD, xmax=M+SD)) +
#           facet_wrap(~gender) +
#           myGgTheme)
# }
# data.sfb.predictors.overview

data.t2.twins = data.t2 %>% filter(gender %>% is.na() == F) %>% 
  select(CASE:age, asi3, pswq, stai_trait, strain, matching.vars) %>% na.omit() %>% 
  mutate(#toMatch = REF %in% {matching.sfb %>% filter(datestamp1 %>% is.na()) %>% pull(REF)}, #only match SFB subjects that have no T1 (but different time point of survey)
    toMatch = QUESTNNR=="sfb"
  ) %>% 
  #filter(QUESTNNR=="base" | toMatch) #only keep studi sample and those of SFB that need twin matching for T1
  filter(REF %in% codes.matched | toMatch) #only keep matched studi sample (i.e., T1 data exist) and those of SFB that need twin matching for T1
data.t2.twins %>% group_by(gender, QUESTNNR, toMatch) %>% summarise(N = n()) %>% arrange(gender, desc(QUESTNNR))

predictors.overview = data.t2.twins %>% select(matching.vars) %>% pivot_longer(cols=-gender, names_to="predictor") %>% 
  group_by(predictor, gender) %>% summarise(M = mean(value, na.rm=T), SD = sd(value, na.rm=T), `2/3 SD` = SD*2/3, `.55*SD` = SD*.55, `SD/2` = SD/2, `SD/3` = SD/3, `SD/4` = SD/4)
#predictors.overview %>% mutate(across(contains("SD"), floor))
for (var in data.t2.twins %>% select(matching.vars) %>% names() %>% setdiff("gender") %>% sort()) {
  print(data.t2.twins %>% ggplot(aes_string(x=var, fill="toMatch")) + 
          geom_histogram(color="black", position="dodge", binwidth=1) + 
          #geom_errorbarh(data=predictors.overview %>% filter(predictor==var), aes(xmin=M-SD, xmax=M+SD)) +
          facet_wrap(~gender) +
          myGgTheme)
}
predictors.overview %>% filter(gender!="diverse") %>% mutate(across(contains("SD"), floor))

#caliper vs. inclusion
data.matched.list = list()
matches.indices.list = list()
matches.sorted.list = list()
matching.list = list()
matches.n.list = list()
for (c in seq(.25, 1, .05)) {
  
  data.matched = tibble()
  matches.indices = tibble()
  matches.sorted = tibble()
  #don't use Matching::Matchby because it offers fewer options => explicitly loop
  for (g in data.t2.twins %>% filter(toMatch) %>% pull(gender) %>% unique()) {
    data.t2.twins.gender = data.t2.twins %>% filter(gender==g) %>% 
      #mutate(asi3.z = scale(asi3) %>% as.numeric(), stai_trait.z = scale(stai_trait) %>% as.numeric(), pswq.z = scale(pswq) %>% as.numeric(), strain.z = scale(strain) %>% as.numeric(), age.z = scale(age) %>% as.numeric()) %>% 
      mutate(across(matching.vars %>% setdiff("gender") %>% c("asi3", "stai_trait", "pswq", "strain", .), 
                    function(x) {x %>% scale() %>% as.numeric()}, .names="{.col}.z"))
    
    
    matching = Matching::Match( #TODO check Matching::GenMatch to calculate weights via genetic search?
      Tr=data.t2.twins.gender %>% pull(toMatch), 
      #X=data.t2.twins.gender %>% select(strain),
      #X=data.t2.twins.gender %>% select(strain.z), #scaling does not seem to make a difference
      #X=data.t2.twins.gender %>% select(asi3, stai_trait, pswq),
      #X=data.t2.twins.gender %>% select(asi3.z, stai_trait.z, pswq.z),
      X = data.t2.twins.gender %>% select(strain, matching.vars %>% setdiff(c("gender"))),
      #exact = c(F, F, F, T),
      #exact = data.t2.twins.gender %>% select(strain, matching.vars %>% setdiff(c("gender"))) %>% names() == "covid_psychological_change_emo",
      caliper=c, #maximum distance (in SDs) for a match
      #TODO stricter caliper for strain?
      ties=T #shall ties be broken deterministically?
    )
    
    #TODO apply matching (provide weights if used above)
    #Matching::MatchBalance(toMatch ~ asi3 + stai_trait + pswq, data=data.t2.twins.gender, match.out=matching)
    
    if (matching %>% is.na() %>% any() == F) {
      matching.df = bind_cols(matching$index.treated, matching$index.control) %>% rename(toMatch = "...1", matched = "...2") %>% 
        group_by(toMatch) %>% mutate(matchNum = paste0("match", 1:n())) %>% 
        pivot_wider(names_from = "matchNum", values_from = "matched")
      matches.indices = matches.indices %>% bind_rows(matching.df)
      
      matches.sorted = matches.sorted %>% bind_rows(data.t2.twins.gender %>% slice(matching.df %>% t() %>% c() %>% na.omit()))
      #View(matches.sorted) #peak into individual values of matches
      #if (writeFiles) matches.sorted %>% write_excel_csv2(., paste0(path, "Matching_", g, ".csv"))
      
      data.matched = data.t2.twins.gender %>% slice(matching.df$toMatch) %>% mutate(
        match = data.t2.twins.gender %>% slice(matching.df$match1) %>% pull(REF)) %>% #only first match (within ties)
        bind_rows(data.matched, .)
    }
  }
  
  matches.n = data.t2.twins %>% filter(toMatch) %>% group_by(gender) %>% summarise(N = n()) %>% mutate(role = "toMatch") %>% 
    bind_rows(matches.sorted %>% group_by(gender) %>% summarise(N = n()/2) %>% mutate(role = "matched")) %>% 
    select(gender, role, everything()) %>% arrange(gender, desc(role)) %>% 
    pivot_wider(names_from=role, values_from=N, values_fill=0) %>% mutate(`%` = matched/toMatch)
  
  data.matched.list[[c %>% as.character()]] = data.matched
  matches.indices.list[[c %>% as.character()]] = matches.indices
  matches.sorted.list[[c %>% as.character()]] = matches.sorted
  matching.list[[c %>% as.character()]] = matching
  matches.n.list[[c %>% as.character()]] = matches.n
}

#caliper vs. average deviation (in SDs)
matching.deviations.list = list()
for (c in data.matched.list %>% names()) {
  data.matched = data.matched.list[[c]] %>% select(CASE:REF, match, everything()) %>% select(-toMatch) #%>% select(-contains(".z"))
  matching.deviations = data.matched %>% 
    left_join(
      data.t2.twins %>% group_by(gender) %>% 
        mutate(across(matching.vars %>% setdiff("gender") %>% c("asi3", "stai_trait", "pswq", "strain", .), 
                      function(x) {x %>% scale() %>% as.numeric()}, .names="{.col}.z")) %>% 
        select(data.matched %>% select(-match) %>% names()),
      by=c("match" = "REF")) %>% 
    mutate(gender = gender.x == gender.y,
           strain.z = abs(strain.z.x - strain.z.y), strain = abs(strain.x - strain.y),
           age.z = abs(age.z.x - age.z.y), age = abs(age.x - age.y), 
           adsk.z = abs(adsk.z.x - adsk.z.y), adsk = abs(adsk.x - adsk.y),
           ius_i.z = abs(ius_i.z.x - ius_i.z.y), ius_i = abs(ius_i.x - ius_i.y), 
           phq2.z = abs(phq2.z.x - phq2.z.y), phq2 = abs(phq2.x - phq2.y),
           covid_psychological_change_emo.z = abs(covid_psychological_change_emo.z.x - covid_psychological_change_emo.z.y),
           covid_psychological_change_emo = abs(covid_psychological_change_emo.x - covid_psychological_change_emo.y)) %>% 
    select(!contains(".x") & !contains(".y"), everything())
  
  matching.deviations.list[[c]] = matching.deviations %>% select(ends_with(".z"), everything()) %>% 
    pivot_longer(strain.z:covid_psychological_change_emo.z, 
                 names_to="variable", values_to="value") %>% select(REF, match, variable, value, everything()) %>% 
    #group_by(gender.x) %>% 
    summarise(across(value, .fns=list(m = mean, sd = sd, min=min, max=max)))
}
matching.deviations.list %>% bind_rows(.id="caliper") %>% mutate(caliper = caliper %>% as.numeric()) %>% 
  ggplot(aes(x=caliper, y=value_m, group=NA)) +
  #matching.deviations.list %>% bind_rows(.id="caliper") %>% ggplot(aes(x=caliper, y=value_m, color=gender.x, group=gender.x)) +
  geom_ribbon(aes(ymin=value_m-value_sd, ymax=value_m+value_sd), alpha=.2) +
  geom_line(size=2) + geom_point(size=4) + 
  scale_color_viridis_d() + myGgTheme + ylab("Average Deviation (SDs)")

#caliper vs. subjects retained (percentage)
matches.n.list %>% bind_rows(.id = "caliper") %>% mutate(caliper = caliper %>% as.numeric()) %>% 
  ggplot(aes(x=caliper, y=`%`, color=gender, group=gender)) +
  geom_line(size=2) + geom_point(size=4) + 
  scale_color_viridis_d() + myGgTheme + ylab("Subjects Retained (%)")
#ggsave("0 Caliper.png" %>% paste0("C:/Users/mar84qk/Dropbox/Arbeit/Corona Psychosocial/Ergebnisse/", .), width=1920, height=1080, units="px")

#caliper vs. subjects retained (total)
# matches.n.list %>% bind_rows(.id = "caliper") %>% ggplot(aes(x=caliper, y=matched, color=gender, group=gender)) +
#   geom_line(size=2) + geom_point(size=4) + 
#   scale_color_viridis_d() + myGgTheme + ylab("Subjects Retained Total")
matches.n.list %>% bind_rows(.id = "caliper") %>% filter(caliper %in% c("0.6", "0.7")) #%>% mutate(lost = toMatch - matched, lost_p = 1 - `%`)

#select results from one caliper
data.matched = data.matched.list[["0.7"]]
matches.indices = matches.indices.list[["0.7"]]
matches.sorted = matches.sorted.list[["0.7"]]
matching = matching.list[["0.7"]]
matches.n = matches.n.list[["0.7"]]

#diagnostics
matches.sorted %>% pivot_longer(age:covid_psychological_change_emo, names_to="variable") %>% 
  group_by(variable, gender, toMatch) %>% summarise(M = mean(value, na.rm=T)) %>% pivot_wider(names_from="variable", values_from="M")

#View(matches.sorted)
#View(data.matched)
#View(matches.indices)

#TODO custom diagnostics function why cases cannot get matched? Rerun without caliper, calculate univariate distances between matches in SDs, count how often distance > threshold. (But it could be that other match would be found with worse average distance but all distances < threshold)
data.matched = data.matched %>% select(CASE:REF, match, everything()) %>% select(-toMatch) #%>% select(-contains(".z"))
matching.deviations = data.matched %>% 
  left_join(
    data.t2.twins %>% group_by(gender) %>% 
      mutate(across(matching.vars %>% setdiff("gender") %>% c("asi3", "stai_trait", "pswq", "strain", .), 
                    function(x) {x %>% scale() %>% as.numeric()}, .names="{.col}.z")) %>% 
      select(data.matched %>% select(-match) %>% names()),
    by=c("match" = "REF")) %>% 
  mutate(gender = gender.x == gender.y,
         strain.z = abs(strain.z.x - strain.z.y), strain = abs(strain.x - strain.y),
         age.z = abs(age.z.x - age.z.y), age = abs(age.x - age.y), 
         adsk.z = abs(adsk.z.x - adsk.z.y), adsk = abs(adsk.x - adsk.y),
         ius_i.z = abs(ius_i.z.x - ius_i.z.y), ius_i = abs(ius_i.x - ius_i.y), 
         phq2.z = abs(phq2.z.x - phq2.z.y), phq2 = abs(phq2.x - phq2.y),
         covid_psychological_change_emo.z = abs(covid_psychological_change_emo.z.x - covid_psychological_change_emo.z.y),
         covid_psychological_change_emo = abs(covid_psychological_change_emo.x - covid_psychological_change_emo.y)) %>% 
  select(!contains(".x") & !contains(".y"), everything())

#matching.deviations %>% summarise(gender = mean(gender==F), across(age:covid_psychological_change_emo, list(m = mean, sd = sd)))
#matching.deviations %>% summarise(gender = mean(gender==F))
matching.deviations %>% pivot_longer(strain.z:covid_psychological_change_emo, names_to="variable", values_to="value") %>% 
  group_by(variable, gender.x) %>% summarise(#gender = mean(gender==F), 
    across(value, .fns=list(m = mean, sd = sd, min=min, max=max))) %>% View("matching.deviations")

matching.deviations %>% select(ends_with(".z"), everything()) %>% 
  pivot_longer(strain.z:covid_psychological_change_emo.z, 
               names_to="variable", values_to="value") %>% select(REF, match, variable, value, everything()) %>% 
  group_by(gender.x) %>% summarise(across(value, .fns=list(m = mean, sd = sd, min=min, max=max))) #%>% View("matching.deviations mean")

matching.deviations %>% select(ends_with(".z"), everything()) %>% 
  pivot_longer(strain.z:covid_psychological_change_emo.z, 
               names_to="variable", values_to="value") %>% 
  ggplot(aes(x=value, fill=gender.x)) + 
  #geom_histogram(color="black", position=position_dodge()) + 
  facet_wrap(~variable) + 
  geom_histogram(aes(y = stat(density)), color="black", position=position_dodge()) + 
  myGgTheme + xlab("deviation (SDs)") #+ scale_y_continuous(labels = scales::percent)

for (v in matching.vars %>% setdiff(("gender"))) {
  cat(v, ": ")
  cor.test(matching.deviations %>% pull(paste0(v, ".x")), matching.deviations %>% pull(paste0(v, ".y"))) %>% apa::cor_apa()
  cat("\n")
}

#data.matched %>% write_rds("data.matched.rds" %>% paste0(path, .))



# Quasi-Longitudinal ------------------------------------------------------
#data.matched = read_rds("data.matched.rds" %>% paste0(path, .))

data.matched.analysis = data.matched %>% select(-CASE, -SERIAL, -QUESTNNR, -strain.z) %>% rename_with(paste0, .cols=-c("REF", "match", "gender"), "_t2") %>% 
  left_join( #add timestamp of t2
    data.t2 %>% select(REF, STARTED) %>% rename(datestamp_t2 = STARTED), 
    by="REF") %>% 
  left_join( #add t0
    data.sfb.t0 %>% select(REF, datestamp, asi3, pswq, stai_trait, strain,
                           spai, lsas, gse, cerq_adapt, cerq_maladapt, cerq_acceptance, ctq, ctq_problem, lte, le) %>% 
                           #contains(ctq) #to analyze subscales of ctq
      rename_with(paste0, .cols=-REF, "_t0"), 
    by="REF") %>% 
  left_join( #add t1
    data.old.studi %>% select(SERIAL, asi3, pswq, stai_trait) %>% rename(match = SERIAL) %>% 
      mutate(strain = bind_cols((asi3-data.sfb.baselines$asi3.m0)/data.sfb.baselines$asi3.sd0, 
                                (pswq-data.sfb.baselines$pswq.m0)/data.sfb.baselines$pswq.sd0, 
                                (stai_trait-data.sfb.baselines$stai_trait.m0)/data.sfb.baselines$stai_trait.sd0) %>% rowMeans()) %>% 
      rename_with(paste0, .cols=-match, "_t1"), by="match") %>% 
  mutate(days = difftime(datestamp_t2, datestamp_t0, units="days") %>% as.double())

#check missing variables
data.matched.analysis %>% filter(if_any(.fns=is.na)) %>% View("Missing Variables T0")
#data.matched.analysis %>% arrange(datestamp_t0) %>% select(contains("t0")) %>% View("Missing Variables T0 fixed")


data.matched.analysis.long = data.matched.analysis %>% pivot_longer(contains("strain"), names_to="time", values_to="strain") %>% 
  mutate(time = time %>% gsub("strain_t", "", .) %>% as.integer()) %>% 
  select(REF, match, time, strain, everything()) %>% arrange(REF, time) %>% 
  mutate(time = case_when(time == 0 ~ "Pre",
                          time == 1 ~ "Peak",
                          time == 2 ~ "Downturn",
                          T ~ "Error") %>% as_factor())

data.matched.analysis.long.scaled = data.matched.analysis.long %>% 
  #z-transform continuous predictors
  mutate(across(.cols=c("age", "adsk", "ius_i", "phq2", "covid_psychological_change_emo", "phq2") %>% paste0("_t2") %>% 
                  c({c("spai", "lsas", "gse", "cerq_adapt", "cerq_maladapt", "cerq_acceptance", "ctq", "lte", "le") %>% paste0("_t0")}) %>% 
                  c("days"), 
                .fns=function(x) {scale(x)[,1]}),
         time = time %>% as.factor())



# Sample Descriptives -----------------------------------------------------
#DVs
data.matched.analysis %>% summarise(across(contains("asi3_") | contains("pswq_") | contains("stai_trait_"), 
                                           mean, na.rm=T, .names="{.col}.m")) %>% select(order(colnames(.)))
data.matched.analysis %>% select(contains("asi3_") | contains("pswq_") | contains("stai_trait_")) %>% select(order(names(.))) %>% select(contains("t0"), contains("t1"), contains("t2")) %>% corrr::correlate()

#T0 predictors
data.matched.analysis %>% group_by(gender) %>% summarise(N = n())
data.matched.analysis %>% summarise(age_t2.m = mean(age_t2, na.rm=T), age_t2.sd = sd(age_t2, na.rm=T))

data.matched.analysis %>% summarise(spai_t0.m = mean(spai_t0*22, na.rm=T), spai_t0.sd = sd(spai_t0*22, na.rm=T), #*22 for sum score instead of mean item score
                                    lsas_t0.m = mean(lsas_t0, na.rm=T), lsas_t0.sd = sd(lsas_t0, na.rm=T))
data.matched.analysis %>% summarise(gse_t0.m = mean(gse_t0, na.rm=T), gse_t0.sd = sd(gse_t0, na.rm=T), z = (gse_t0.m - 29.43) / 5.36)
data.matched.analysis %>% summarise(across(contains("cerq"), mean, na.rm=T, .names="{.col}.m"))
data.matched.analysis %>% summarise(ctq_t0.m = mean(ctq_t0, na.rm=T), ctq_t0.sd = sd(ctq_t0, na.rm=T),
                                    ctq_problem = mean(ctq_problem_t0))
data.matched.analysis %>% summarise(lte_t0.m = mean(lte_t0, na.rm=T), lte_t0.sd = sd(lte_t0, na.rm=T))
data.matched.analysis %>% summarise(le_t0.m = mean(le_t0, na.rm=T), le_t0.sd = sd(le_t0, na.rm=T))


#T2 predictor
data.matched.analysis %>% summarise(phq2_t2.m = mean(phq2_t2), phq2_t2.sd = sd(phq2_t2), cutoff = mean(phq2_t2 >= 3))



# Selective Attrition -----------------------------------------------------
df = data.sfb.t0 %>% mutate(stayed = REF %in% {data.matched.analysis %>% pull(REF)})
df %>%  group_by(stayed) %>% summarise(n = n(),
                                       across(c("asi3", "pswq", "stai_trait", "spai", "lsas", "gse", "cerq_maladapt", "cerq_adapt", "cerq_acceptance", "ctq", "lte", "le"),
                                              .fn=mean, na.rm=T))
                                              #.fn=list(m = "mean", sd = "sd"), na.rm=T)) %>% View("T0 Descriptives")
attrition = tibble(var = c("asi", "pswq", "stai_trait", "spai", "lsas", "gse", "cerq_maladapt", "cerq_adapt", "cerq_acceptance", "ctq", "lte", "le")) %>% 
  rowwise() %>% 
  mutate(test = t.test(df %>% filter(stayed==T) %>% pull(var),
                       df %>% filter(stayed==F) %>% pull(var), 
                       var.equal=T) %>% list(),
         p = test$p.value,
         sig = case_when(p < .001 ~ "***", p < .01 ~ "**", p < .05 ~ "*", p < .1 ~ "#", p < .2 ~ "##", T ~ ""),
         apa = test %>% apa::t_apa(print=F, es_ci=T)) %>% ungroup()
attrition
attrition %>% filter(p == min(p))

# Analyses ----------------------------------------------------------------
#time main effect
data.matched.analysis.long.scaled %>% 
  afex::aov_ez(id="REF", dv="strain", .,
               within=c("time"),
               between=c("gender", "age_t2", "days"),
               observed=c("gender", "age_t2", "days"),
               factorize=F) %>% apa::anova_apa(force_sph_corr=T)
# data.matched.analysis.long %>% group_by(time, gender) %>% summarise(strain.se = se(strain), strain = mean(strain)) %>% 
#   ggplot(aes(x=time, y=strain, color=gender, group=gender)) + 
#   geom_hline(yintercept=0, linetype="dashed") +
#   geom_violin(data=data.matched.analysis.long, aes(group=interaction(time, gender))) +
#   geom_errorbar(aes(ymin=strain-strain.se, ymax=strain+strain.se), size=2, width=.25, position=position_dodge(width=.9)) +
#   geom_line(size=2, position=position_dodge(width=.9)) + geom_point(size=5, position=position_dodge(width=.9)) + myGgTheme +
#   scale_x_discrete("", labels = c("Pre","Peak","Downturn")) + scale_color_discrete("Gender") + ylab("Strain")
data.matched.analysis.long %>% group_by(time, gender) %>% summarise(strain.se = se(strain), strain = mean(strain)) %>% 
  ggplot(aes(x=time, y=strain, color=gender, group=gender)) + 
  geom_hline(yintercept=0, linetype="dashed") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1, data=data.matched.analysis.long, aes(group=interaction(time, gender)), position=position_dodge(width=.9)) +
  #geom_violin(aes(group=interaction(time, gender)), data=data.matched.analysis.long) +
  geom_violin(aes(group=interaction(time, gender), fill=gender), color="black", alpha=.5, data=data.matched.analysis.long) +
  geom_boxplot(width=0.1, position=position_dodge(width=.9), data=data.matched.analysis.long, aes(group=interaction(time, gender))) +
  geom_errorbar(aes(ymin=strain-strain.se*1.96, ymax=strain+strain.se*1.96), size=2, width=.25) +
  geom_line(size=2) + geom_point(size=5) + myGgTheme +
  scale_x_discrete("", labels = c("Pre","Peak","Downturn")) + ylab("Strain") + scale_color_discrete("Gender") + scale_fill_discrete("Gender")
#ggsave("1 Time & Gender.png" %>% paste0("C:/Users/mar84qk/Dropbox/Arbeit/Corona Psychosocial/Ergebnisse/", .), width=1920, height=1080, units="px")
with(data.matched.analysis, t.test(x=strain_t1, y=strain_t0, paired=T) %>% apa::t_apa(es_ci=T))
with(data.matched.analysis, t.test(x=strain_t2, y=strain_t1, paired=T) %>% apa::t_apa(es_ci=T))
with(data.matched.analysis, t.test(x=strain_t2, y=strain_t0, paired=T) %>% apa::t_apa(es_ci=T))
data.matched.analysis.long %>% group_by(gender) %>% summarise(strain = mean(strain))

#Social Phobia and Anxiety Inventory (SPAI)
data.matched.analysis.long.scaled %>% 
  afex::aov_ez(id="REF", dv="strain", .,
               within=c("time"),
               between=c("gender", "age_t2", "days", "spai_t0"),
               observed=c("gender", "age_t2", "days", "spai_t0"),
               factorize=F) %>% apa::anova_apa(force_sph_corr=T)
data.matched.analysis.long = data.matched.analysis.long %>% mutate(spai_group = ifelse(spai_t0 > median(spai_t0), "high", "low") %>% factor(levels=c("high", "low")))
plot.spai = data.matched.analysis.long %>% 
  group_by(time, spai_group) %>% summarise(strain.se = se(strain), strain = mean(strain)) %>% 
  ggplot(aes(x=time, y=strain, color=spai_group, group=spai_group)) + 
  geom_hline(yintercept=0, linetype="dashed") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1, data=data.matched.analysis.long, aes(group=interaction(time, spai_group)), position=position_dodge(width=.9)) +
  #geom_violin(aes(group=interaction(time, spai_group)), data=data.matched.analysis.long) +
  geom_violin(aes(group=interaction(time, spai_group), fill=spai_group), color="black", alpha=.5, data=data.matched.analysis.long) +
  geom_boxplot(width=0.1, position=position_dodge(width=.9), data=data.matched.analysis.long, aes(group=interaction(time, spai_group))) +
  geom_errorbar(aes(ymin=strain-strain.se*1.96, ymax=strain+strain.se*1.96), size=2, width=.25) +
  geom_line(size=2) + geom_point(size=5) + myGgTheme + ggtitle("SPAI") +
  scale_x_discrete(labels = c("Pre","Peak","Downturn")) + scale_color_discrete("SPAI") + scale_fill_discrete("SPAI") + ylab("Strain") + xlab("")
plot.spai
#ggsave("2.1 SPAI.png" %>% paste0("C:/Users/mar84qk/Dropbox/Arbeit/Corona Psychosocial/Ergebnisse/", .), width=1920, height=1080, units="px")
data.matched.analysis %>% select(spai_t0, strain_t0, strain_t1, strain_t2) %>% corrr::correlate()
#data.matched.analysis.long %>% group_by(gender) %>% summarise(spai_t0 = mean(spai_t0)) #gender main effect
data.matched.analysis.long %>% group_by(REF, gender, spai_t0) %>% summarise(strain = mean(strain)) %>% 
  group_by(gender) %>% summarise(r = cor.test(strain, spai_t0) %>% apa::cor_apa())
#data.matched.analysis %>% transmute(spai_t0 = spai_t0, slope01 = strain_t1-strain_t0, slope12 = strain_t2-strain_t1) %>% corrr::correlate()
with(data.matched.analysis %>% mutate(slope = strain_t1-strain_t0), cor.test(spai_t0, slope) %>% apa::cor_apa(r_ci=T))
with(data.matched.analysis %>% mutate(slope = strain_t2-strain_t1), cor.test(spai_t0, slope) %>% apa::cor_apa(r_ci=T))

#Liebowitz Social Anxiety Scale (LSAS)
data.matched.analysis.long.scaled %>% 
  afex::aov_ez(id="REF", dv="strain", .,
               within=c("time"),
               between=c("gender", "age_t2", "days", "lsas_t0"),
               observed=c("gender", "age_t2", "days", "lsas_t0"),
               factorize=F) %>% apa::anova_apa(force_sph_corr=T)
data.matched.analysis.long = data.matched.analysis.long %>% mutate(lsas_group = ifelse(lsas_t0 > median(lsas_t0), "high", "low") %>% factor(levels=c("high", "low")))
plot.lsas = data.matched.analysis.long %>% 
  group_by(time, lsas_group) %>% summarise(strain.se = se(strain), strain = mean(strain)) %>% 
  ggplot(aes(x=time, y=strain, color=lsas_group, group=lsas_group)) + 
  geom_hline(yintercept=0, linetype="dashed") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1, data=data.matched.analysis.long, aes(group=interaction(time, lsas_group)), position=position_dodge(width=.9)) +
  #geom_violin(aes(group=interaction(time, lsas_group)), data=data.matched.analysis.long) +
  geom_violin(aes(group=interaction(time, lsas_group), fill=lsas_group), color="black", alpha=.5, data=data.matched.analysis.long) +
  geom_boxplot(width=0.1, position=position_dodge(width=.9), data=data.matched.analysis.long, aes(group=interaction(time, lsas_group))) +
  geom_errorbar(aes(ymin=strain-strain.se*1.96, ymax=strain+strain.se*1.96), size=2, width=.25) +
  geom_line(size=2) + geom_point(size=5) + myGgTheme + ggtitle("LSAS") +
  scale_x_discrete(labels = c("Pre","Peak","Downturn")) + scale_color_discrete("LSAS") + scale_fill_discrete("LSAS") + ylab("Strain") + xlab("")
plot.lsas
#ggsave("2.2 LSAS.png" %>% paste0("C:/Users/mar84qk/Dropbox/Arbeit/Corona Psychosocial/Ergebnisse/", .), width=1920, height=1080, units="px")
data.matched.analysis %>% select(lsas_t0, strain_t0, strain_t1, strain_t2) %>% corrr::correlate()
data.matched.analysis.long %>% group_by(gender) %>% summarise(lsas_t0 = mean(lsas_t0)) #n.s.
with(data.matched.analysis %>% mutate(slope = strain_t1-strain_t0), cor.test(lsas_t0, slope) %>% apa::cor_apa(r_ci=T))
with(data.matched.analysis %>% mutate(slope = strain_t2-strain_t1), cor.test(lsas_t0, slope) %>% apa::cor_apa(r_ci=T))

with(data.matched.analysis, cor.test(spai_t0, lsas_t0) %>% apa::cor_apa(r_ci=T))

#General Self Efficacy (GSE)
data.matched.analysis.long.scaled %>% 
  afex::aov_ez(id="REF", dv="strain", .,
               within=c("time"),
               between=c("gender", "age_t2", "days", "gse_t0"),
               observed=c("gender", "age_t2", "days", "gse_t0"),
               factorize=F) %>% apa::anova_apa(force_sph_corr=T)
#data.matched.analysis.long.scaled %>% summarise(missings.gse = mean(gse_t0 %>% is.na()))
#data.matched.analysis %>% filter(gse_t0 %>% is.na() == F) %>% nrow()
data.matched.analysis.long = data.matched.analysis.long %>% mutate(gse_group = ifelse(gse_t0 > median(gse_t0, na.rm=T), "high", "low") %>% factor(levels=c("low", "high")))
plot.gse = data.matched.analysis.long %>% filter(gse_group %>% is.na() == F) %>% 
  group_by(time, gse_group) %>% summarise(strain.se = se(strain), strain = mean(strain)) %>% 
  ggplot(aes(x=time, y=strain, color=gse_group, group=gse_group)) + 
  geom_hline(yintercept=0, linetype="dashed") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1, data=data.matched.analysis.long, aes(group=interaction(time, gse_group)), position=position_dodge(width=.9)) +
  #geom_violin(aes(group=interaction(time, gse_group)), data=data.matched.analysis.long) +
  geom_violin(aes(group=interaction(time, gse_group), fill=gse_group), color="black", alpha=.5, data=data.matched.analysis.long %>% filter(gse_group %>% is.na() == F)) +
  geom_boxplot(width=0.1, position=position_dodge(width=.9), data=data.matched.analysis.long %>% filter(gse_group %>% is.na() == F), aes(group=interaction(time, gse_group))) +
  geom_errorbar(aes(ymin=strain-strain.se*1.96, ymax=strain+strain.se*1.96), size=2, width=.25) +
  geom_line(size=2) + geom_point(size=5) + myGgTheme + ggtitle("GSE") +
  scale_x_discrete(labels = c("Pre","Peak","Downturn")) + scale_color_discrete("GSE") + scale_fill_discrete("GSE") + ylab("Strain") + xlab("")
plot.gse
#ggsave("3.4 GSE.png" %>% paste0("C:/Users/mar84qk/Dropbox/Arbeit/Corona Psychosocial/Ergebnisse/", .), width=1920, height=1080, units="px")
data.matched.analysis %>% select(gse_t0, strain_t0, strain_t1, strain_t2) %>% corrr::correlate()
with(data.matched.analysis %>% mutate(slope = strain_t1-strain_t0), cor.test(gse_t0, slope) %>% apa::cor_apa(r_ci=T))
with(data.matched.analysis %>% mutate(slope = strain_t2-strain_t1), cor.test(gse_t0, slope) %>% apa::cor_apa(r_ci=T))



# Cognitive Emotion Regulation Questionnaire (CERQ)
#CERQ maladaptive
data.matched.analysis.long.scaled %>% 
  afex::aov_ez(id="REF", dv="strain", .,
               within=c("time"),
               between=c("gender", "age_t2", "days", "cerq_maladapt_t0"),
               observed=c("gender", "age_t2", "days", "cerq_maladapt_t0"),
               factorize=F) %>% apa::anova_apa(force_sph_corr=T)
data.matched.analysis.long = data.matched.analysis.long %>% mutate(cerq_maladapt_group = ifelse(cerq_maladapt_t0 > median(cerq_maladapt_t0), "high", "low") %>% factor(levels=c("high", "low")))
plot.mal = data.matched.analysis.long %>% 
  group_by(time, cerq_maladapt_group) %>% summarise(strain.se = se(strain), strain = mean(strain)) %>% 
  ggplot(aes(x=time, y=strain, color=cerq_maladapt_group, group=cerq_maladapt_group)) + 
  geom_hline(yintercept=0, linetype="dashed") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1, data=data.matched.analysis.long, aes(group=interaction(time, cerq_maladapt)), position=position_dodge(width=.9)) +
  #geom_violin(aes(group=interaction(time, cerq_maladapt_group)), data=data.matched.analysis.long) +
  geom_violin(aes(group=interaction(time, cerq_maladapt_group), fill=cerq_maladapt_group), color="black", alpha=.5, data=data.matched.analysis.long) +
  geom_boxplot(width=0.1, position=position_dodge(width=.9), data=data.matched.analysis.long, aes(group=interaction(time, cerq_maladapt_group))) +
  geom_errorbar(aes(ymin=strain-strain.se*1.96, ymax=strain+strain.se*1.96), size=2, width=.25) +
  geom_line(size=2) + geom_point(size=5) + myGgTheme + ggtitle("CERQ-mal") +
  scale_x_discrete(labels = c("Pre","Peak","Downturn")) + scale_color_discrete("CERQ-mal") + scale_fill_discrete("CERQ-mal") + ylab("Strain") + xlab("")
plot.mal
#ggsave("3.1 CERQ mal.png" %>% paste0("C:/Users/mar84qk/Dropbox/Arbeit/Corona Psychosocial/Ergebnisse/", .), width=1920, height=1080, units="px")
data.matched.analysis %>% select(cerq_maladapt_t0, strain_t0, strain_t1, strain_t2) %>% corrr::correlate()
with(data.matched.analysis %>% mutate(slope = strain_t1-strain_t0), cor.test(cerq_maladapt_t0, slope) %>% apa::cor_apa(r_ci=T))
with(data.matched.analysis %>% mutate(slope = strain_t2-strain_t1), cor.test(cerq_maladapt_t0, slope) %>% apa::cor_apa(r_ci=T))

cowplot::plot_grid(plot.spai, plot.lsas, plot.gse, plot.mal, ncol=2, labels="AUTO") #Figure 3
#ggsave("Figure 3. Risk Factors.png" %>% paste0("C:/Users/mar84qk/Dropbox/Arbeit/Corona Psychosocial/Ergebnisse/", .), width=1920*2, height=1080*2, units="px")

#CERQ adaptive
data.matched.analysis.long.scaled %>% 
  afex::aov_ez(id="REF", dv="strain", .,
               within=c("time"),
               between=c("gender", "age_t2", "days", "cerq_adapt_t0"),
               observed=c("gender", "age_t2", "days", "cerq_adapt_t0"),
               factorize=F) %>% apa::anova_apa(force_sph_corr=T)
data.matched.analysis %>% select(cerq_adapt_t0, strain_t0, strain_t1, strain_t2) %>% corrr::correlate()

data.matched.analysis.long = data.matched.analysis.long %>% mutate(cerq_adapt_group = ifelse(cerq_adapt_t0 > median(cerq_adapt_t0), "high", "low") %>% factor(levels=c("low", "high")))
plot.adapt.simple = data.matched.analysis.long %>% 
  group_by(time, cerq_adapt_group) %>% summarise(strain.se = se(strain), strain = mean(strain)) %>% 
  ggplot(aes(x=time, y=strain, color=cerq_adapt_group, group=cerq_adapt_group)) + 
  geom_hline(yintercept=0, linetype="dashed") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1, data=data.matched.analysis.long, aes(group=interaction(time, cerq_adapt_group)), position=position_dodge(width=.9)) +
  #geom_violin(aes(group=interaction(time, cerq_adapt_group)), data=data.matched.analysis.long) +
  geom_violin(aes(group=interaction(time, cerq_adapt_group), fill=cerq_adapt_group), color="black", alpha=.5, data=data.matched.analysis.long) +
  geom_boxplot(width=0.1, position=position_dodge(width=.9), data=data.matched.analysis.long, aes(group=interaction(time, cerq_adapt_group))) +
  geom_errorbar(aes(ymin=strain-strain.se*1.96, ymax=strain+strain.se*1.96), size=2, width=.25) +
  geom_line(size=2) + geom_point(size=5) + myGgTheme + ggtitle("CERQ-adapt") +
  scale_x_discrete(labels = c("Pre","Peak","Downturn")) + scale_color_discrete("CERQ-adapt") + scale_fill_discrete("CERQ-adapt") + ylab("Strain") + xlab("")
plot.adapt.simple
#ggsave("3.2 CERQ adapt.png" %>% paste0("C:/Users/mar84qk/Dropbox/Arbeit/Corona Psychosocial/Ergebnisse/", .), width=1920, height=1080, units="px")

#interactions
data.matched.analysis.long = data.matched.analysis.long %>% mutate(days_group = ifelse(days > median(days), "lag high", "lag low") %>% as_factor())
plot.adapt = data.matched.analysis.long %>% 
  group_by(time, cerq_adapt_group, days_group, gender) %>% summarise(strain.se = se(strain), strain = mean(strain)) %>% 
  ggplot(aes(x=time, y=strain, color=cerq_adapt_group, group=cerq_adapt_group)) + #, shape=days_group)) + 
  #facet_wrap(vars(gender), ncol=1) +
  facet_grid(rows=vars(days_group), cols=vars(gender)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_violin(aes(group=interaction(time, cerq_adapt_group), fill=cerq_adapt_group), color="black", alpha=.5, data=data.matched.analysis.long) +
  geom_boxplot(width=0.1, position=position_dodge(width=.9), data=data.matched.analysis.long, aes(group=interaction(time, cerq_adapt_group))) +
  geom_errorbar(aes(ymin=strain-strain.se*1.96, ymax=strain+strain.se*1.96), size=2, width=.25) +
  geom_line(size=2) + geom_point(size=5) + myGgTheme + ggtitle("CERQ-adapt") +
  scale_x_discrete(labels = c("Pre","Peak","Downturn")) + scale_color_discrete("CERQ-adapt") + scale_fill_discrete("CERQ-adapt") + ylab("Strain") + xlab("")
plot.adapt
#ggsave("3.2.2 CERQ adapt gender days.png" %>% paste0("C:/Users/mar84qk/Dropbox/Arbeit/Corona Psychosocial/Ergebnisse/", .), width=1920*2, height=1080*2, units="px")

#male
data.matched.analysis.long.scaled %>% filter(gender=="male") %>% mutate(across(.cols=c("age_t2", "days", "cerq_adapt_t0"), scale)) %>% 
  afex::aov_ez(id="REF", dv="strain", .,
               within=c("time"),
               between=c("age_t2", "days", "cerq_adapt_t0"),
               observed=c("age_t2", "days", "cerq_adapt_t0"),
               factorize=F) %>% apa::anova_apa(force_sph_corr=T)
data.matched.analysis.long %>% mutate(across(.cols=c("age_t2", "days", "cerq_adapt_t0"), scale),
                                      cerqXdays = cerq_adapt_t0 * days) %>% 
  group_by(gender, time) %>% summarise(cor = cor.test(strain, cerqXdays) %>% apa::cor_apa(r_ci=T, print=F))

#female
data.matched.analysis.long.scaled %>% filter(gender=="female") %>% mutate(across(.cols=c("age_t2", "days", "cerq_adapt_t0"), scale)) %>% 
  afex::aov_ez(id="REF", dv="strain", .,
               within=c("time"),
               between=c("age_t2", "days", "cerq_adapt_t0"),
               observed=c("age_t2", "days", "cerq_adapt_t0"),
               factorize=F) %>% apa::anova_apa(force_sph_corr=T)
data.matched.analysis.long.scaled %>% filter(gender=="female") %>% mutate(across(.cols=c("age_t2", "days", "cerq_adapt_t0"), scale)) %>% 
  mutate(cerq_adapt_group = ifelse(cerq_adapt_t0 > median(cerq_adapt_t0), "high", "low") %>% as_factor(),
         age_group = ifelse(age_t2 > median(age_t2), "age high", "age low") %>% as_factor()) %>% 
  group_by(time, cerq_adapt_group, age_group) %>% summarise(strain.se = se(strain), strain = mean(strain)) %>% 
  ggplot(aes(x=time, y=strain, color=cerq_adapt_group, group=cerq_adapt_group)) +
  facet_wrap(vars(age_group), ncol=1) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_errorbar(aes(ymin=strain-strain.se, ymax=strain+strain.se), size=2) +
  geom_line(size=2) + geom_point(size=5) + myGgTheme +
  scale_x_discrete(labels = c("Pre","Peak","Downturn")) + scale_color_discrete("CERQ-adapt") + ylab("Strain") + xlab("")
#cerq_adapt low & age low = risk group for females (i.e. lower recovery)
#ggsave("3.2.3 CERQ adapt female age.png" %>% paste0("C:/Users/mar84qk/Dropbox/Arbeit/Corona Psychosocial/Ergebnisse/", .), width=1920, height=1080*2, units="px")

with(data.matched.analysis, cor.test(days, age_t2) %>% apa::cor_apa(r_ci=T))
data.matched.analysis %>% group_by(gender) %>% summarise(cor = cor.test(days, age_t2) %>% apa::cor_apa(r_ci=T, print=F))


#CERQ acceptance
data.matched.analysis.long.scaled %>% 
  afex::aov_ez(id="REF", dv="strain", .,
               within=c("time"),
               between=c("gender", "age_t2", "days", "cerq_acceptance_t0"),
               observed=c("gender", "age_t2", "days", "cerq_acceptance_t0"),
               factorize=F) %>% apa::anova_apa(force_sph_corr=T)
data.matched.analysis.long = data.matched.analysis.long %>% mutate(cerq_accept_group = ifelse(cerq_acceptance_t0 > median(cerq_acceptance_t0), "high", "low") %>% factor(levels=c("low", "high")))
plot.accept = data.matched.analysis.long %>% 
  group_by(time, cerq_accept_group, gender) %>% summarise(strain.se = se(strain), strain = mean(strain)) %>% 
  mutate(cerq_accept_group = cerq_accept_group %>% factor(levels=c("low", "high"))) %>% 
  ggplot(aes(x=time, y=strain, color=cerq_accept_group, group=cerq_accept_group)) + 
  facet_wrap(vars(gender), ncol=2) +
  geom_hline(yintercept=0, linetype="dashed") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1, data=data.matched.analysis.long, aes(group=interaction(time, cerq_accept_group)), position=position_dodge(width=.9)) +
  #geom_violin(aes(group=interaction(time, cerq_accept_group)), data=data.matched.analysis.long) +
  geom_violin(aes(group=interaction(time, cerq_accept_group), fill=cerq_accept_group), color="black", alpha=.5, data=data.matched.analysis.long) +
  geom_boxplot(width=0.1, position=position_dodge(width=.9), data=data.matched.analysis.long, aes(group=interaction(time, cerq_accept_group))) +
  geom_errorbar(aes(ymin=strain-strain.se*1.96, ymax=strain+strain.se*1.96), size=2, width=.25) +
  geom_line(size=2) + geom_point(size=5) + myGgTheme + ggtitle("CERQ-accept") +
  scale_x_discrete(labels = c("Pre","Peak","Downturn")) + scale_color_discrete("CERQ-accept") + scale_fill_discrete("CERQ-accept") + ylab("Strain") + xlab("")
plot.accept
#ggsave("3.3 CERQ accept.png" %>% paste0("C:/Users/mar84qk/Dropbox/Arbeit/Corona Psychosocial/Ergebnisse/", .), width=1920*2, height=1080, units="px")
data.matched.analysis %>% select(cerq_acceptance_t0, strain_t0, strain_t1, strain_t2) %>% corrr::correlate()

cowplot::plot_grid(plot.adapt, plot.accept, rel_heights = c(2, 1), ncol=1, labels="AUTO") #Figure 4
#ggsave("Figure 4. Risk Factors Gender.png" %>% paste0("C:/Users/mar84qk/Dropbox/Arbeit/Corona Psychosocial/Ergebnisse/", .), width=1920*2, height=1080*3, units="px")


#Childhood Trauma Questionnaire (CTQ)
data.matched.analysis.long.scaled %>% 
  afex::aov_ez(id="REF", dv="strain", .,
               within=c("time"),
               between=c("gender", "age_t2", "days", "ctq_t0"),
               observed=c("gender", "age_t2", "days", "ctq_t0"),
               factorize=F) %>% apa::anova_apa(force_sph_corr=T)
data.matched.analysis.long = data.matched.analysis.long %>% mutate(ctq_group = ifelse(ctq_t0 > median(ctq_t0), "high", "low") %>% factor(levels=c("high", "low")))
data.matched.analysis.long %>% 
  group_by(time, ctq_group) %>% summarise(strain.se = se(strain), strain = mean(strain)) %>% 
  mutate(ctq_group = ctq_group %>% factor(levels=c("high", "low"))) %>% 
  ggplot(aes(x=time, y=strain, color=ctq_group, group=ctq_group)) + 
  geom_hline(yintercept=0, linetype="dashed") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1, data=data.matched.analysis.long, aes(group=interaction(time, ctq_group)), position=position_dodge(width=.9)) +
  #geom_violin(aes(group=interaction(time, ctq_group)), data=data.matched.analysis.long) +
  geom_violin(aes(group=interaction(time, ctq_group), fill=ctq_group), color="black", alpha=.5, data=data.matched.analysis.long) +
  geom_boxplot(width=0.1, position=position_dodge(width=.9), data=data.matched.analysis.long, aes(group=interaction(time, ctq_group))) +
  geom_errorbar(aes(ymin=strain-strain.se*1.96, ymax=strain+strain.se*1.96), size=2, width=.25) +
  geom_line(size=2) + geom_point(size=5) + myGgTheme +
  scale_x_discrete(labels = c("Pre","Peak","Downturn")) + scale_color_discrete("CTQ") + scale_fill_discrete("CTQ") + ylab("Strain") + xlab("")
#ggsave("4.1 CTQ.png" %>% paste0("C:/Users/mar84qk/Dropbox/Arbeit/Corona Psychosocial/Ergebnisse/", .), width=1920, height=1080, units="px")
#data.matched.analysis %>% select(ctq_t0, strain_t0, strain_t1, strain_t2) %>% corrr::correlate()
with(data.matched.analysis, cor.test(ctq_t0, strain_t0) %>% apa::cor_apa(r_ci=T))
with(data.matched.analysis, cor.test(ctq_t0, strain_t1) %>% apa::cor_apa(r_ci=T))
with(data.matched.analysis, cor.test(ctq_t0, strain_t2) %>% apa::cor_apa(r_ci=T))

with(data.matched.analysis %>% mutate(slope = strain_t1-strain_t0), cor.test(ctq_t0, slope) %>% apa::cor_apa(r_ci=T))
with(data.matched.analysis %>% mutate(slope = strain_t2-strain_t1), cor.test(ctq_t0, slope) %>% apa::cor_apa(r_ci=T))


#Life Threatening Events (LTE)
data.matched.analysis.long.scaled %>% 
  afex::aov_ez(id="REF", dv="strain", .,
               within=c("time"),
               between=c("gender", "age_t2", "days", "lte_t0"),
               observed=c("gender", "age_t2", "days", "lte_t0"),
               factorize=F) %>% apa::anova_apa(force_sph_corr=T)
data.matched.analysis.long = data.matched.analysis.long %>% mutate(lte_group = ifelse(lte_t0 > median(lte_t0, na.rm=T), "high", "low") %>% factor(levels=c("high", "low")))
data.matched.analysis.long %>% filter(lte_t0 %>% is.na() == F) %>% 
  group_by(time, lte_group) %>% summarise(strain.se = se(strain), strain = mean(strain)) %>% 
  ggplot(aes(x=time, y=strain, color=lte_group, group=lte_group)) + 
  geom_hline(yintercept=0, linetype="dashed") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1, data=data.matched.analysis.long, aes(group=interaction(time, lte_group)), position=position_dodge(width=.9)) +
  #geom_violin(aes(group=interaction(time, lte_group)), data=data.matched.analysis.long) +
  geom_violin(aes(group=interaction(time, lte_group), fill=lte_group), color="black", alpha=.5, data=data.matched.analysis.long %>% filter(lte_t0 %>% is.na() == F)) +
  geom_boxplot(width=0.1, position=position_dodge(width=.9), data=data.matched.analysis.long %>% filter(lte_t0 %>% is.na() == F), aes(group=interaction(time, lte_group))) +
  geom_errorbar(aes(ymin=strain-strain.se*1.96, ymax=strain+strain.se*1.96), size=2, width=.25) +
  geom_line(size=2) + geom_point(size=5) + myGgTheme +
  scale_x_discrete(labels = c("Pre","Peak","Downturn")) + scale_color_discrete("LTE") + scale_fill_discrete("LTE") + ylab("Strain") + xlab("")
#ggsave("4.2 LTE.png" %>% paste0("C:/Users/mar84qk/Dropbox/Arbeit/Corona Psychosocial/Ergebnisse/", .), width=1920, height=1080, units="px")
data.matched.analysis %>% select(lte_t0, strain_t0, strain_t1, strain_t2) %>% corrr::correlate()

#Adverse Life Events ([A]LE)
data.matched.analysis.long.scaled %>% 
  afex::aov_ez(id="REF", dv="strain", .,
               within=c("time"),
               between=c("gender", "age_t2", "days", "le_t0"),
               observed=c("gender", "age_t2", "days", "le_t0"),
               factorize=F) %>% apa::anova_apa(force_sph_corr=T)
data.matched.analysis.long = data.matched.analysis.long %>% mutate(le_group = ifelse(le_t0 > median(le_t0), "high", "low") %>% factor(levels=c("high", "low")))
data.matched.analysis.long %>% filter(le_t0 %>% is.na() == F) %>% 
  mutate(le_group = ifelse(le_t0 > median(le_t0), "high", "low") %>% as_factor()) %>% 
  group_by(time, le_group) %>% summarise(strain.se = se(strain), strain = mean(strain)) %>% 
  ggplot(aes(x=time, y=strain, color=le_group, group=le_group)) + 
  geom_hline(yintercept=0, linetype="dashed") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=1, data=data.matched.analysis.long, aes(group=interaction(time, le_group)), position=position_dodge(width=.9)) +
  #geom_violin(aes(group=interaction(time, le_group)), data=data.matched.analysis.long) +
  geom_violin(aes(group=interaction(time, le_group), fill=le_group), color="black", alpha=.5, data=data.matched.analysis.long %>% filter(le_t0 %>% is.na() == F)) +
  geom_boxplot(width=0.1, position=position_dodge(width=.9), data=data.matched.analysis.long %>% filter(le_t0 %>% is.na() == F), aes(group=interaction(time, le_group))) +
  geom_errorbar(aes(ymin=strain-strain.se*1.96, ymax=strain+strain.se*1.96), size=2, width=.25) +
  geom_line(size=2) + geom_point(size=5) + myGgTheme +
  scale_x_discrete(labels = c("Pre","Peak","Downturn")) + scale_color_discrete("ALE") + scale_fill_discrete("ALE") + ylab("Strain") + xlab("")
#ggsave("4.3 LE.png" %>% paste0("C:/Users/mar84qk/Dropbox/Arbeit/Corona Psychosocial/Ergebnisse/", .), width=1920, height=1080, units="px")
data.matched.analysis %>% select(le_t0, strain_t0, strain_t1, strain_t2) %>% corrr::correlate()

#5x interaction
# data.le %>% group_by(time, gender, age_group) %>% mutate(n_group = n()) %>% 
#   group_by(time, gender, le_group, age_group, lag_group) %>% 
#   summarise(strain.se = se(strain), strain = mean(strain), n = n(), gsize = n / n_group) %>% 
#   ggplot(aes(x=time, y=strain, color=le_group, shape=lag_group, group=interaction(le_group, lag_group))) + 
#   facet_grid(rows=vars(age_group), cols=vars(gender)) +
#   geom_hline(yintercept=0, linetype="dashed") +
#   geom_errorbar(aes(ymin=strain-strain.se, ymax=strain+strain.se), size=2, position=position_dodge(width=.5)) +
#   #geom_line(size=2, position=position_dodge(width=.5)) + geom_point(size=5, position=position_dodge(width=.5)) + myGgTheme
#   geom_line(size=2, position=position_dodge(width=.5)) + geom_point(aes(size=gsize), position=position_dodge(width=.5)) + myGgTheme

#deconstruction of effects
data.le = data.matched.analysis.long.scaled %>% 
  mutate(le_group = ifelse(le_t0 > median(le_t0), "high", "low") %>% as_factor(),
         age_group = ifelse(age_t2 > median(age_t2), "age high", "age low"),
         lag_group = ifelse(days > median(days), "lag high", "lag low"))

for (a in data.le %>% pull(age_group) %>% unique()) {
  cat(a, ":\n")
  data.le %>% filter(age_group == a) %>% 
    afex::aov_ez(id="REF", dv="strain", .,
                 within=c("time"),
                 between=c("gender", "days", "le_t0"),
                 observed=c("gender", "days", "le_t0"),
                 factorize=F) %>% apa::anova_apa(force_sph_corr=T)
} #below median: few effects => only look at above median
a = data.le %>% filter(age_group %>% grepl("high", .)) %>% pull(age_group) %>% unique()
for (g in data.le %>% pull(gender) %>% unique()) {
  cat(g, " (& ", a, "):\n")
  data.le %>% filter(age_group == a, gender == g) %>% 
    afex::aov_ez(id="REF", dv="strain", .,
                 within=c("time"),
                 between=c("days", "le_t0"),
                 observed=c("days", "le_t0"),
                 factorize=F) %>% apa::anova_apa(force_sph_corr=T)
}
g = "male"
data.le %>% group_by(time, gender, age_group) %>% mutate(n_group = n()) %>%
  group_by(time, gender, le_group, age_group, lag_group) %>%
  summarise(strain.se = se(strain), strain = mean(strain), n = n(), gsize = n / n_group) %>%
  filter(gender == g, age_group == a) %>% 
  ggplot(aes(x=time, y=strain, color=le_group, group=le_group)) +
  facet_wrap(vars(lag_group)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_errorbar(aes(ymin=strain-strain.se, ymax=strain+strain.se), size=2, position=position_dodge(width=.5)) +
  geom_line(size=2, position=position_dodge(width=.5)) + geom_point(size=5, position=position_dodge(width=.5)) + myGgTheme
  #geom_line(size=2, position=position_dodge(width=.5)) + geom_point(aes(size=gsize), position=position_dodge(width=.5)) + myGgTheme
for (l in data.le %>% pull(lag_group) %>% unique()) {
  cat(l, " (& ", a, " & ", g, "):\n")
  with(data.le %>% filter(gender == g, age_group == a, lag_group == l) %>% group_by(REF) %>% summarise(strain = mean(strain), le_t0 = mean(le_t0)), 
       cor.test(strain, le_t0)) %>% apa::cor_apa()
}


#PHQ2 (T2! => hard to interpret)
# data.matched.analysis.long.scaled %>% 
#   afex::aov_ez(id="REF", dv="strain", .,
#                within=c("time"),
#                between=c("gender", "age_t2", "days", "phq2_t2"),
#                observed=c("gender", "age_t2", "days", "phq2_t2"),
#                factorize=F) %>% apa::anova_apa(force_sph_corr=T)
# data.matched.analysis.long %>% filter(phq2_t2 %>% is.na() == F) %>% 
#   mutate(phq_group = ifelse(phq2_t2 > median(phq2_t2), "high", "low") %>% as_factor()) %>% 
#   group_by(time, phq_group) %>% summarise(strain.se = se(strain), strain = mean(strain)) %>% 
#   ggplot(aes(x=time, y=strain, color=phq_group)) + 
#   geom_hline(yintercept=0, linetype="dashed") +
#   geom_errorbar(aes(ymin=strain-strain.se, ymax=strain+strain.se), size=2) +
#   geom_line(size=2) + geom_point(size=5) + myGgTheme
# data.matched.analysis %>% select(phq2_t2, strain_t0, strain_t1, strain_t2) %>% corrr::correlate()


#Interactions?
data.matched.analysis %>% select(spai_t0, cerq_maladapt_t0, gse_t0, ctq_t0, phq2_t2) %>% corrr::correlate()
data.matched.analysis.long.scaled %>% 
  afex::aov_ez(id="REF", dv="strain", .,
               within=c("time"),
               # between=c("gender", "age_t2", "days", "spai_t0", "cerq_maladapt_t0", "gse_t0", "ctq_t0", "phq2_t2"),
               # observed=c("gender", "age_t2", "days", "spai_t0", "cerq_maladapt_t0", "gse_t0", "ctq_t0", "phq2_t2"),
               # between=c("gender", "age_t2", "spai_t0", "cerq_maladapt_t0", "gse_t0", "ctq_t0", "phq2_t2"),
               # observed=c("gender", "age_t2", "spai_t0", "cerq_maladapt_t0", "gse_t0", "ctq_t0", "phq2_t2"),
               between=c("gender", "age_t2", "spai_t0", "cerq_maladapt_t0", "gse_t0", "ctq_t0"),
               observed=c("gender", "age_t2", "spai_t0", "cerq_maladapt_t0", "gse_t0", "ctq_t0"),
               factorize=F) %>% apa::anova_apa(force_sph_corr=T)
