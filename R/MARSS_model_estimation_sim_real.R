library(tidyverse)
library(ggplot2)
library(ggforce)
library(atsar)
library(rstan)
library(loo)
library(magrittr)
# remotes::install_github("twangss/atsa-seals")
library(atsar)

source("R/marss_stan_functions.R")

pinn_df=read.csv("data/Pinniped database for publication.csv")
turt_df_long <- pinn_df %>% arrange(ID) %>%
  filter(Use == "yes") %>%
  mutate(CommonName = Common.Name, Site = Location.of.population, RMU = Subspecies) %>%
  mutate(ts_id = paste(CommonName,"#",Site)) %>% 
  pivot_longer(
    cols = starts_with("y"),
    names_to = "Year",
    names_prefix = "y",
    values_to = "Count",
    values_drop_na = TRUE) %>%
  mutate(Count=as.numeric(Count),Year=as.numeric(Year)) %>%
  group_by(ID) %>%
  drop_na(Count) %>%
  mutate(n=n()) %>% filter(n>3) # MUST BE THREE 


start_date = Sys.Date()
N_sims = 10 # if modeling simulated data
true_state="SIM" # if modeling simulated data
true_state="real" # if modeling real data
sim_i=0
for (true_state in c("Ueq","GP")){
for (sim_i in 1:N_sims){

species_names <- unique(turt_df_long$CommonName)
N_species <- length(species_names)

# get the distance matrix
distance_matrix_list = readRDS("data/site_distance/seal_distance_matrix_subsp_2024_1218.RDS")

Obs_method = T # Obs param est by survey methods toggle

mod_sel_grid=NULL
error_params=NULL
pred_data=NULL
loo_all = NULL
gaussian_rel = NULL
# sim_i=1
taxa = paste("seal_real",true_state,sim_i,sep="")
model_type_name = "GPcorr_prior"
mod_sel_file_name = paste("data/simulation_GP_MA1/sim_outputs/",taxa,"_mod_sel_grid_",model_type_name,"_",start_date,".RDS",sep="")
pred_data_file_name = paste("data/simulation_GP_MA1/sim_outputs/",taxa,"_pred_data_",model_type_name,"_",start_date,".RDS",sep="")
error_params_file_name = paste("data/simulation_GP_MA1/sim_outputs/",taxa,"_error_params_",model_type_name,"_",start_date,".RDS",sep="")
loo_file_name = paste("data/simulation_GP_MA1/sim_outputs/",
                      taxa,"_loo_",model_type_name,"_",start_date,".RDS",sep="")
gaussian_file_name = paste("data/simulation_GP_MA1/sim_outputs/",
                      taxa,"_gp_",model_type_name,"_",start_date,".RDS",sep="")

options(mc.cores = parallel::detectCores())

for (i in 1:N_species){
  species_i <- species_names[i]
  cat(species_i,"\n")
  species_i_pre_df <- turt_df_long %>% 
    filter(CommonName == species_i) %>% arrange(ID) %>%
    mutate(log.spawner = log(Count+1)) %>% 
    group_by(ID) %>% 
    mutate(n=sum(!is.na(log.spawner))) 
  
  species_i_df <- species_i_pre_df %>% ungroup() %>% 
    dplyr::select(ID,Year,Site, log.spawner) %>%  # get just the columns that I need
    # drop_na(log.spawner) %>%
    group_by(ID,Site) %>% 
    complete(Year = min(Year):max(Year),
             fill = list(log.spawner = NA)) %>%
    ungroup() %>%
    pivot_wider(names_from = c("Site","ID"), values_from = "log.spawner") %>% 
    arrange(Year) %>% 
    column_to_rownames(var = "Year") %>% # make the years rownames
    # filter(if_any(everything(), ~ !is.na(.)))%>%
    as.matrix() %>% # turn into a matrix with year down the rows
    t() # make time across the columns
  
  # seals
  # id_vec for seals which has multiple mapping onto 1
  id_vec <- species_i_pre_df %>% subset(CommonName == species_i) %>% ungroup() %>% 
    # filter(Indiv.sum==1) %>% # need this for correct mapping right? but you commented it out...maybe for sims? 
    distinct(ID,Indiv.pop.map) %>% pull(Indiv.pop.map)
  
  add_up_vec <- species_i_pre_df %>% subset(CommonName == species_i) %>% ungroup() %>% 
    # filter(Indiv.sum==1) %>% # need this for correct mapping right? but you commented it out...maybe for sims? 
    distinct(ID,Indiv.sum,Scaling) %>% mutate(add_up = Indiv.sum*Scaling) %>% pull(add_up)
  
  RMU_vec = species_i_pre_df %>% subset(CommonName == species_i) %>% ungroup() %>% filter(Indiv.sum==1) %>% distinct(ID,RMU) %>%
    pull(RMU)

  #== Set up for Stan
  # seal
  n_states <- max(as.numeric(fct_inorder(id_vec))) # site_vec
  id_states = as.numeric(fct_inorder(id_vec))
  
  # all
  r_unequal <- seq(1, nrow(species_i_df))
  r_equal <- rep(1, nrow(species_i_df))
  uq_unequal <- seq(1, n_states)
  uq_equal <- rep(1, n_states)
  stocks <- (as.numeric(factor(RMU_vec))) # site_vec
  
      cat("Species",species_i,"\n")

      # set up mcmc options
      mcmc_list  = list(n_mcmc = 4000, n_burn = 1000, n_chain = 3, n_thin = 1,step_size=0.4,adapt_delta=0.9)
      
      if (Obs_method == T){
        method_vec = species_i_pre_df %>% subset(CommonName == species_i) %>% 
          distinct(ID,Sampling.method) %>%
          pull(Sampling.method)
        method_vec_in = as.numeric(fct_inorder(method_vec))
      } else { method_vec_in = r_equal
      method_vec = "sigma_obs"}
      
      #= MARSS specifications
      marss_list = list(states = id_states, #site_vec
                        obsVariances = method_vec_in, #as.numeric(factor(units_vec)), OR r_equal
                        proVariances = uq_equal,
                        trends = uq_equal,
                        stocks = uq_unequal) # stocks = stocks, this was 
      
      # set up data_list
      data_list_tmp <- setup_data(y = species_i_df,
                                  est_nu = FALSE,
                                  est_trend = FALSE,
                                  family = "gaussian", 
                                  # mcmc_list = mcmc_list,
                                  marss = marss_list)
      
      
      #== Spatial covariance matrix
      data_list_tmp$data$Dist = distance_matrix_list[[species_i]]
      data_list_tmp$data$add_up = add_up_vec
      
      if (nrow(species_i_df) == 1) {data_list_tmp$data$Dist = matrix(0,1,1) }
      if (species_i == "New Zealand sea lion") {
        
        data_list_tmp$data$Dist = matrix(1,4,4)  
      }
      
      model_vec = c("Ueq","GP")

    for (model_i in model_vec){
      cat("true_state",true_state,"\n")
      cat("Model",model_i,"\n")
      model_type = model_i
      
      if (model_type == "Ueq") {cmd_file_name = "marss_cmd.stan"}
      if (model_type == "GP") {cmd_file_name = "marss_gp_cmd_old.stan"}

      
        # Priors, these will only be useful for cmd stan files with "_prior" in the name, for prior sensitivity
      # Define paired values for each parameter
      param1 <- list(c(0.05, 0.2), c(0.1, 0.1))
      param2 <- list(c(0.05, 0.05), c(0.01, 0.01))
      param3 <- list(c(0, 80000), c(0, 160000))
      param4 <- list(c(1, 2), c(1, 1))
      
      # Generate all combinations of the paired values
      prior_combo <- expand.grid(param1, param2, param3, param4)
      
      # Flatten the list columns into separate columns
      prior_combo <- data.frame(
        Param1_1 = sapply(prior_combo$Var1, "[", 1),
        Param1_2 = sapply(prior_combo$Var1, "[", 2),
        Param2_1 = sapply(prior_combo$Var2, "[", 1),
        Param2_2 = sapply(prior_combo$Var2, "[", 2),
        Param3_1 = sapply(prior_combo$Var3, "[", 1),
        Param3_2 = sapply(prior_combo$Var3, "[", 2),
        Param4_1 = sapply(prior_combo$Var4, "[", 1),
        Param4_2 = sapply(prior_combo$Var4, "[", 2)
        )
      
      for (prior_i in 1:nrow(prior_combo)){
        data_list_tmp$data$sigma_obs_priors = c(prior_combo[prior_i,1],prior_combo[prior_i,2])
        data_list_tmp$data$sigma_proc_priors = c(prior_combo[prior_i,3],prior_combo[prior_i,4])
        data_list_tmp$data$rho_gp_priors = c(prior_combo[prior_i,5],prior_combo[prior_i,6])
        data_list_tmp$data$alpha_gp_priors = c(prior_combo[prior_i,7],prior_combo[prior_i,8])
        
      
      library(cmdstanr)
      file <- file.path(cmdstan_path(), "pinnipeds", cmd_file_name)
      mod <- cmdstan_model(file)
      
      fit2 <- mod$sample(
        data = data_list_tmp$data,
        seed = 124,
        chains = mcmc_list$n_chain,
        parallel_chains = mcmc_list$n_chain,
        iter_warmup = mcmc_list$n_burn,
        iter_sampling = mcmc_list$n_mcmc,
        adapt_delta=0.99,
        step_size=0.4,
        refresh = 500 # print update every 500 iters
      )
      cat("fitting is done \n")
      # stan_est = extract(fit2)
      fit2_sum = fit2$summary()
      # save model diags
      diag_summ = fit2$diagnostic_summary()
      divergent_perc = sum(diag_summ$num_divergent)/(mcmc_list$n_chain*mcmc_list$n_mcmc)

      eff_sample_size = fit2_sum$ess_tail/(mcmc_list$n_chain*mcmc_list$n_mcmc)
      n_eff0.01_tmp = sum(eff_sample_size<0.01,na.rm=T)/length(eff_sample_size)
      
      rhat1.1_tmp = sum(fit2_sum$rhat>1.01,na.rm=T)/length(fit2_sum$rhat)
      
      loo_tmp = loo(fit2$draws("log_lik"))
      
      tot_sum_x = fit2$draws("tot_sum_x",format = "matrix")

      tot_sum = data.frame(tot.mean = apply(tot_sum_x, 2, mean),
                                 tot.sd = apply(tot_sum_x, 2, sd),
                                 tot.q05 = apply(tot_sum_x,2,quantile, probs = c(0.05)),
                                 tot.q10 = apply(tot_sum_x,2,quantile, probs = c(0.10)),
                                 tot.q25 = apply(tot_sum_x,2,quantile, probs = c(0.25)),
                                 tot.q50 = apply(tot_sum_x,2,quantile, probs = c(0.50)),
                                 tot.q75 = apply(tot_sum_x,2,quantile, probs = c(0.75)),
                                 tot.q90 = apply(tot_sum_x,2,quantile, probs = c(0.9)),
                                 tot.q95 = apply(tot_sum_x,2,quantile, probs = c(0.95)),
                                 var="total_abundance_norm") %>%
        mutate(year = as.numeric(colnames(species_i_df)))
      
      
      if (model_type == "GP") {
      gaussian_rel_out = fit2$draws("gaussian_rel",format = "matrix")
      gaussian_rel_tmp = data.frame(tot.q10 = apply(gaussian_rel_out,2,quantile, probs = c(0.10)),
                           tot.q50 = apply(gaussian_rel_out,2,quantile, probs = c(0.50)),
                           tot.q90 = apply(gaussian_rel_out,2,quantile, probs = c(0.90)),
                           var="gaussian_rel",species = species_i) %>%
        mutate(x = 1:10000)
      gaussian_rel = rbind(gaussian_rel,gaussian_rel_tmp) 
      saveRDS(gaussian_rel,gaussian_file_name)
      }
      
      loo_list_name = paste(species_i,model_type,sep="-")
      loo_all[[loo_list_name]] = loo_tmp
      
      saveRDS(loo_all,loo_file_name)
      
      model_diags = data.frame(species = species_i, 
                               rhat1.1=rhat1.1_tmp,n_eff0.01=n_eff0.01_tmp,div_p=divergent_perc,
                               loo = loo_tmp$estimates["looic","Estimate"],
                               loo_se = loo_tmp$estimates["looic","SE"]) %>% 
        mutate(model_type = model_type,prior = prior_i)
      mod_sel_grid = rbind(mod_sel_grid,model_diags) 

      # Save Process/Obs Sigmas
      if (model_type == "Ueq") {proc_variable_names = c("sigma_process_real")}
      if (model_type == "GP") {proc_variable_names = c("sigma_process_real","alpha_real","rho_gp")}

    
      sigma_proc_draws <- fit2$draws(proc_variable_names,format = "matrix")
      sigma_proc_tmp = sigma_proc_draws %>% data.frame() %>% pivot_longer(everything(), names_to = "units") %>%
        group_by(units) %>%
        summarise(mean = mean(value),
                  sd = sd(value),
                  q025 = quantile(value, probs = c(0.025)),
                  q50 = quantile(value, probs = c(0.5)),
                  q975 = quantile(value, probs = c(0.975)),
                  rhat = fit2$summary(proc_variable_names)$rhat) %>%
        mutate(var="sigma_proc")

      sigma_obs_draws <- fit2$draws("sigma_obs",format = "matrix")
      sigma_obs_tmp = data.frame(units = levels(fct_inorder(method_vec)), #"sigma_obs", #unique_units_vec, #levels(fct_inorder(method_vec))
                                 mean = apply(sigma_obs_draws, 2, mean),
                                 sd = apply(sigma_obs_draws, 2, sd),
                                 q025 = apply(sigma_obs_draws,2,quantile, probs = c(0.025)),
                                 q50 = quantile(sigma_obs_draws, probs = c(0.5)),
                                 q975 = apply(sigma_obs_draws,2,quantile, probs = c(0.975)),
                                 rhat = fit2$summary("sigma_obs")$rhat,
                                 var="sigma_obs")
    
      
      error_params_tmp = rbind(sigma_proc_tmp,sigma_obs_tmp) %>% mutate(model_type = model_type,prior = prior_i)
      error_params_tmp$species = species_i 
      error_params = rbind(error_params,error_params_tmp)
      
      # join raw data to plot 
      raw_data = species_i_df %>% data.frame() %>% 
        set_colnames(colnames(species_i_df)) %>%
        rownames_to_column("ts") %>% 
        pivot_longer(!ts,
                     names_to="year", 
                     values_to = "data") %>% 
        mutate(year=as.numeric(year)) 
      
      
      # extract PRED  from commandstanr
      preds = fit2$summary(variables = "pred")
      preds_matrix = cmdstanr_preds_sd(preds,species_i_df)
      pred_data_tmp = preds_matrix %>% left_join(raw_data) %>% mutate(species = species_i, model_type = model_type,prior = prior_i) %>%
        left_join(tot_sum)
      
      p = pred_data_tmp %>% ggplot() + 
        geom_ribbon(aes(x=year,ymin = mean-1.96*sd, ymax = mean+1.96*sd), alpha=0.2) +
        # geom_line(aes(x=year,y=state_fit)) +
        geom_line(aes(x=year,y=(mean))) +
        geom_point(aes(x=year,y=(data)),color="red")+
        facet_wrap(~ts,scales="free") + theme_bw()
      print(p)
      
      pred_data=rbind(pred_data,pred_data_tmp)
      
      saveRDS(mod_sel_grid,mod_sel_file_name)
      saveRDS(error_params,error_params_file_name)
      saveRDS(pred_data,pred_data_file_name)
      
      
      
  }
}
}
}}


# seals PRED
pred_data_file_name = "model/seal_realreal0_pred_data_GPcorr_2025-03-17.RDS"
pred_data = readRDS(pred_data_file_name) 


## Model selection
loo_tbl_pre = NULL
i=1
for (i in seq(1,length(loo_all),by=2)){ 
  speciescut  =  names(loo_all)[i]
  species_i = sub("\\-.*", "", speciescut)
  
  species_i = str_replace(speciescut, "-GP", "")
  species_i = str_replace(species_i, "-Ueq", "")
  
  loo_compare_tbl = loo_compare(loo_all[[i]],loo_all[[i+1]])
  loo_compare_tbl_sort = loo_compare_tbl[order((row.names(loo_compare_tbl))),]
  
  loo_upper = loo_compare_tbl[2,"elpd_diff"] + 1.645*loo_compare_tbl[2,"se_diff"] 
  loo_lower = loo_compare_tbl[2,"elpd_diff"] - 1.645*loo_compare_tbl[2,"se_diff"] 
  if (loo_upper < -4) {
    if(row.names(loo_compare_tbl)[1] == "model2"){best_model_i =  "GP"}
    if(row.names(loo_compare_tbl)[1] == "model1"){best_model_i =  "Ueq"}
  } else {best_model_i = c("GP","Ueq")}
  
  loo_tbl_i = data.frame(species = species_i,model_type = best_model_i)
  loo_tbl_pre = rbind(loo_tbl_pre, loo_tbl_i)
}

loo_tbl = loo_tbl_pre %>% mutate(keep = 1)  %>% 
  group_by(species) %>% 
  mutate(priority = case_when(model_type == "GP" ~ 2,
                              model_type == "Ueq" ~ 1,
                              .default = NA)) %>%
  filter(priority == min(priority))

mod_sel_grid = mod_sel_grid %>% 
  mutate(priority = case_when(model_type == "GP" ~ 2,
                                             model_type == "Ueq" ~ 1,
                                             .default = NA))

model_grid_best = mod_sel_grid %>% 
  mutate(priority = case_when(model_type == "GP" ~ 2,
                              model_type == "Ueq" ~ 1,
                              .default = NA)) %>%
  left_join(loo_tbl) %>%
  group_by(species) %>% 
  filter(keep == 1) %>%
  filter(priority == min(priority))

best_pass2_model_fits = pred_data %>% 
  inner_join(model_grid_best, by = c("species","model_type")) 


# Plot 
#=== fix the errors to be constant, find the distance (# of years) to the nearest data value
best_pass2_model_fits_dist = best_pass2_model_fits %>% 
  group_by(ts) %>%
  mutate(distance = NA)

ts_fits_distance = NULL
ts_vec_loop = unique(best_pass2_model_fits_dist$ts)
for (i in 1:length(ts_vec_loop)){
  ts_tmp = ts_vec_loop[i]
  ts_fits_tmp = best_pass2_model_fits_dist %>% filter(ts == ts_tmp)
  
  a <- which(!is.na(ts_fits_tmp$data)) 
  b <- which(is.na(ts_fits_tmp$data))
  
  for (i in 1:length(b)){
    dist_tmp = min(abs(a - b[i]),na.rm=T)
    ts_fits_tmp$distance[b[i]] = dist_tmp
  }
  
  ts_fits_distance = rbind(ts_fits_distance,ts_fits_tmp)
}

# fix the ballooning errors
best_pass2_model_fits_error_fix = ts_fits_distance %>% ungroup() %>% # 
  mutate(sd_tmp = sd) %>%
  mutate(sd_tmp = case_when(distance > 5 ~ NA,
                             .default = sd_tmp)) %>%
  dplyr::group_by(ts) %>%
  # this linearly interpolates NA values INBETWEEN data years, na.rm keeps the ends NA
  # mutate(lo_diff_approx = zoo::na.approx(lo_diff, na.rm=FALSE)) %>% 
  # mutate(hi_diff_approx = zoo::na.approx(hi_diff, na.rm=FALSE)) %>%
  # toggle below on if we want to extend the TS to all years all time
  # ungroup() %>%
  # complete(ts,year, # fill in missing years with NA
  #          fill = list(hi_diff_approx = NA)) %>%
  # group_by(ts) %>%
  # this fills in the ends
  # fill(lo_diff_approx) %>% 
  fill(c(sd_tmp), .direction = "downup") #change to up only if necessary


# add unit scaling to model fits
best_pass2_model_fits_error_fix_id = best_pass2_model_fits_error_fix %>%
  # mutate(ID = as.numeric(str_extract(ts, "[^_]+")))
  mutate(ID = as.numeric(stringr::word(ts, 2, sep = "_")))

scaling_tbl = turt_df_long %>% 
  ungroup() %>%
  dplyr::select(CommonName,ID,Scaling,RMU, Units,Indiv.sum) %>% # 
  distinct()

best_pass2_model_fits_error_fix_scale = best_pass2_model_fits_error_fix_id %>% 
  left_join(scaling_tbl) %>% 
  drop_na(RMU) # the dataset was edited between model_type runs so dropping previously included sites that are now excluded

# saveRDS(best_pass2_model_fits_error_fix_scale, "model/outputs/turtle_best_model_fits_error_fix_scale_20240716.RDS")
# best_pass2_model_fits_error_fix_scale = readRDS("model/outputs/turtle_best_model_fits_error_fix_scale_20240716.RDS")

best_pass2_model_fits_error_fix_scale_sd = best_pass2_model_fits_error_fix_scale %>% 
  ungroup() %>% 
  # group_by(year) %>%
  mutate(norm_mean = exp((mean)), # just look at mode, if mean then: +0.5*sd_tmp^2
         normal_sd = exp(mean+(sd_tmp^2)/2)*sqrt(exp(sd_tmp^2)-1))

#=== plot the Stan fits
best_pass2_model_fits_error_fix_scale$label <- paste(best_pass2_model_fits_error_fix_scale$ts,"\n","units = ", best_pass2_model_fits_error_fix_scale$Units,sep="")
best_pass2_species = unique(best_pass2_model_fits_error_fix_scale$species)
# pdf("model/plots/seals_best_species_marss_stan_ts_fits_20250313.pdf", width = 11, height = 8.5)
ncol = 5
nrow=4
for (species_i in best_pass2_species){
  print(species_i)
  # species_i = "Loggerhead"
  species_best_pass2_fits = best_pass2_model_fits_error_fix_scale %>% filter(species == species_i)
  # proc_state_vec = unique(species_best_pass2_fits$proc_state)
  RMU_vec = unique(species_best_pass2_fits$RMU)
  # RMU_i="North Pacific"
  for(RMU_i in RMU_vec){
  # for (proc_state_tmp in proc_state_vec){
    species_best_pass2_fits_proc = species_best_pass2_fits %>% 
      # filter(proc_state == proc_state_tmp) %>%
      filter(RMU == RMU_i)
  page_length = ceiling(length(unique(species_best_pass2_fits_proc$ts))/(ncol*nrow))
  for (page_i in 1:page_length){
    # page_i=1
    p = species_best_pass2_fits_proc %>% ggplot() +
      geom_ribbon(aes(x=year,ymin = mean-sd_tmp, ymax = mean+sd_tmp), alpha=0.2) +
      geom_line(aes(x=year,y=mean)) +
      geom_point(aes(x=year,y=data),color="red")+
      facet_wrap_paginate(~label,scales="free", ncol = ncol, nrow = nrow, page = page_i) +
      # facet_wrap(~ts,scales="free") + 
      theme_bw() +
      ggtitle(paste(species_i,RMU_i,
                    # species_best_pass2_fits_proc$proc_state[1],
                    sep=" - "))
    print(p)
  #}
  }}}
dev.off()


#==== INDEX BY SPECIES
# plot species trends
RMU_species_number_ts = turt_df_long %>% ungroup()%>%
  # dplyr::select(CommonName,ID,Units_clean) %>% 
  # distinct(CommonName,ID,Units_clean) %>%
  group_by(CommonName) %>% 
  mutate(species_n = n_distinct(ID)) %>% 
  filter(Indiv.sum == 1) %>%
  filter(Scaling==1) %>% # to get the common denominator untis
  dplyr::select(CommonName,species_n,Units_clean) %>% distinct()

# join n summaries to previous table
index_df_tmp_n_species = index_df_tmp_species %>% 
  left_join(RMU_species_number_ts, by = c("species" = "CommonName"))

# POSTERIOR PRED
index_df_tmp_n_species = best_pass2_model_fits_error_fix_scale_sd %>%
  select(year,species,tot.q10, tot.q25,tot.q50,tot.q75,tot.q90) %>%
  distinct() %>% left_join(RMU_species_number_ts, by = c("species" = "CommonName"))

# read in the IUCN pop for coverage sorting in the plotting
seal_IUCN=readxl::read_excel("data/seal_IUCN_gen_length.xlsx") %>%
  dplyr::select(Our_name, IUCN.pop,IUCN.pop.year,Trend,IUCN.cut,Climate_Zone) %>% mutate(IUCN.pop = as.numeric(IUCN.pop))

cov_species_tbl = index_df_tmp_n_species %>% 
  inner_join(seal_IUCN, by = c("species" = "Our_name", "year" = "IUCN.pop.year")) %>%
  mutate(coverage = tot.q50/IUCN.pop*100,
         coverage = pmin(coverage,100),
         coverage = ceiling(coverage)) %>% # x.tot.norm if no posterior pred
  mutate(coverage_cat = case_when(coverage <= 0.33 ~ "low",
                                  coverage > 0.33 & coverage <= 0.67 ~ "med",
                                  coverage >0.67 ~ "high")) %>%
  ungroup()%>%
  dplyr::select(species,coverage_cat,coverage,Trend,IUCN.cut,Climate_Zone)

index_df_tmp_n_species_cov = index_df_tmp_n_species %>% left_join(cov_species_tbl, by = "species") %>%
  mutate(Units_clean = case_when(Units_clean == "All Individuals" ~ "I",
                                 Units_clean == "Pups" ~ "P",
                                 Units_clean == "Adults" ~ "A",
                                 Units_clean == "Non-pups" ~ "NP"))

index_df_tmp_n_species_cov$species_label <- paste(index_df_tmp_n_species_cov$species,", ",
                                              "p = ",  index_df_tmp_n_species_cov$species_n,", ",
                                              " \n ",
                                              "c = ",  index_df_tmp_n_species_cov$coverage_cat,", ",
                                              "units = ", index_df_tmp_n_species_cov$Units_clean,
                                              sep="")

index_df_tmp_n_species_cov$species_label <- paste(index_df_tmp_n_species_cov$species,", ",
                                                  index_df_tmp_n_species_cov$Units_clean,
                                                  sep="")

best_data_means = best_pass2_model_fits_error_fix_scale_sd %>%
  mutate(data_norm = exp(data)) %>%
  group_by(ID) %>% 
  mutate(data_mean = mean(data_norm,na.rm=T)) %>%
  dplyr::select(species,year,ts,data_norm,data_mean) %>%
  drop_na()

index_means = index_df_tmp_n_species_cov %>% 
  group_by(species) %>%
  mutate(index_mean = mean(tot.q50,na.rm=T)) %>%  # x.tot.norm
  dplyr::select(species,index_mean) %>% distinct()

best_data_means_trans = best_data_means %>% left_join(index_means) %>%
  mutate(data_norm_trans = data_norm*(index_mean/data_mean))

index_df_tmp_n_species_wdata = best_data_means_trans %>% right_join(index_df_tmp_n_species_cov)

# ordering for ggplot
index_df_tmp_n_species_wdata$coverage_cat <- factor(index_df_tmp_n_species_wdata$coverage_cat, levels = c("high", "med","low"))
index_df_tmp_n_species_wdata_species_order = index_df_tmp_n_species_wdata %>% 
  # arrange(desc(species_n)) %>%
  arrange(desc(coverage))

index_df_tmp_n_species_wdata$Trend <- factor(index_df_tmp_n_species_wdata$Trend, levels = c("Increasing", "Decreasing","Stable or Unknown"))
index_df_tmp_n_species_wdata_species_order = index_df_tmp_n_species_wdata %>% 
  # arrange(desc(species_n)) %>%
  arrange(Trend,desc(coverage))

index_df_tmp_n_species_wdata_species_order = index_df_tmp_n_species_wdata %>% 
  # arrange(desc(species_n)) %>%
  arrange(species)

# you need to break the following 2 steps up so that the order can be right
index_df_tmp_n_species_wdata_species_order = index_df_tmp_n_species_wdata_species_order %>%
  mutate(species_label = factor(species_label, levels = unique(index_df_tmp_n_species_wdata_species_order$species_label)))

unique(index_df_tmp_n_species_wdata_species_order$species_label)

# plot it
y_axis_divide = 1000
data_in = index_df_tmp_n_species_wdata_species_order %>%
  group_by(species) %>%
  mutate(ymax_val = max(tot.q75/y_axis_divide)) 

see = data_in %>%
  group_by(species) %>%
  filter(year >= IUCN.cut) %>%
  select(species,year,tot.q50) %>%
  distinct() %>%
  group_by(species) %>%
  summarise(r = summary(lm(log(tot.q50)~year))$coefficients[2,1],
            pval = summary(lm(log(tot.q50)~year))$coefficients[2,4]) %>%
  mutate(trend_cat = case_when(r < -0.0075 & pval < 0.01 ~ "Decreasing",
                               r > 0.0075 & pval < 0.01 ~ "Increasing",
                               .default = "Stable or Unknown")) %>%
  right_join(data_in) %>%
  mutate(trend_cat = replace_na(trend_cat, "Stable or Unknown"))

p= see %>%
  ggplot() +

  # geom_ribbon(aes(x=year,ymin = pmax((x.tot.norm-year_species_sd)/y_axis_divide,0), ymax = (x.tot.norm+year_species_sd)/y_axis_divide), alpha=0.4) +
  # geom_line(aes(x=year,y=(x.tot.norm/y_axis_divide)),lwd=0.9) +
  
  # geom_ribbon(aes(x=year,ymin = tot.q01/y_axis_divide, ymax = tot.q90/y_axis_divide), alpha=0.2) +
  geom_vline(aes(xintercept=IUCN.cut),alpha=0.5,linetype="dotted") +
  geom_ribbon(aes(x=year,ymin = tot.q25/y_axis_divide, ymax = tot.q75/y_axis_divide,fill=trend_cat), alpha=0.4) +
  geom_line(aes(x=year,y=(tot.q50/y_axis_divide),color=trend_cat),lwd=0.9) +
  scale_color_manual(values = c( "#610B0B","#0E5917","#4D4F70"))+
  scale_fill_manual(values = c("#610B0B","#0E5917","#4D4F70"))+

  scale_x_continuous(breaks=scales::pretty_breaks(n=4))+
  scale_y_continuous(breaks=scales::pretty_breaks(n=3))+
  labs(x="Year",y=paste("Monitored Abundance (1000s units)",sep="")) +
  # geom_hline(yintercept=0, linetype = "dashed", alpha=0.3) +
  facet_wrap(~species_label, scales = "free", ncol = 4) +
  # expand_limits(y = c(0,data_in$ymax_val)) + 
  coord_cartesian(y=c(0,NA)) + 
  # xlim(c(1965,2022)) +
  # theme_bw() +
  theme(strip.background = element_rect(colour=NA, fill=NA),
        strip.text = ggfacet::element_textbox_highlight(
          size = 9, face = "plain",
          fill = NA, box.color = NA, color = "gray10",
          halign = .5, linetype = 1, r = unit(0, "pt"), width = unit(1, "npc"),
          padding = margin(2, 0, 1, 0), margin = margin(0, 1, 3, 1)
          ),
        panel.border = element_rect(colour = "gray10", fill=NA, size=0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color = "gray10"),
        axis.text.y = element_text(color = "gray10"),
        axis.title.x = element_text(color = "gray10"),
        axis.title.y = element_text(color = "gray10"),
        panel.background = element_rect(fill = "transparent",colour = NA),
        legend.position = c(0.75,0.05),legend.direction="horizontal",
        # plot.title = element_text(hjust = .01, vjust=-7),
        plot.margin = unit(c(0, 0.5, 0, 0.5), "cm")
        )
p
pdf("model/plots/seals_species_trends_marss_portrait.pdf", width = 8.5, height = 7)
print(p)
dev.off()

