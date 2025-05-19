# run after MARSS_model_estimation_sim_real
compare_model_fits = pred_data 

compare_model_fits_dist = compare_model_fits %>% 
  group_by(ts,model_type) %>%
  mutate(distance = NA)

ts_fits_distance = NULL
modeltype_vec_loop = unique(compare_model_fits_dist$model_type)
for (j in 1:length(modeltype_vec_loop)){
  modeltype_tmp = modeltype_vec_loop[j]
  compare_model_fits_dist_model_type = compare_model_fits_dist %>% filter(model_type == modeltype_tmp)
ts_vec_loop = unique(compare_model_fits_dist_model_type$ts)
for (i in 1:length(ts_vec_loop)){
  ts_tmp = ts_vec_loop[i]
  ts_fits_tmp = compare_model_fits_dist_model_type %>% filter(ts == ts_tmp)
  
  a <- which(!is.na(ts_fits_tmp$data)) 
  b <- which(is.na(ts_fits_tmp$data))
  
  for (i in 1:length(b)){
    dist_tmp = min(abs(a - b[i]),na.rm=T)
    ts_fits_tmp$distance[b[i]] = dist_tmp
  }
  
  ts_fits_distance = rbind(ts_fits_distance,ts_fits_tmp)
}
}

# fix the ballooning errors
compare_model_fits_error_fix = ts_fits_distance %>% ungroup() %>% # 
  mutate(sd_tmp = sd) %>%
  mutate(sd_tmp = case_when(distance > 5 ~ NA,
                            .default = sd_tmp)) %>%
  dplyr::group_by(ts) %>%
  fill(c(sd_tmp), .direction = "downup") #change to up only if necessary


compare_model_species_index = compare_model_fits_error_fix %>%
  mutate(ID = as.numeric(stringr::word(ts, 2, sep = "_"))) %>%
  left_join(scaling_tbl) %>%
  filter(Indiv.sum == 1) %>%
  mutate(norm_mean = exp((mean)), #mode now, make mean: +0.5*sd_tmp^2
         normal_sd = norm_mean*sqrt(exp(sd_tmp^2)-1)) %>% 
  mutate(state_fit_norm = norm_mean*Scaling,
         state_sd = normal_sd*Scaling) %>%
  group_by(year,species,model_type) %>%
  summarise(year=mean(year),
            x.tot.norm=sum(state_fit_norm,na.rm=T),
            # lo.tot.norm = sum(state_lo_norm,na.rm=T),
            # hi.tot.norm=sum(state_hi_norm,na.rm=T),
            year_species_sd=sqrt(sum(state_sd^2,na.rm=T))) %>%
  group_by(species,model_type) %>%
  mutate(x.tot.norm.mean = x.tot.norm/mean(x.tot.norm))

# SWITCH TO POSTERIOR PRED
compare_model_species_index = pred_data %>%
  select(model_type,year,species,tot.q50) %>%
  distinct() %>%
  group_by(species,model_type) %>%
  mutate(meanscale = tot.q50/mean(tot.q50)) 

lpi_tbl_all = readRDS("model/outputs/lpi_tbl_realdata.RDS")
lpi_to_join_with_MARSS = lpi_tbl_all %>%
  filter(LPI_final > 0) %>%
  filter(weights == 1) %>%
  group_by(species) %>%
  mutate(meanscale = LPI_final/mean(LPI_final)) %>%
  mutate(model_type = "LPI",year=as.numeric(year)) %>%
  select(species,year,meanscale,model_type)

MARSS_LPI_species_trend = rbind(compare_model_species_index,lpi_to_join_with_MARSS) %>%
  mutate(model_type = case_when(model_type == "GP" ~ "Synchrony Model",
                                   model_type == "Ueq" ~ "Independent Model",
                                   model_type == "LPI" ~ "LPI"))
  

cols <- c("Synchrony Model" = "#186689", "LPI" = "forestgreen", 
          "Independent Model" = "#e23b0d")


p = MARSS_LPI_species_trend %>%
  ggplot(aes(color = model_type, linetype = model_type)) +
  # scaled units
  # geom_point(aes(x=year,y=data_norm_trans),alpha=0.2)+
  # geom_line(aes(x=year,y=data_norm_trans,group=ID),alpha=0.1)+
  # geom_ribbon(aes(x=year,ymin = pmax((x.tot.norm-year_species_sd)/y_axis_divide,0), ymax = (x.tot.norm+year_species_sd)/y_axis_divide), alpha=0.4) +
  geom_line(aes(x=year,y=(meanscale)),lwd=0.9) +
  scale_x_continuous(breaks=scales::pretty_breaks())+
  scale_y_continuous(breaks=scales::pretty_breaks())+
  labs(x="Year",y=paste("Mean scaled abundance",sep="")) +
  scale_color_manual(values = cols) + 
  # geom_hline(yintercept=0, linetype = "dashed", alpha=0.3) +
  facet_wrap(~species, scales = "free", ncol = 4) +
  labs(color = "Model Type", linetype = "Model Type") +
  # expand_limits(y = c(0)) + 
  # xlim(c(1965,2022)) +
  theme_classic() +
  theme(legend.position = c(0.8, 0.05))
p

pdf("model/plots/seal_compare_model_species_trends.pdf", width = 8.5, height = 11)
print(p)
dev.off()
