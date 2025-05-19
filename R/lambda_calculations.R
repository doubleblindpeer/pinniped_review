#=== plot lambda trends across years
# join to generation length data, IUCN table is the middle step to join
seal_IUCN=readxl::read_excel("data/seal_IUCN_gen_length.xlsx")
seal_gen_length=readxl::read_excel("data/pinniped_generation_lengths.xlsx")

IUCN_genlength = seal_IUCN %>% 
  left_join(seal_gen_length, by = c("GL_name"="Scientific_name")) %>%
  mutate(GenerationLength_y = GenerationLength_d/365) %>%
  select(Our_name,GenerationLength_y,Realm,Taxa_Group, IUCN_Status,Climate_Zone,`Current Population Trend`)

lambda_year = index_df_tmp_n_species %>% #best_pass2_model_fits from model sel 2
  left_join(IUCN_genlength,by=c("species"="Our_name")) %>% 
  group_by(species,GenerationLength_y,Realm,Taxa_Group, IUCN_Status,Climate_Zone) %>% 
  arrange(year) %>%
  mutate(lambda = (tot.q50-lag(tot.q50))/lag(tot.q50))

pdf("model/plots/lambda_years.pdf", width = 5, height = 3)
lambda_year_boxplot_n = lambda_year %>% group_by(year) %>%
  dplyr::mutate(n_species = n()) %>% filter(lambda<0.3 & lambda > -0.3)

p=ggplot(lambda_year_boxplot_n, aes(y = lambda*100, x = year)) +
  coord_cartesian(ylim=c(-5,9),xlim=c(1950,2020))+ # using coord cartesian to avoid changing boxplots
  # xlim(c(1950,2024))+
  geom_hline(yintercept=0,color="grey50") + 
  geom_boxplot(aes(group=year,fill=n_species),outlier.shape=NA, coef = 0,linewidth=0.2,color="white") +
  viridis::scale_fill_viridis(option="mako",direction=-1) +
  stat_summary(aes(group=year),fun.y=mean, geom="point", shape=20, size=1.5, color="#E9724C") +
  geom_smooth(color="#C5283D",fill="#C5283D",se=T,level = 0.95,lwd=1.5,alpha=0.2,
              method = loess, method.args = list(family = "symmetric")) +
  theme_classic() + 
  labs(y="Annual growth rate (%)",x="Year",fill="Number of species") +
  # guides(guide_legend("Number of species"))+
  theme(
    legend.position = c(0,1), # top left position
    legend.justification = c(-.1, 1.1), # top left justification
    legend.key.size = unit(0.3, 'cm')
    # legend.box.margin = margin(5, l = 5, unit = "mm") # small margin
  )
print(p)
dev.off()
p

#==== plot species lambdas overall

# function for fitting r
years=why$year
model.trend.scales = why$tot.q50
fit_r = function(years,model.trend.scales){
  # get variables
  Nyear = last(years) - first(years) + 1
  N1_start = first(model.trend.scales)
  N1_start = mean(model.trend.scales[1:round(length(model.trend.scales)/5)])
  
  r_start = (log(last(model.trend.scales)) - log(first(model.trend.scales)))/Nyear
  
  output_year = first(years):last(years)
  
  # like function
  Likelihood = function(pars){ 
    N1 = pars[1]
    r = pars[2]
    sigma = pars[3]
    
    output <- rep(NA, Nyear)
    output[1] <- N1
    for (i in 1:(Nyear-1)) {
      n.t <- output[i]
      n.t.plus.1 <- n.t + n.t * r 
      n.t.plus.1 <- max(1e-5,n.t.plus.1)
      output[i+1] <- n.t.plus.1
    }
    
    nll  <- -sum(dnorm(model.trend.scales,output[which(output_year%in%years)],sd=sigma),na.rm=T,log=T)
    
    if (is.na(nll)){nll=1e7}
    return(nll)
  }
  start_pars=c(N1_start,r_start,sd(model.trend.scales)/4)
  # start_pars=c(N1_start*2,r_start)
  result = nlminb(start=start_pars,Likelihood,
                  control=list(eval.max=1e5, iter.max=1e5))
  if (result$convergence ==0){
  std <- sqrt(diag(solve(numDeriv::hessian(Likelihood, result$par))))}
  else{ std = 100}
  r_std_vec = c(result$par[2],std[2])
  
  return(r_std_vec)
  
}

temp = lambda_species %>% group_by(species) %>%
  summarise(ymax = max(year), ymin = min(year)) 

lambda_species = lambda_year %>%
  arrange(year) %>%
  group_by(species,IUCN_Status)%>%
  # filter(year > max(year)-3*GenerationLength_y)%>%
  # filter(year < max(year))%>%
  # mutate(n_year = n()) %>%
  # filter(n_year >9) %>%
  arrange(species,year) %>%
  summarise(
            r_fit = fit_r(year,tot.q50)[1],
            r_se = fit_r(year,tot.q50)[2],
            mean_lambda = mean(lambda,na.rm=T)*100,
            se_lambda = sd(lambda,na.rm=T)/sqrt(sum(!is.na(lambda)))*100,
            year_tot = max(year)-min(year)+1,
            year_min= min(year),
            year_max = max(year)) %>%
  mutate(trend_cat = case_when((mean_lambda + 1.28*se_lambda) > 0 & (mean_lambda - 1.28*se_lambda) > 0  ~ "Increasing",
         (mean_lambda + 1.28*se_lambda) < 0 & (mean_lambda - 1.28*se_lambda) < 0  ~ "Decreasing",
        .default = "Unknown"))

lambda_pops=best_pass2_model_fits_error_fix_scale_sd %>% #index_df_tmp
  left_join(IUCN_genlength,by=c("species"="Our_name")) %>% 
  group_by(ts,GenerationLength_y,Realm,Taxa_Group, IUCN_Status,Climate_Zone) %>% 
  arrange(year) %>%
  mutate(lambda = (exp(mean)-lag(exp(mean)))/lag(exp(mean))) %>%
  group_by(species,ts)%>%
  filter(data > -99) %>%
  # filter(year > max(year)-3*GenerationLength_y)%>%
  # mutate(n_year = n()) %>%
  # filter(n_year >9) %>%
  summarise(
            # r_fit = fit_r(year,tot.q50)*100,
            mean_lambda = median(lambda,na.rm=T)*100,
            se_lambda = sd(lambda,na.rm=T)/sqrt(sum(!is.na(lambda)))*100,
            year_tot = max(year)-min(year)+1,
            year_max = max(year)) %>%
  filter(abs(mean_lambda) < 15)

lambda_pops %>% filter(species == "Hawaiian monk seal")

lambda_species_plot = ggplot(data =lambda_species) + 
  geom_histogram(aes(x=mean_lambda)) + 
  geom_vline(xintercept = 0, color = "black",alpha=0.8, linetype = "dashed") +
  labs(x = "Total species abundance growth rate (%)",y = "Count") + theme_classic() +
  theme_classic()   
lambda_pops_plot = ggplot(data =lambda_pops) + 
  geom_vline(xintercept = 0, color = "black",alpha=0.8, linetype = "dashed") + 
  geom_boxplot(aes(x=mean_lambda,y=reorder(species, mean_lambda, median)),lwd=0.5) + 
  theme_classic() +  
  labs(x = "Population abundance growth rate (%)",y = "Species")

lambda_plots = cowplot::plot_grid(lambda_species_plot, lambda_pops_plot, ncol = 1, 
                                  rel_widths = c(1, 1),
                                  rel_heights = c(0.5, 1),
                                  labels=c("a","b"))

pdf("model/plots/lambda_species_population.pdf", width = 6, height = 6)
print(lambda_plots)
dev.off()


lambda_species$IUCN_Status <- factor(lambda_species$IUCN_Status, levels = c("Least Concern", "Near Threatened", "Vulnerable", "Endangered"))
lambda_species$`Current Population Trend` <- factor(lambda_species$`Current Population Trend`, 
                                                    levels = c( "Unknown","Increasing","Stable", "Decreasing"))


#== by IUCN Status
pdf("model/plots/lambda_IUCN_status.pdf", width = 8, height = 6)
p=ggplot(data=lambda_species,aes(x=mean_lambda,y=IUCN_Status, color=IUCN_Status)) + 
  geom_vline(xintercept = 0, color = "black",alpha=0.8, linetype = "dashed") + 
  geom_boxplot(alpha=1,outlier.shape = NA) + 
  # geom_jitter(aes(size=year_tot),alpha=0.8) + 
  geom_jitter(size=3,alpha=0.8) +
  theme_classic() +  
  scale_colour_manual(name = "IUCN Status",values = c("#9ACD9B","#016766","#CD9A02","#cd6630")) +
  # scale_radius(name = "Data Year Range",range=c(1,10)) +
  theme(
    legend.position = c(1,1), # top left position
    legend.justification = c(1.2, 1.1)) + 
  labs(x = "Total species abundance growth rate (%)",y = "IUCN Status")
print(p)
dev.off()

hist(lambda_species$mean_lambda)

pdf("model/plots/lambda_IUCN_trend.pdf", width = 8, height = 6)
p=ggplot(data=lambda_species,aes(x=mean_lambda,y=`Current Population Trend`, color=`Current Population Trend`)) + 
  geom_boxplot(alpha=1,outlier.shape = NA) + 
  geom_jitter(aes(size=year_tot),alpha=0.8) + theme_classic() +  
  scale_colour_manual(name = "IUCN Trend",values = c("#d6ccc2","#a7c957","#e9c46a","#bc4749")) +
  scale_radius(name = "Data Year Range",range=c(1,10)) +
  theme(
    legend.position = c(1,1), # top left position
    legend.justification = c(1.2, 1.1)) + 
  labs(x = "Species growth rate",y = "IUCN Current Trend")
print(p)
dev.off()