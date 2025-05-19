library(sf)
library(gdistance)
library(tidyverse)
library(stars)
library(terra)
library(raster)

# read in data
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


# set up the distance shapefiles & rasters
sf_use_s2(T)
ocean = st_read("data/Shapefile/ne_10m_ocean/ne_10m_ocean.shp") %>%
  st_set_crs(4326) %>%
  mutate(value = 1)

r <- rast(ext(ocean), nrow = 1000, ncol = 1000)
ocean_conductive_rast <- raster::rasterize(vect(ocean), r, field=1,background = 0) 
ocean_rastlayer = raster(ocean_conductive_rast)
plot(ocean_rastlayer)
ocean_transition <- gdistance::transition(ocean_rastlayer, transitionFunction = mean, directions = 8,
                                          symm= T)
ocean_transition <- geoCorrection(ocean_transition, type="c")

#==== Calculate shortest ocean path length
distance_df = NULL
length_shapefiles = NULL
Name_RMU = turt_IDs %>% dplyr::select(CommonName,RMU) %>%
  distinct()
i=20
for (i in 1:nrow(Name_RMU)){ 
  CommonName_tmp = Name_RMU$CommonName[i]
  RMU_tmp = Name_RMU$RMU[i]
  cat(CommonName_tmp,"\n",RMU_tmp,"\n")
  species_i_coords = turt_IDs %>% 
    filter(CommonName == CommonName_tmp,
           RMU == RMU_tmp)
  
  nsites = nrow(species_i_coords)
  
  turt_df_species_IDs_tmp = species_i_coords %>% pull(ID) %>% 
    # sort for easier matching in for loop later
    sort() %>% unique()
  nsites = length(turt_df_species_IDs_tmp)
  
  # distance_matrix = matrix(data = NA, nrow = nsites, ncol = nsites)
  
  for(site_i in 1:nsites){
    # distance_matrix[site_i,site_i] = 0
    if(site_i == nsites){next}
    for (site_j in (site_i + 1):nsites) {
    
      site_i_ID = turt_df_species_IDs_tmp[site_i]
      site_j_ID = turt_df_species_IDs_tmp[site_j]
      
      turt_df_long_i = species_i_coords %>% filter(ID == site_i_ID) 
      turt_df_long_j = species_i_coords %>% filter(ID == site_j_ID) 
    
      sf_use_s2(T)
      site_i_coords = turt_df_long_i %>% dplyr::select(Longitude,Latitude) %>% as.matrix()
      site_j_coords = turt_df_long_j %>% dplyr::select(Longitude,Latitude) %>% as.matrix()
      
      distance <- shortestPath(ocean_transition, 
                               site_i_coords,site_j_coords, 
                               output = "SpatialLines")
      
      distance_sf = st_as_sf(distance) %>% st_set_crs(4326) %>%
        mutate(row_id = site_i_ID,col=site_j_ID)
      distance_sf = st_wrap_dateline(distance_sf) 
      
      distance_tmp = st_length(distance_sf)
      
      distance_df_tmp = data.frame(row_id = site_i_ID, col_id=site_j_ID, 
                              CommonName =  CommonName_tmp, RMU = RMU_tmp, distance_m = distance_tmp)
      
      length_shapefiles = rbind(length_shapefiles,distance_sf)
      
      distance_df=rbind(distance_df,distance_df_tmp)
      
    }
  }
}

# saveRDS(distance_df, "data/site_distance/seal_distance_long_subsp_2024_1218.RDS")
# st_write(length_shapefiles, "data/site_distance/seal_length_shapefiles_subsp_2024_1218.shp", append=F)

#=== Look at empirical evidence for spatial covariance
turt_df_IDs_cc_rmu_tmp =  turt_df_long %>% pull(ID) %>% unique() %>% sort()
nsites = length(turt_df_IDs_cc_rmu_tmp)
cor_df = NULL

for(site_i in 1:nsites){
  cat(site_i,"\n")
  # distance_matrix[site_i,site_i] = 0
  # if(site_i == nsites){next}
  for (site_j in (site_i + 1):nsites) {
    
    site_i_ID = turt_df_IDs_cc_rmu_tmp[site_i]
    site_j_ID = turt_df_IDs_cc_rmu_tmp[site_j]
    
    turt_df_long_i = turt_df_long %>% filter(ID == site_i_ID) %>%
      dplyr::select(ID, CommonName, RMU, Year, Count)
    
    turt_df_long_j = turt_df_long %>% filter(ID == site_j_ID) %>%
      dplyr::select(ID, CommonName, RMU,Year, Count)
    
    turt_df_long_i_j = turt_df_long_i %>% left_join(turt_df_long_j, by = "Year") %>%
      drop_na()
    
    ncount_tmp = nrow(turt_df_long_i_j)
    if(ncount_tmp<2){next}
    
    lm_tmp = lm(data = turt_df_long_i_j, log(Count.x+0.01)~log(Count.y+0.01))
    
    corcoeff_tmp = cor(log(turt_df_long_i_j$Count.x+0.01),log(turt_df_long_i_j$Count.y+0.01))
    
    cov_tmp = cov(log(turt_df_long_i_j$Count.x+0.01),log(turt_df_long_i_j$Count.y+0.01))
    
    summ_lm_tmp = summary(lm_tmp)
    
    if(nrow(summ_lm_tmp$coefficients)<2){next}
    
    cor_tmp = summ_lm_tmp$coeff[2,1]
    pval_tmp = summ_lm_tmp$coeff[2,4]
    
    
    cor_df_tmp = data.frame(row_id = site_i_ID, col_id=site_j_ID, 
                            cor = cor_tmp, pval = pval_tmp, ncount = ncount_tmp,
                            corcoeff = corcoeff_tmp, cov= cov_tmp,
                            CommonName_row = unique(turt_df_long_i$CommonName),
                            RMU_row = unique(turt_df_long_i$RMU), 
                            CommonName_col = unique(turt_df_long_j$CommonName),
                            RMU_col = unique(turt_df_long_j$RMU))
    cor_df = rbind(cor_df,cor_df_tmp)
  }
}


cor_df = readRDS("data/site_distance/seal_cor_df_2024_1218.RDS")
cor_df_filter = cor_df %>% drop_na(pval) %>% filter(ncount > 0)

distance_df = readRDS("data/site_distance/seal_distance_long_subsp_2024_1218.RDS")

distance_cor_long = distance_df %>% left_join(cor_df_filter, by = c("row_id","col_id"))

distance_cor_long_cleans = distance_cor_long %>% drop_na(pval) %>%
  mutate(distance_km = as.numeric(distance_m)/1000)

# like function
Likelihood = function(pars){ 
  log_alpha = pars[1]
  log_rho = pars[2]
  log_sigma = pars[3]
  
  alpha = exp(log_alpha)
  rho = exp(log_rho)
  sigma = exp(log_sigma)
  
  pred_covar = (alpha)*exp(-1/2*(1/rho)*distance_cor_long_cleans$distance_km^2)
  
  nll  <- -sum(dnorm(distance_cor_long_cleans$cov,pred_covar,sigma,log=T),na.rm=T)
  
  if (is.na(nll)){nll=1e7}
  return(nll)
}

Likelihood(c(log(0.5),log(4000),log(0.2)))
start_pars = c(log(0.5),log(4000),log(0.2))
result = nlminb(start=start_pars,Likelihood,
                control=list(eval.max=1e5, iter.max=1e5))

alpha_hat = exp(result$par[1])
rho_hat = exp(result$par[2])
sigma_hat = exp(result$par[3])

# inputs for the simulation
alpha_hat = 0.9
rho_hat = 200000
sd_hat = 0.11

variance_ID = turt_df_long %>% 
  group_by(ID) %>%
  summarise(variance_ID = var(log(Count+0.01)),
            n=n()) %>%
filter(n>5)
  
#=== Set up Distance Matrix format for MARSS AND set up simulation covariance values
distance_df = readRDS("data/site_distance/seal_distance_long_subsp_2024_1218.RDS")
CommonName_vec = sort(unique(distance_df$CommonName))
distance_matrix_list = NULL
covariance_matrix_list = NULL
CommonName_i=4
for (CommonName_i in 1:length(CommonName_vec)){
  CommonName_tmp = CommonName_vec[CommonName_i]
  cat(CommonName_tmp,"\n")
  
  species_i_coords = turt_IDs %>% 
    filter(CommonName == CommonName_tmp)
  nsites = nrow(species_i_coords)
  
  species_i_IDs = species_i_coords %>% pull(ID) %>% sort()
  covariance_matrix = matrix(data = 0, nrow = nsites, ncol = nsites)
  rownames(covariance_matrix) = species_i_IDs
  colnames(covariance_matrix) = species_i_IDs
  
  distance_matrix = matrix(data = -99, nrow = nsites, ncol = nsites)
  rownames(distance_matrix) = species_i_IDs
  colnames(distance_matrix) = species_i_IDs
  
  for(site_i in 1:nsites){
    covariance_matrix[site_i,site_i] = 1
    distance_matrix[site_i,site_i] = 0
    if(site_i == nsites){next}
    for (site_j in (site_i + 1):nsites) {
  
      site_i_ID = species_i_IDs[site_i]
      site_j_ID = species_i_IDs[site_j]
      
      distance_m_tmp = distance_df %>% filter(row_id == min(site_i_ID,site_j_ID),
                             col_id == max(site_i_ID,site_j_ID)) %>%
        pull(distance_m)
      
      distance_km_tmp = as.numeric(distance_m_tmp/1000)
      if(length(distance_m_tmp) == 0){next} #if theres no distance data skip
      covariance_matrix[site_i,site_j] = alpha_hat*exp((-1/2)*(1/rho_hat)*distance_km_tmp^2)
      covariance_matrix[site_j,site_i] = covariance_matrix[site_i,site_j]
      
      distance_matrix[site_i,site_j] = distance_km_tmp
      distance_matrix[site_j,site_i] = distance_matrix[site_i,site_j]
  
    }
  }
  distance_matrix[site_i,site_i] = 0
  # turn the correlation matrix into a covariance matrix
  covariance_matrix_list[[CommonName_i]] = covariance_matrix*sd_hat*sd_hat
  names(covariance_matrix_list)[CommonName_i] = CommonName_tmp
  
  distance_matrix_list[[CommonName_i]] = distance_matrix
  names(distance_matrix_list)[CommonName_i] = CommonName_tmp
}

# saveRDS(covariance_matrix_list, "data/site_distance/seal_covariance_matrix_subsp_2025_0131.RDS")
covariance_matrix_list = readRDS("data/site_distance/seal_covariance_matrix_subsp_2025_0131.RDS")


# start the simulation
start_pops_year1 = turt_df_long %>% 
  group_by(ID) %>%
  summarise(log_Countmean_y1 = log(mean(Count)))

N_sims = 10
sim_i = 1
phi = 0 #0.7
GP= T
true_state_name = "GP"
for (sim_i in 1:N_sims){
  set.seed(123+sim_i)
Year_seq = min(turt_df_long$Year):max(turt_df_long$Year)
# Year_seq = 1:2000
sim_pop_matrix = NULL
for ( i in 1:length(covariance_matrix_list)){
  print(i)
  covariance_matrix = covariance_matrix_list[[i]]
  covariance_matrix_r = covariance_matrix

  # TOGGLE ON FOLLOWING TO TURN OFF COVARIANCE
  if (GP==F){
  covariance_matrix_r = diag(diag(covariance_matrix))
  }
  
  # year 1
  row_IDs = as.numeric(rownames(covariance_matrix)) %>% sort()
  sim_pops_year_tmp = start_pops_year1 %>% filter(ID %in% row_IDs) %>% pull(log_Countmean_y1) %>%
    data.frame()
  
  growth_rates_tmp_tminus1 = as.vector(mvtnorm::rmvnorm(n=1,sigma = covariance_matrix_r))  
  # year 2 and on
  for (y_tmp in 1:(length(Year_seq)-1)){
    growth_rates_tmp = as.vector(mvtnorm::rmvnorm(n=1,
                                                  mean = rep(0,nrow(covariance_matrix_r)),
                                                  sigma = covariance_matrix_r))  
  
  sim_pops_year_tmp[,y_tmp+1] = sim_pops_year_tmp[,y_tmp]+growth_rates_tmp
  growth_rates_tmp_tminus1 = growth_rates_tmp 
  }
  
  sim_pop_matrix_tmp = as.matrix(sim_pops_year_tmp)
  rownames(sim_pop_matrix_tmp) = row_IDs
  colnames(sim_pop_matrix_tmp) = Year_seq
  
  sim_pop_matrix = rbind(sim_pop_matrix, sim_pop_matrix_tmp)
}

sim_pop_matrix_long = sim_pop_matrix %>% as.data.frame() %>%
  rownames_to_column("ID") %>%
  pivot_longer(!ID, names_to = "Year", values_to = "Count_log") %>%
  mutate(ID = as.numeric(ID),
         Year = as.numeric(Year),
         Count_log = as.numeric(Count_log))

turt_df_long_attribute = turt_df_long %>% 
  dplyr::select(ID, Common.Name, Subspecies,ts_id,Site,Indiv.pop.map,Indiv.sum) %>% distinct()
sim_pop_matrix_long_true = sim_pop_matrix_long %>% left_join(turt_df_long_attribute)

# plot it to see all the RW MA1
sim_pop_matrix_long_true %>% 
  ggplot(aes(x=Year, y=Count_log,color=ID,group=ID)) + 
  geom_line(alpha=0.8) +
  facet_wrap(~Common.Name,scales="free")

# join real observation years to simulation data
Obs_years_df = turt_df_long %>% dplyr::select(ID, Year) %>%
  mutate(Obs_year = T) 

sim_pop_matrix_long_true_obs = sim_pop_matrix_long_true %>% 
  left_join(Obs_years_df, by = join_by(ID, Year)) %>%
  ungroup()%>%
  rowwise()%>%
  mutate(Count_log_obs = rnorm(n=1,mean=Count_log,sd=0.3))

# plot it to see all the observation error and RW MA1
p = sim_pop_matrix_long_true_obs %>% 
  ggplot(aes(x=Year, y=Count_log_obs,color=ID,group=ID)) + 
  geom_line(alpha=0.8) +
  facet_wrap(~Common.Name,scales="free")

print(p)

save_file_name = paste("data/simulation_GP_MA1/seal_sim_true_obsV5_",true_state_name,"_sim",sim_i,".RDS",sep="")
saveRDS(sim_pop_matrix_long_true_obs,save_file_name)
}




## PLOT EXAMPLE OF SEAL OCEAN DISTANCE COVARIANCE
# set up your example for harbor seals
ID_POI = 70
Subspecies_POI = turt_IDs %>% filter(ID == ID_POI) %>% pull(Subspecies)
vitulina_data = turt_IDs %>% filter(Subspecies == Subspecies_POI)

# get all points and POI
vitulina_coords_sf = vitulina_data %>%
  drop_na() %>%
  st_as_sf(coords = c("Longitude", "Latitude")) %>%
  st_set_crs(4326)

POI_dutch_wadden = vitulina_coords_sf %>% filter(ID == ID_POI)

# get shapefiles
seal_length_shapefiles = st_read("data/site_distance/seal_length_shapefiles_subsp_2024_0928.shp")
POI_lengths_dutch_wadden = seal_length_shapefiles %>% filter(row_id == ID_POI | col == ID_POI)

# get distances to shade color and plot 
distance_POI = distance_df %>% filter(row_id == ID_POI | col_id == ID_POI) %>%
  mutate(distance_km = as.numeric(distance_m/1000))
POI_lengths_dutch_wadden_distances = POI_lengths_dutch_wadden %>% left_join(distance_POI, by =c("col"= "col_id","row_id"))

other_distances_by_ID = distance_POI %>% mutate(other_id = case_when(row_id == ID_POI ~ col_id,
                                             col_id == ID_POI ~ row_id)) %>%
  dplyr::select(distance_km,other_id)

vitulina_coords_sf_distances = vitulina_coords_sf %>% left_join(other_distances_by_ID,by=c("ID" = "other_id"))

library(rnaturalearth)
world <- ne_countries(scale = "medium", returnclass = "sf") %>% 
  filter(admin != "Antarctica")

bbox = vitulina_coords_sf %>% st_transform(3832) %>% st_bbox()

POI_world_p = ggplot() +
  # geom_sf(data = ocean50, fill = "#114B5F") +
  geom_sf(data = world, fill = "grey60",color=NA) +
  geom_sf(data = st_jitter(POI_lengths_dutch_wadden_distances,amount=0.05), 
          aes(color = distance_km), alpha = 1,lwd=2.5,linetype = "2212") +
  geom_sf(data = st_jitter(vitulina_coords_sf_distances,amount=0), 
          aes(fill = distance_km), alpha = 1,size=5,shape=21,color="black")+
  geom_sf(data = POI_dutch_wadden, fill = "#A21111",color="black", alpha = 1,size=6,shape=21)+
  # geom_sf(data = species_i_coords_shp, color = "#F45B69", alpha = 1.5) +
  # coord_sf() +
  coord_sf(crs = st_crs(3832),xlim = c(bbox[1], bbox[3]),ylim = c(bbox[2], bbox[4]))+
  scale_colour_gradient2(low = "#3F5A3D", mid = "#6A994E",high = "#DACC3E",midpoint = 2.5,trans = scales::log_trans(base = 10))+
  scale_fill_gradient2(low = "#3F5A3D", mid = "#6A994E",high = "#DACC3E",midpoint = 2.5,trans = scales::log_trans(base = 10))+
  theme_bw()+
  theme(legend.position = c(0.2, 0.25),legend.key.size = unit(0.4, "cm"),legend.box.margin = margin(0, 0, 0, 0),
        panel.grid.major = element_blank(),
        axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()
        )
  
print(POI_world_p)

# compare the synchrony in mean trends

POI_species_i_pre_df <- turt_df_long %>% 
  filter(Subspecies == Subspecies_POI) %>%
  # filter(RMU == rmu_i) %>% # get only this ESU
  mutate(log.spawner = log(Count+1)) %>% # create a column called log.spawner
  group_by(ID) %>% 
  mutate(n=sum(!is.na(log.spawner))) %>%
  filter(Indiv.sum == 1) %>%
  dplyr::select(ID,Year,log.spawner) %>%
  mutate(log.spawner_meanscale = log.spawner/mean(log.spawner)) %>%
  filter(!ID %in% c(75,76,77,78,73,74,72))
 

other_species_i_pre_df_distances = POI_species_i_pre_df %>% left_join(other_distances_by_ID,by=c("ID" = "other_id"))%>%
  mutate(distance_km = distance_km+1) %>% drop_na()
POI_species_i_df= POI_species_i_pre_df %>% filter(ID == ID_POI)

# make nice colors
library(wesanderson)
library(ochRe)

pal = wes_palettes[["Zissou1Continuous"]]
pal = ochre_palettes[["williams_pilbara"]]
pal = ochre_palettes[["nolan_ned"]]
pal_in = colorRampPalette(pal)(4)

# here is the change:
# other_species_i_pre_df_distances$distance_km <- factor(
#   other_species_i_pre_df_distances$distance_km,
#   levels = sort(unique(other_species_i_pre_df_distances$distance_km),decreasing=T)) 

#create my color ramp function
YlGnBu <- colorRampPalette(RColorBrewer::brewer.pal(9, 'YlGnBu'))
scale_color_manual(values = YlGnBu(n_distinct(other_species_i_pre_df_distances$distance_km)))


POI_trends  = ggplot(data = other_species_i_pre_df_distances) + 
  geom_line(aes(x=Year,y=log.spawner_meanscale,group=ID,color=distance_km),alpha=0.9,size = 1,linetype="91") +
  geom_point(aes(x=Year,y=log.spawner_meanscale,group=ID,color=distance_km),alpha=0.9,size = 2) +
  geom_line(data = POI_species_i_df,
            aes(x=Year,y=log.spawner_meanscale),color="#A21111",lwd=2) +
  geom_point(data = POI_species_i_df,
            aes(x=Year,y=log.spawner_meanscale),color="#A21111",size=3) +
  labs(y="Log Count Mean-scaled")+
  # scale_color_gradientn(colours = rev(pal_in),trans = "log") +
  # scale_color_manual(values = YlGnBu(n_distinct(other_species_i_pre_df_distances$distance_km)))+
  # scale_alpha_continuous(trans = "log")+
  # scale_color_ochre(palette="williams_pilbara", discrete=FALSE) +
  # scale_color_gradient(trans = scales::log_trans(base = 3))
  # scale_color_distiller(trans = scales::log_trans(base = 3),palette = "Spectral")+
  scale_colour_gradient2(low = "#3F5A3D", mid = "#6A994E",high = "#DACC3E",midpoint = 2.5,trans = scales::log_trans(base = 10))+
  # scale_color_continuous(trans = scales::log_trans(base = 10),type = "viridis") + 
  theme_bw() +   theme(legend.position = "none") 
POI_trends
## CREATE COVARIANCE FUNCTION EXAMPLE
cov_plot3 = distance_cor_long_cleans_clip %>% filter(CommonName == "California sea lion") %>%
  filter(cor >= -0.5 & cor < 1) %>%
  ggplot() + 
  geom_point(aes(x=distance_km, y=cor),alpha=1) + 
  geom_function(data = data.frame(distance_km = 0, cov = 0, CommonName = "Australian sea lion"), lwd=2,
                fun = GP_fn,args = list(alpha = 0.333, rho = 81995, sigma = (0+1)^2), color = "#C99CD0") + 
  labs(y="Spatial covariance",x="Distance (km)")+
  theme(legend.position = "none") +
  theme_bw()


POI_compposite = cowplot::plot_grid(POI_world_p, POI_trends,cov_plot3, ncol = 1, rel_heights = c(1, 1,1),labels=c("a","b","c"))
POI_compposite

cowplot::plot_grid(POI_trends,cov_plot3, ncol = 1, rel_heights = c(1,1))


POI_compposite_horizontal = cowplot::plot_grid(POI_world_p, POI_trends,cov_plot3, ncol = 3, rel_widths = c(1, 1,1),labels=c("a","b","c"))
POI_compposite_horizontal
