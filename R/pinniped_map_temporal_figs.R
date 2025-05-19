library(tidyverse)
library(ggplot2)
library(MARSS)
library(sf)
library(ggthemes) 
library(rnaturalearth) # for map
library(here)
library(ggspatial) #for scale bar and map arrow

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

counts_per_year = merged_long %>% drop_na(Count) %>%
  group_by(Year,Collected.by) %>%
  summarise(n_count = n())

#==== Temporal coverage plot
length(unique(merged_long$CommonName))/2

p = merged_long %>% drop_na(Count) %>% 
  ggplot(aes(CommonName, Year,group=ID,color=CommonName,fill=CommonName)) +
  geom_line(alpha=0.4,position=position_dodge(width = 1)) +
  geom_point(alpha=0.9,position=position_dodge(width = 1),pch=21, size=1.2) + 
  theme_bw() +
  coord_flip() +
  scale_x_discrete(limits=rev)+
  scale_color_manual(values=rep(c("#493843", "#61988E"),length(unique(merged_long$CommonName))/2))+
  scale_fill_manual(values=rep(c("#493843", "#61988E"),length(unique(merged_long$CommonName))/2))+
  ylim(c(1950,2020)) + 
  labs(x="Species")+
  theme(panel.spacing = unit(0, units = "cm"), # removes space between panels
       strip.placement = "outside", # moves the states down
       strip.background = element_rect(fill = NA),legend.position="none")# removes the background from the state names

pdf("model/plots/pinniped_temporal_coverage.pdf", width = 8.5, height = 8)
print(p)
dev.off()



#==== Map plot 
merged_long <- pinn_df %>% mutate(ID = row_number()) %>%
  pivot_longer(
    cols = starts_with("y"),
    names_to = "Year",
    names_prefix = "y",
    values_to = "Count",
    values_drop_na = TRUE) %>%
  mutate(Count=as.numeric(Count),Year=as.numeric(Year)) %>% drop_na(Count)

Combo_long = merged_long %>%  
  group_by(ID, Family,Collected.by,Use, Longitude,Latitude) %>%
  summarise(n=n()) %>%
  mutate(Lon = Longitude,
         Lat = Latitude)

Combo_notUsed = Combo_long %>%
  filter(Use == "no")

Combo_Used = Combo_long %>%
  na.omit() %>%
  filter(Use == "yes") 

# Make into shape files 

latlong_bind_sf_notUsed = Combo_notUsed %>% 
  st_as_sf(coords = c("Lon","Lat"))%>% 
  st_set_crs(4326) #sets the coordinates to WGS 84

latlong_bind_sf_Used = Combo_Used %>% 
  st_as_sf(coords = c("Lon","Lat"))%>% 
  st_set_crs(4326) #sets the coordinates to WGS 84

# get world map
world <- ne_countries(scale = "medium", returnclass = "sf")

# Kovacs shapefile
kovacs_sf <- st_read(here("data", "Shapefile", "Kovacs.shp")) # Replace with your actual path

# ensure that Kovacs is in the same CRS as world
kovacs_sf <- st_transform(kovacs_sf, st_crs(world))

#Hacky way to get the legend to work
kovacs_sf$fillcol <- factor("Kovacs 2011 Global Pinniped Extent", 
                         levels = "Kovacs 2011 Global Pinniped Extent")

p <- ggplot() +
  geom_sf(data = kovacs_sf, aes(fill = "#dcecf7"), color = NA, alpha = 1) + # Kovacs sf with legend MUST COME FIRST
  geom_sf(data = world, fill = "grey70", color = NA) + # World map
  geom_sf(data = latlong_bind_sf_Used, aes(size = n, color = Family), alpha = 0.7) + 
  geom_sf(data = latlong_bind_sf_notUsed, aes(shape = "Data not used"), color = "black", alpha = 0.5) + 
  theme_bw() +
  theme(
    # plot.margin=unit(c(0,0,0,0), "pt"),
    legend.position = "bottom",
    legend.text = element_text(size = 8),
    # legend.box = "vertical", # Set legend box arrangement
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
        )+
          #legend.box.background = element_rect(color = "grey", size = 0.5)) +
  coord_sf(crs = st_crs(4326)) +
  #set map limits
  coord_sf(xlim = c(-165, 170),ylim = c(-82, 80))+
  # annotation_scale(location = "bl", width_hint = 0.3) + #I think this looked weird.
  # annotation_north_arrow(location = "bl", which_north = "true", 
  #                        pad_x = unit(2, "cm"), pad_y = unit(2, "cm"),
  #                        height = unit(1, "cm"), width = unit(1, "cm"),
  #                        style = north_arrow_fancy_orienteering) +
  scale_size_continuous(name = "Number of observations") +
  scale_fill_manual(values = "#c7e3f0", name = "", labels = "Global pinniped extent") +
  labs(color = "", size = "Count of Seals") +
  scale_shape_manual(values = c(4), labels = "Unused data",name="") +
  scale_color_manual(values = c("#4F359B","#054A29", "#BE5A38")) + # Replace with your actual family names and desired colors
  guides(colour = guide_legend(override.aes = list(alpha = 1,size=5))) +
  theme(panel.grid.major = element_blank())
p

pdf("model/plots/pinniped_world_map_2024_1014.pdf", width = 11, height = 8.5)
print(p)
dev.off()
