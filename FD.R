library(ade4)
library(RVAideMemoire)
library(cluster)
library(clue)
library(sf)
library(ggplot)
library(plyr)
library(dplyr)

# load bird data

bird_data <- readRDS("Documents/bioagora/output/bird_data_clean.rds")

### SFI

# load trait-species matrix
trait <- read.csv("Téléchargements/life_history_bird2.csv",header = TRUE)

# species name consistency
trait$Species <- gsub("_"," ",trait$Species)
trait$Species[which(trait$Species == "Picus viridis")] <- "Picus viridis / Picus sharpei"
trait$Species[which(trait$Species == "Sylvia cantillans")] <- "Sylvia cantillans / Sylvia subalpina"
trait$Species[which(trait$Species == "Dendrocoptes medius")] <- "Leiopicus medius"
trait$Species[which(trait$Species == "Luscinia svecica")] <- "Cyanecula svecica"
trait$Species[which(trait$Species == "Cyanopica cyanus")] <- "Cyanopica cooki"

# select the same species as for the other indices
trait <- trait[which(trait$Species %in% unique(bird_data$sci_name_out)),]

# remove taxonomic info
trait <- trait[,-c(1:3)]

# remove traits not well covered
table.na <- apply(trait,2,function(x){sum(is.na(x))}) 
trait2 <- trait[,-which(table.na>12)]

# select traits
# all traits that could be relevant
trait_selected <- c("Species","LengthU_MEAN","WingU_MEAN","TailU_MEAN","BillU_MEAN","TarsusU_MEAN","WeightU_MEAN","Sexual.dimorphism",
                    "Clutch_MEAN","Broods.per.year","EggL_MEAN","EggW_MEAN","Egg_MASS","Young","Association.during.nesting","Nest.type",
                    "Nest.building","Mating.system","Incubation.period","Incubation.sex","Territoriality","Sedentary",
                    "Short.distance.migrant","Long.distance.migrant","Deciduous.forest","Coniferous.forest","Woodland","Shrub","Savanna",
                    "Tundra","Grassland","Mountain.meadows","Reed","Swamps","Desert","Freshwater","Marine","Rocks","Human.settlements","Folivore_B",
                    "Frugivore_B","Granivore_B","Arthropods_B","Other.invertebrates_B","Fish_B","Other.vertebrates_B","Carrion_B","Omnivore_B")
# the one used in Carmona et al. 2021 https://doi.org/10.1126/sciadv.abf2675
trait_selected <- c("Species","Clutch_MEAN","Broods.per.year","WeightU_MEAN","Incubation.period","Age.of.first.breeding","Fledging.period",
                    "Egg_MASS","LengthU_MEAN")
# the same but with diet as in Godet et al.
trait_selected <- c("Species","Clutch_MEAN","Broods.per.year","WeightU_MEAN","Incubation.period","Age.of.first.breeding","Fledging.period",
                    "Egg_MASS","LengthU_MEAN","Frugivore_B","Granivore_B","Arthropods_B","Other.invertebrates_B","Fish_B","Other.vertebrates_B","Carrion_B","Omnivore_B")
                    
trait3 <- trait2[,trait_selected]

# specify non noncontinuous traits
nc_traits <- c("Sexual.dimorphism","Young","Association.during.nesting","Nest.type",
               "Nest.building","Mating.system","Incubation.sex","Territoriality","Sedentary",
               "Short.distance.migrant","Long.distance.migrant","Deciduous.forest","Coniferous.forest","Woodland",
               "Shrub","Savanna","Tundra","Grassland","Mountain.meadows","Reed","Swamps","Desert","Freshwater","Marine","Rocks","Human.settlements",
               "Folivore_B","Frugivore_B","Granivore_B","Arthropods_B","Other.invertebrates_B","Fish_B","Other.vertebrates_B","Carrion_B","Omnivore_B")
nc_trait_selected <- trait_selected[which(trait_selected %in% nc_traits)]

trait3[,nc_trait_selected] <- 
  lapply(trait3[,nc_trait_selected], factor)

# transformation in distance matrix (gower)
asymm_trait <- c("Territoriality","Sedentary","Short.distance.migrant","Long.distance.migrant","Deciduous.forest","Coniferous.forest","Woodland","Shrub",
                 "Savanna","Tundra","Grassland","Mountain.meadows","Reed","Swamps","Desert","Freshwater","Marine","Rocks","Human.settlements",
                 "Folivore_B","Frugivore_B","Granivore_B","Arthropods_B","Other.invertebrates_B","Fish_B","Other.vertebrates_B","Carrion_B","Omnivore_B")
asymm_trait_selected <- which(trait3 %in% asymm_trait)-1

if(length(asymm_trait_selected) == 0){
  mat.dist <- daisy(trait3[,2:ncol(trait3)], metric="gower")
}else{
  mat.dist <- daisy(trait3[,2:ncol(trait3)], metric="gower", type=list(asymm=asymm_trait_selected))
}
attr(mat.dist, "Labels") <- trait3$Species

data_sfi <- data.frame(species = trait3$Species, sfi = apply(mat.dist,1,mean))


### SAI => should be rename SOI species occurrence index?

data_sai <- ddply(bird_data, .(sci_name_out), .fun = 
                    function(x){
                      # sites where the species occurs
                      site_occ <- unique(x$siteID[which(x$count>0)])
                      # focus on occurrence frequency within the sites where the species occurs
                      x_occ <- x[which(x$siteID %in% site_occ),]
                      # get occurrence instead of abundance
                      x_occ$occ <- ifelse(x_occ$count > 0, 1, 0)
                      # compute the frequency of occurence per site
                      x_sai <- data.frame(x_occ |> dplyr::group_by(siteID) |> dplyr::summarise(mean_oc = mean(occ)))
                      # return the mean of occurence frequency
                      return(data.frame(sai = mean(x_sai$mean_oc)))
                    }
                  )

### SRI, calcuated using atlas in Godet et al. 2015, here using convex hull (minimum polygon encompassing all sites of occurence)

# site location
site_mainland_sf_reproj <- readRDS("Documents/bioagora/output/site_mainland_sf_reproj.rds")

# map of European countries involved in bird monitoring
grid_eu_mainland_outline <- st_read("Documents/bioagora/output/grid_eu_mainland_outline.gpkg")

# get the max convex hull (i.e. mimic the case when a species was present in all site, to serve as a reference)

bird_convh <- st_convex_hull(st_union(site_mainland_sf_reproj))
bird_convh_inter <- st_intersection(bird_convh,grid_eu_mainland_outline)

ggplot(bird_convh_inter) + geom_sf()

# a few example of species with different distribution, to check the approach
TURMER_convh <- st_convex_hull(st_union(site_mainland_sf_reproj[which(site_mainland_sf_reproj$siteID %in% bird_data$siteID[which(bird_data$sci_name_out=="Turdus merula" & bird_data$count > 0)]),]))
TURMER_convh_inter <- st_intersection(TURMER_convh,grid_eu_mainland_outline)
as.numeric(st_area(TURMER_convh_inter))/as.numeric(st_area(bird_convh_inter))
ggplot(TURMER_convh_inter) + geom_sf()

SYLMEL_convh <- st_convex_hull(st_union(site_mainland_sf_reproj[which(site_mainland_sf_reproj$siteID %in% bird_data$siteID[which(bird_data$sci_name_out=="Sylvia melanocephala" & bird_data$count > 0)]),]))
SYLMEL_convh_inter <- st_intersection(SYLMEL_convh,grid_eu_mainland_outline)
as.numeric(st_area(SYLMEL_convh_inter))/as.numeric(st_area(bird_convh_inter))
ggplot(SYLMEL_convh_inter) + geom_sf()

CYACOO_convh <- st_convex_hull(st_union(site_mainland_sf_reproj[which(site_mainland_sf_reproj$siteID %in% bird_data$siteID[which(bird_data$sci_name_out=="Cyanopica cooki" & bird_data$count > 0)]),]))
CYACOO_convh_inter <- st_intersection(CYACOO_convh,grid_eu_mainland_outline)
as.numeric(st_area(CYACOO_convh_inter))/as.numeric(st_area(bird_convh_inter))
ggplot(CYACOO_convh_inter) + geom_sf()

# calculate for all species
data_sri <- ddply(bird_data, .(sci_name_out), .fun = 
                    function(x){
                      sp_convh <- st_convex_hull(st_union(site_mainland_sf_reproj[which(site_mainland_sf_reproj$siteID %in% x$siteID[which(x$count > 0)]),]))
                      sp_convh_inter <- st_intersection(sp_convh,grid_eu_mainland_outline)
                      sp_sri <- as.numeric(st_area(sp_convh_inter))/as.numeric(st_area(bird_convh_inter))
                      return(data.frame(sri = sp_sri))
                    },
                  .progress = "text"
                  ) 


### SGIc, not computed as it is based on species habitat which is not relevant for other taxa

### merge the three indices
sxi <- merge(data_sri,data_sai[,c("sci_name_out","sai")], by="sci_name_out")
sxi <- merge(sxi,data_sfi, by.x="sci_name_out", by.y="species",all.x=TRUE)

### check correlation
test_multicor <- round(cor(na.omit(sxi[,2:4])),3)
get_upper_tri <- function(test_multicor){
  test_multicor[lower.tri(test_multicor)]<- NA
  return(test_multicor)
}
test_multicor <- get_upper_tri(test_multicor)
test_multicor <- reshape2::melt(test_multicor, na.rm = TRUE)
ggplot(data = test_multicor, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1),
        axis.title = element_blank())+
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  coord_fixed()

ggplot(sxi) + 
  geom_point(aes(x = sri, y = sfi))
ggplot(sxi) + 
  geom_point(aes(x = sai, y = sfi))
ggplot(sxi) + 
  geom_point(aes(x = sri, y = sai))



### CWM community weighted mean by abundance and CM community mean for occurrence

cwm_cm <- ddply(bird_data, .(siteID,year), .fun = 
                    function(x){
                      ex_cwm <- merge(x,sxi, by = "sci_name_out",all.x=TRUE)
                      ex_cwm$sum_sfi <- ex_cwm$count*ex_cwm$sfi
                      ex_cwm$sum_sai <- ex_cwm$count*ex_cwm$sai
                      ex_cwm$sum_sri <- ex_cwm$count*ex_cwm$sri
                      
                      return(data.frame(cwm_sfi = sum(ex_cwm$sum_sfi)/sum(ex_cwm$count),
                                        cwm_sai = sum(ex_cwm$sum_sai)/sum(ex_cwm$count),
                                        cwm_sri = sum(ex_cwm$sum_sri)/sum(ex_cwm$count),
                                        cm_sfi = mean(ex_cwm$sfi),
                                        cm_sai = mean(ex_cwm$sai),
                                        cm_sri = mean(ex_cwm$sri)
                                        ))

                    },
                  .progress = "text"
) 


# compare cwm vs cm


ggplot(cwm_cm[which(cwm_cm$year==2011),]) + 
  geom_point(aes(x = cwm_sri, y = cm_sri))
ggplot(cwm_cm[which(cwm_cm$year==2011),]) + 
  geom_point(aes(x = cwm_sai, y = cm_sai))
ggplot(cwm_cm[which(cwm_cm$year==2011),]) + 
  geom_point(aes(x = cwm_sfi, y = cm_sfi))

ggplot(cwm_cm[which(cwm_cm$year==2011),]) + 
  geom_point(aes(x = cm_sri, y = cm_sai))
ggplot(cwm_cm[which(cwm_cm$year==2011),]) + 
  geom_point(aes(x = cm_sfi, y = cm_sai))
ggplot(cwm_cm[which(cwm_cm$year==2011),]) + 
  geom_point(aes(x = cm_sfi, y = cm_sri))
