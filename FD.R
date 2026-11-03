library(ade4)
library(RVAideMemoire)
library(cluster)
library(clue)
library(sf)
library(ggplot2)
library(plyr)
library(dplyr)
library(terra)
library(exactextractr)
library(mgcv)
library(gratia)

# Load bird data

bird_data <- readRDS("data/bird_data_clean.rds")

### SFI

# load trait-species matrix
trait <- read.csv("data/life_history_bird2.csv", header = TRUE)

# species name consistency
trait$Species <- gsub("_", " ", trait$Species)
trait$Species[which(
  trait$Species == "Picus viridis"
)] <- "Picus viridis / Picus sharpei"
trait$Species[which(
  trait$Species == "Sylvia cantillans"
)] <- "Sylvia cantillans / Sylvia subalpina"
trait$Species[which(
  trait$Species == "Dendrocoptes medius"
)] <- "Leiopicus medius"
trait$Species[which(trait$Species == "Luscinia svecica")] <- "Cyanecula svecica"
trait$Species[which(trait$Species == "Cyanopica cyanus")] <- "Cyanopica cooki"

# select the same species as for the other indices
trait <- trait[which(trait$Species %in% unique(bird_data$sci_name_out)), ]

# remove taxonomic info
trait <- trait[, -c(1:3)]

# remove traits not well covered
table.na <- apply(trait, 2, function(x) {
  sum(is.na(x))
})
trait2 <- trait[, -which(table.na > 12)]

# select traits
# all traits that could be relevant
trait_selected <- c(
  "Species",
  "LengthU_MEAN",
  "WingU_MEAN",
  "TailU_MEAN",
  "BillU_MEAN",
  "TarsusU_MEAN",
  "WeightU_MEAN",
  "Sexual.dimorphism",
  "Clutch_MEAN",
  "Broods.per.year",
  "EggL_MEAN",
  "EggW_MEAN",
  "Egg_MASS",
  "Young",
  "Association.during.nesting",
  "Nest.type",
  "Nest.building",
  "Mating.system",
  "Incubation.period",
  "Incubation.sex",
  "Territoriality",
  "Sedentary",
  "Short.distance.migrant",
  "Long.distance.migrant",
  "Deciduous.forest",
  "Coniferous.forest",
  "Woodland",
  "Shrub",
  "Savanna",
  "Tundra",
  "Grassland",
  "Mountain.meadows",
  "Reed",
  "Swamps",
  "Desert",
  "Freshwater",
  "Marine",
  "Rocks",
  "Human.settlements",
  "Folivore_B",
  "Frugivore_B",
  "Granivore_B",
  "Arthropods_B",
  "Other.invertebrates_B",
  "Fish_B",
  "Other.vertebrates_B",
  "Carrion_B",
  "Omnivore_B"
)
# the one used in Carmona et al. 2021 https://doi.org/10.1126/sciadv.abf2675
trait_selected <- c(
  "Species",
  "Clutch_MEAN",
  "Broods.per.year",
  "WeightU_MEAN",
  "Incubation.period",
  "Age.of.first.breeding",
  "Fledging.period",
  "Egg_MASS",
  "LengthU_MEAN"
)
# the same but with diet as in Godet et al.
trait_selected <- c(
  "Species",
  "Clutch_MEAN",
  "Broods.per.year",
  "WeightU_MEAN",
  "Incubation.period",
  "Age.of.first.breeding",
  "Fledging.period",
  "Egg_MASS",
  "LengthU_MEAN",
  "Frugivore_B",
  "Granivore_B",
  "Arthropods_B",
  "Other.invertebrates_B",
  "Fish_B",
  "Other.vertebrates_B",
  "Carrion_B",
  "Omnivore_B"
)

trait3 <- trait2[, trait_selected]

# specify non noncontinuous traits
nc_traits <- c(
  "Sexual.dimorphism",
  "Young",
  "Association.during.nesting",
  "Nest.type",
  "Nest.building",
  "Mating.system",
  "Incubation.sex",
  "Territoriality",
  "Sedentary",
  "Short.distance.migrant",
  "Long.distance.migrant",
  "Deciduous.forest",
  "Coniferous.forest",
  "Woodland",
  "Shrub",
  "Savanna",
  "Tundra",
  "Grassland",
  "Mountain.meadows",
  "Reed",
  "Swamps",
  "Desert",
  "Freshwater",
  "Marine",
  "Rocks",
  "Human.settlements",
  "Folivore_B",
  "Frugivore_B",
  "Granivore_B",
  "Arthropods_B",
  "Other.invertebrates_B",
  "Fish_B",
  "Other.vertebrates_B",
  "Carrion_B",
  "Omnivore_B"
)
nc_trait_selected <- trait_selected[which(trait_selected %in% nc_traits)]

trait3[, nc_trait_selected] <-
  lapply(trait3[, nc_trait_selected], factor)

# transformation in distance matrix (gower)
asymm_trait <- c(
  "Territoriality",
  "Sedentary",
  "Short.distance.migrant",
  "Long.distance.migrant",
  "Deciduous.forest",
  "Coniferous.forest",
  "Woodland",
  "Shrub",
  "Savanna",
  "Tundra",
  "Grassland",
  "Mountain.meadows",
  "Reed",
  "Swamps",
  "Desert",
  "Freshwater",
  "Marine",
  "Rocks",
  "Human.settlements",
  "Folivore_B",
  "Frugivore_B",
  "Granivore_B",
  "Arthropods_B",
  "Other.invertebrates_B",
  "Fish_B",
  "Other.vertebrates_B",
  "Carrion_B",
  "Omnivore_B"
)
asymm_trait_selected <- which(trait3 %in% asymm_trait) - 1

if (length(asymm_trait_selected) == 0) {
  mat.dist <- daisy(trait3[, 2:ncol(trait3)], metric = "gower")
} else {
  mat.dist <- daisy(
    trait3[, 2:ncol(trait3)],
    metric = "gower",
    type = list(asymm = asymm_trait_selected)
  )
}
attr(mat.dist, "Labels") <- trait3$Species

data_sfi <- data.frame(species = trait3$Species, sfi = apply(mat.dist, 1, mean))


### SAI => should be rename SOI species occurrence index?

data_sai <- ddply(bird_data, .(sci_name_out), .fun = function(x) {
  # sites where the species occurs
  site_occ <- unique(x$siteID[which(x$count > 0)])
  # focus on occurrence frequency within the sites where the species occurs
  x_occ <- x[which(x$siteID %in% site_occ), ]
  # get occurrence instead of abundance
  x_occ$occ <- ifelse(x_occ$count > 0, 1, 0)
  # compute the frequency of occurence per site
  x_sai <- data.frame(
    x_occ |> dplyr::group_by(siteID) |> dplyr::summarise(mean_oc = mean(occ))
  )
  # return the mean of occurence frequency
  return(data.frame(sai = mean(x_sai$mean_oc)))
})

### SRI, calcuated using atlas in Godet et al. 2015, here using convex hull (minimum polygon encompassing all sites of occurence)

# site location
site_mainland_sf_reproj <- readRDS(
  "data/site_mainland_sf_reproj.rds"
)

# map of European countries involved in bird monitoring
grid_eu_mainland_outline <- st_read(
  "data/grid_eu_mainland_outline.gpkg"
)

# get the max convex hull (i.e. mimic the case when a species was present in all site, to serve as a reference)

bird_convh <- st_convex_hull(st_union(site_mainland_sf_reproj))
bird_convh_inter <- st_intersection(bird_convh, grid_eu_mainland_outline)

ggplot(bird_convh_inter) + geom_sf()

# a few example of species with different distribution, to check the approach
TURMER_convh <- st_convex_hull(st_union(site_mainland_sf_reproj[
  which(
    site_mainland_sf_reproj$siteID %in%
      bird_data$siteID[which(
        bird_data$sci_name_out == "Turdus merula" & bird_data$count > 0
      )]
  ),
]))
TURMER_convh_inter <- st_intersection(TURMER_convh, grid_eu_mainland_outline)
as.numeric(st_area(TURMER_convh_inter)) / as.numeric(st_area(bird_convh_inter))
ggplot(TURMER_convh_inter) + geom_sf()

SYLMEL_convh <- st_convex_hull(st_union(site_mainland_sf_reproj[
  which(
    site_mainland_sf_reproj$siteID %in%
      bird_data$siteID[which(
        bird_data$sci_name_out == "Sylvia melanocephala" & bird_data$count > 0
      )]
  ),
]))
SYLMEL_convh_inter <- st_intersection(SYLMEL_convh, grid_eu_mainland_outline)
as.numeric(st_area(SYLMEL_convh_inter)) / as.numeric(st_area(bird_convh_inter))
ggplot(SYLMEL_convh_inter) + geom_sf()

CYACOO_convh <- st_convex_hull(st_union(site_mainland_sf_reproj[
  which(
    site_mainland_sf_reproj$siteID %in%
      bird_data$siteID[which(
        bird_data$sci_name_out == "Cyanopica cooki" & bird_data$count > 0
      )]
  ),
]))
CYACOO_convh_inter <- st_intersection(CYACOO_convh, grid_eu_mainland_outline)
as.numeric(st_area(CYACOO_convh_inter)) / as.numeric(st_area(bird_convh_inter))
ggplot(CYACOO_convh_inter) + geom_sf()

# calculate for all species
data_sri <- ddply(
  bird_data,
  .(sci_name_out),
  .fun = function(x) {
    sp_convh <- st_convex_hull(st_union(site_mainland_sf_reproj[
      which(site_mainland_sf_reproj$siteID %in% x$siteID[which(x$count > 0)]),
    ]))
    sp_convh_inter <- st_intersection(sp_convh, grid_eu_mainland_outline)
    sp_sri <- as.numeric(st_area(sp_convh_inter)) /
      as.numeric(st_area(bird_convh_inter))
    return(data.frame(sri = sp_sri))
  },
  .progress = "text"
)


### SGIc, not computed as it is based on species habitat which is not relevant for other taxa

### merge the three indices
sxi <- merge(
  data_sri,
  data_sai[, c("sci_name_out", "sai")],
  by = "sci_name_out"
)
sxi <- merge(
  sxi,
  data_sfi,
  by.x = "sci_name_out",
  by.y = "species",
  all.x = TRUE
)

### check correlation
test_multicor <- round(cor(na.omit(sxi[, 2:4])), 3)
get_upper_tri <- function(test_multicor) {
  test_multicor[lower.tri(test_multicor)] <- NA
  return(test_multicor)
}
test_multicor <- get_upper_tri(test_multicor)
test_multicor <- reshape2::melt(test_multicor, na.rm = TRUE)
ggplot(data = test_multicor, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1, 1),
    space = "Lab",
    name = "Pearson\nCorrelation"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
    axis.title = element_blank()
  ) +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  coord_fixed()

ggplot(sxi) +
  geom_point(aes(x = sri, y = sfi))
ggplot(sxi) +
  geom_point(aes(x = sai, y = sfi))
ggplot(sxi) +
  geom_point(aes(x = sri, y = sai))

# write.csv(sxi,"outputs/sxi_bird.csv", row.names = FALSE)


### CWM community weighted mean by abundance and CM community mean for occurrence

cwm_cm <- ddply(
  bird_data,
  .(siteID, year),
  .fun = function(x) {
    ex_cwm <- merge(x, sxi, by = "sci_name_out", all.x = TRUE)
    ex_cwm$sum_sfi <- ex_cwm$count * ex_cwm$sfi
    ex_cwm$sum_sai <- ex_cwm$count * ex_cwm$sai
    ex_cwm$sum_sri <- ex_cwm$count * ex_cwm$sri

    return(data.frame(
      cwm_sfi = sum(ex_cwm$sum_sfi) / sum(ex_cwm$count),
      cwm_sai = sum(ex_cwm$sum_sai) / sum(ex_cwm$count),
      cwm_sri = sum(ex_cwm$sum_sri) / sum(ex_cwm$count),
      cm_sfi = mean(ex_cwm$sfi),
      cm_sai = mean(ex_cwm$sai),
      cm_sri = mean(ex_cwm$sri)
    ))
  },
  .progress = "text"
)


# compare cwm vs cm

ggplot(cwm_cm[which(cwm_cm$year == 2011), ]) +
  geom_point(aes(x = cwm_sri, y = cm_sri))
ggplot(cwm_cm[which(cwm_cm$year == 2011), ]) +
  geom_point(aes(x = cwm_sai, y = cm_sai))
ggplot(cwm_cm[which(cwm_cm$year == 2011), ]) +
  geom_point(aes(x = cwm_sfi, y = cm_sfi))

ggplot(cwm_cm[which(cwm_cm$year == 2011), ]) +
  geom_point(aes(x = cm_sri, y = cm_sai))
ggplot(cwm_cm[which(cwm_cm$year == 2011), ]) +
  geom_point(aes(x = cm_sfi, y = cm_sai))
ggplot(cwm_cm[which(cwm_cm$year == 2011), ]) +
  geom_point(aes(x = cm_sfi, y = cm_sri))

# write.csv(cwm_cm,"outputs/community_indices_bird.csv", row.names = FALSE)

### Explanatory variables

#### HFI

# load HFI data from Mu et al. 2022

hfi2018 <- rast("data/hfp2018.tif")

# crop and reproject hfi for Europe

grid_eu_mainland_outline_moll <- st_transform(grid_eu_mainland_outline, st_crs(hfi2018))

hfi2018_europe <- mask(hfi2018,grid_eu_mainland_outline_moll) # terrestrial cells

hfi2018_europe <- crop(hfi2018_europe,ext(grid_eu_mainland_outline_moll)) # remove oversea territories

hfi2018_europe_ETRS89 <- project(hfi2018_europe, "epsg:3035")

# select a buffer around each monitoring site

site_mainland_buffer <- st_buffer(site_mainland_sf_reproj, dist = 2500) # here a buffer with a 2.5km radius that can be changed

# get average hfi in each buffer

site_hfi <- exact_extract(hfi2018_europe_ETRS89,site_mainland_buffer, fun="mean")


#### Climate

# to download from https://www.chelsa-climate.org/datasets/chelsa_bioclim

# reproject and crop for Europe

for(i in c("01", # Mean annual temperature calculated as the average of mean monthly temperatures over the year
           "04", # Temperature seasonality given by the standard deviation of mean monthly temperatures
           "07", # Annual temperature range calculated as bio05 − bio06; measures amplitude between warmest and coldest months
           "12", # Sum of monthly precipitation totals across the year
           "15" # Coefficient of variation (100 × SD ÷ mean) of monthly precipitation totals
           )){
  
  past <- rast(paste0("data/CHELSA_bio",i,"_1981-2010_V.2.1.tif"))
  past <- terra::crop(past, grid_eu_mainland_outline_wgs84, snap="near", extend=FALSE)
  past <- terra::project(past,crs(grid_eu_mainland_outline))
  assign(paste0("bio",i,"_1981_2010_eu"),past)
  writeRaster(x = past, 
              filename = paste0("data/bio",i,"_1981_2010_eu.tif"),
              overwrite = TRUE)
}

# load the new climate maps

for(i in c("01","04","07","12","15")){
  past <- rast(paste0("data/bio",i,"_1981_2010_eu.tif"))
  past <- terra::mask(past,grid_eu_mainland_outline)
  assign(paste0("bio",i,"_1981_2010_eu"),past)
}

# extract climate data within each buffer

bio01_past <- exact_extract(bio01_1981_2010_eu,site_mainland_buffer, fun="mean") # Mean annual temperature calculated as the average of mean monthly temperatures over the year
bio04_past <- exact_extract(bio04_1981_2010_eu,site_mainland_buffer) # Temperature seasonality given by the standard deviation of mean monthly temperatures
bio04_past <- sqrt(laply(bio04_past, function(x){return(c(weighted.mean(x$value^2,w = x$coverage_fraction)))}))
bio07_past <- exact_extract(bio07_1981_2010_eu,site_mainland_buffer, fun="mean") # Annual temperature range calculated as bio05 − bio06; measures amplitude between warmest and coldest months
bio12_past <- exact_extract(bio12_1981_2010_eu,site_mainland_buffer, fun="mean") # Sum of monthly precipitation totals across the year
bio15_past <- exact_extract(bio15_1981_2010_eu,site_mainland_buffer, fun="mean") # Coefficient of variation (100 × SD ÷ mean) of monthly precipitation totals

### Analysis

#### Prepare the dataset

# combine response and explanatory variables

site_mainland_buffer$hfi <- site_hfi
site_mainland_buffer$bio01_past <- bio01_past
site_mainland_buffer$bio04_past <- bio04_past
site_mainland_buffer$bio07_past <- bio07_past
site_mainland_buffer$bio12_past <- bio12_past
site_mainland_buffer$bio15_past <- bio15_past

site_mainland_buffer_nosf <- site_mainland_buffer # get an non spatial object
st_geometry(site_mainland_buffer_nosf) <- NULL

cwm_cm_hfi <- merge(cwm_cm,site_mainland_buffer_nosf,by="siteID", all.x=TRUE)

# write.csv(cwm_cm_hfi,"outputs/community_indices_bird_with_covariates.csv", row.names = FALSE)

# check correlation (bio04 and bio07 very similar, keep just bio07 as more easy to interpret)

test_multicor <- round(cor(na.omit(site_mainland_buffer_nosf[,c("hfi","bio01_past","bio04_past","bio07_past","bio12_past","bio15_past")])), 3)
get_upper_tri <- function(test_multicor) {
  test_multicor[lower.tri(test_multicor)] <- NA
  return(test_multicor)
}
test_multicor <- get_upper_tri(test_multicor)
test_multicor <- reshape2::melt(test_multicor, na.rm = TRUE)
ggplot(data = test_multicor, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1, 1),
    space = "Lab",
    name = "Pearson\nCorrelation"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
    axis.title = element_blank()
  ) +
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) +
  coord_fixed()

#### Run the analysis

# get a first look over the relationship

ggplot(cwm_cm_hfi) +
  geom_point(aes(y = cm_sfi, x=hfi))
ggplot(cwm_cm_hfi) +
  geom_point(aes(y = cm_sai, x=hfi))
ggplot(cwm_cm_hfi) +
  geom_point(aes(y = cm_sri, x=hfi))

# check distribution of the response

ggplot(cwm_cm_hfi) +
  geom_histogram(aes(x = cm_sfi))
ggplot(cwm_cm_hfi) +
  geom_histogram(aes(x = cm_sai))
ggplot(cwm_cm_hfi) +
  geom_histogram(aes(x = cm_sri))

# select just one year per site to avoid providing too much weight to site sampled several years and get independence between values

cwm_cm_hfi_oneyearpersite <- ddply(
  cwm_cm_hfi,
  .(siteID),
  .fun = function(x) {
    set.seed(123) # for reproducibility
    x_sample <- sample_n(x, 1)
    
    return(x_sample)
  },
  .progress = "text"
)

# write.csv(cwm_cm_hfi_oneyearpersite,"outputs/bird_data_to_run_analysis.csv", row.names = FALSE)

##### SFI

formula_gam <- "cm_sfi ~ s(hfi,k=4) + s(bio01_past,k=4) + s(bio07_past,k=4) + s(bio12_past,k=4) + s(bio15_past,k=4)"

global_mod_sfi_gauss <- bam(as.formula(paste(formula_gam,sep=" + ",paste(c("area_sampled_m2", # area monitored, i.e. sampling effort, should explain part of the change across community
                                                                           "te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=10)" # tensor product to handle spatial autocorrelation
                                                                           ), collapse = " + "))),
                  data=cwm_cm_hfi_oneyearpersite)

summary(global_mod_sfi_gauss)
gam.check(global_mod_sfi_gauss) # seems that variance is stable but gaussian assumptium not the best choice, res vs. fitted is fine
plot(global_mod_sfi_gauss, residuals = TRUE) # almost linear increase (first plot) of sfi with hfi



global_mod_sfi_gamma <- bam(as.formula(paste(formula_gam,sep=" + ",paste(c("area_sampled_m2","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=10)"), collapse = " + "))),
                            family=Gamma(link=log),data=cwm_cm_hfi_oneyearpersite)

summary(global_mod_sfi_gamma)
gam.check(global_mod_sfi_gamma) # no big improvement compare to gaussian -> stay with gaussian to start as more easy to interpret
plot(global_mod_sfi_gamma, residuals = TRUE)

# try with polynomials for hfi to interpret more easily the relationship SFI vs. HFI

formula_gam <- "cm_sfi ~ hfi + I(hfi^2) + I(hfi^3) + s(bio01_past,k=4) + s(bio07_past,k=4) + s(bio12_past,k=4) + s(bio15_past,k=4)"

global_mod_sfi_gauss <- bam(as.formula(paste(formula_gam,sep=" + ",paste(c("area_sampled_m2","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=10)"), collapse = " + "))),
                            data=cwm_cm_hfi_oneyearpersite)

summary(global_mod_sfi_gauss)
gam.check(global_mod_sfi_gauss)
plot(global_mod_sfi_gauss, residuals = TRUE)

ds <- data_slice(global_mod_sfi_gauss, hfi = evenly(hfi, n = 100))
fv <- fitted_values(global_mod_sfi_gauss, data = ds)
ggplot(fv,aes(x = hfi, y = .fitted)) +
  geom_point(data = cwm_cm_hfi_oneyearpersite, mapping = aes(y = cm_sfi), size = 0.5) +
  geom_function(fun = function(x){summary(global_mod_sfi_gauss)$p.coeff[1]+x*summary(global_mod_sfi_gauss)$p.coeff[2]+x^2*summary(global_mod_sfi_gauss)$p.coeff[3]+x^3*summary(global_mod_sfi_gauss)$p.coeff[4]}, colour = "blue", xlim=c(0,50))

# there is a small decrease for very low of hfi, perhaps we can use exp or piecewise, to check


##### SAI

formula_gam <- "cm_sai ~ s(hfi,k=4) + s(bio01_past,k=4) + s(bio07_past,k=4) + s(bio12_past,k=4) + s(bio15_past,k=4)"

global_mod_sai_gauss <- bam(as.formula(paste(formula_gam,sep=" + ",paste(c("area_sampled_m2","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=10)"), collapse = " + "))),
                            data=cwm_cm_hfi_oneyearpersite)

summary(global_mod_sai_gauss)
gam.check(global_mod_sai_gauss)
plot(global_mod_sai_gauss, residuals = TRUE)


formula_gam <- "cm_sai ~ hfi + I(hfi^2) + I(hfi^3) + s(bio01_past,k=4) + s(bio07_past,k=4) + s(bio12_past,k=4) + s(bio15_past,k=4)"

global_mod_sai_gauss <- bam(as.formula(paste(formula_gam,sep=" + ",paste(c("area_sampled_m2","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=10)"), collapse = " + "))),
                            data=cwm_cm_hfi_oneyearpersite)

summary(global_mod_sai_gauss)
gam.check(global_mod_sai_gauss)
plot(global_mod_sai_gauss, residuals = TRUE)

ds <- data_slice(global_mod_sai_gauss, hfi = evenly(hfi, n = 100))
fv <- fitted_values(global_mod_sai_gauss, data = ds)
ggplot(fv,aes(x = hfi, y = .fitted)) +
  geom_point(data = cwm_cm_hfi_oneyearpersite, mapping = aes(y = cm_sai), size = 0.5) +
  geom_function(fun = function(x){summary(global_mod_sai_gauss)$p.coeff[1]+x*summary(global_mod_sai_gauss)$p.coeff[2]+x^2*summary(global_mod_sai_gauss)$p.coeff[3]+x^3*summary(global_mod_sai_gauss)$p.coeff[4]}, colour = "blue", xlim=c(0,50))

# better fit with s(hfi), here also exp could be checked instead of poly

formula_gam <- "cm_sri ~ s(hfi,k=4) + s(bio01_past,k=4) + s(bio07_past,k=4) + s(bio12_past,k=4) + s(bio15_past,k=4)"

global_mod_sri_gauss <- bam(as.formula(paste(formula_gam,sep=" + ",paste(c("area_sampled_m2","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=10)"), collapse = " + "))),
                            data=cwm_cm_hfi_oneyearpersite)

summary(global_mod_sri_gauss)
gam.check(global_mod_sri_gauss)
plot(global_mod_sri_gauss, residuals = TRUE)


formula_gam <- "cm_sri ~ hfi + I(hfi^2) + I(hfi^3) + s(bio01_past,k=4) + s(bio07_past,k=4) + s(bio12_past,k=4) + s(bio15_past,k=4)"

global_mod_sri_gauss <- bam(as.formula(paste(formula_gam,sep=" + ",paste(c("area_sampled_m2","te(Long_LAEA,Lat_LAEA,bs='tp',fx=TRUE,k=10)"), collapse = " + "))),
                            data=cwm_cm_hfi_oneyearpersite)

summary(global_mod_sri_gauss)
gam.check(global_mod_sri_gauss)
plot(global_mod_sri_gauss, residuals = TRUE)

ds <- data_slice(global_mod_sri_gauss, hfi = evenly(hfi, n = 100))
fv <- fitted_values(global_mod_sri_gauss, data = ds)
ggplot(fv,aes(x = hfi, y = .fitted)) +
  geom_point(data = cwm_cm_hfi_oneyearpersite, mapping = aes(y = cm_sri), size = 0.5) +
  geom_function(fun = function(x){summary(global_mod_sri_gauss)$p.coeff[1]+x*summary(global_mod_sri_gauss)$p.coeff[2]+x^2*summary(global_mod_sri_gauss)$p.coeff[3]+x^3*summary(global_mod_sri_gauss)$p.coeff[4]}, colour = "blue", xlim=c(0,50))

