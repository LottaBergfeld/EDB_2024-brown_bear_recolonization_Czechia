library(RangeShiftR)
library(terra)
library(viridis)
library(RColorBrewer)
library(paletteer)
library(grid)
library(gridExtra)
library(tidyverse)
library(ggplot2) 

setwd("C:/Users/lotta/OneDrive/Dokumente/CLEWS_Uni_Potsdam/4_Semester/Ecosystem_Dynamics_and_Biodiversity/Seminar_paper/R-Code")

source("make_our_patches.R")

dirpath = "C:/Users/lotta/OneDrive/Dokumente/CLEWS_Uni_Potsdam/4_Semester/Ecosystem_Dynamics_and_Biodiversity/Seminar_paper/R-Code/"

dir.create(paste0(dirpath,"Inputs"), showWarnings = TRUE)
dir.create(paste0(dirpath,"Outputs"), showWarnings = TRUE)
dir.create(paste0(dirpath,"Output_Maps"), showWarnings = TRUE)


#1. Landscape Module
# read in land cover data for Czechia
landsc <- terra::rast(paste0(dirpath, "Inputs/Corine_2018_4km_Czechia.asc"))
landsc.f <- as.factor(landsc)
# add the land cover classes to the raster attribute table
(rat <- levels(landsc.f)[[1]])
rat[["Corine_2018_4km_Czechia"]] <- c("Artificial surfaces", "Arable land", "Pastures and heterogeneous agricultural areas", "Forests", "Grasslands, moors and heathland", "Wetlands and water bodies")
levels(landsc.f) <- rat
# Plot
par(mfrow=c(1,1))
plot(landsc.f, axes = F, col=paletteer_d("wesanderson::Cavalcanti1"))

# with the make_our_patches script we made the forests raster, where 1 = forest
# Read in forest patch files
forest_patches <- terra::rast(paste0(dirpath,'Inputs/woodland_patchIDs_2018_4km_Czechia.asc'))
# Makes spatial polygons
forest_patches_pol <- as.polygons(forest_patches, dissolve=T)

# Plot forest patches and add labels
par(mfrow=c(1,1))
plot(forest_patches,axes=F, legend=F, col = c('grey',rep(brewer.pal(n = 10, name = "Paired"),4)))
forest_patches_pol <- as.polygons(forest_patches, dissolve=T) # Makes spatial polygons
text(forest_patches_pol,labels=values(forest_patches_pol)[,1])

#define landscape parameter
land <- ImportedLandscape(LandscapeFile = "Corine_2018_4km_Czechia.asc",
                          PatchFile = "woodland_patchIDs_2018_4km_Czechia.asc", 
                          Resolution = 4000,
                          Nhabitats = 7, # 7 types of habitats
                          K_or_DensDep = c(0, 0, 0, 0.0006, 0, 0,0)
)


# 2. Demography Module
# Define Transition matrix (in our case it is rather an age matrix than a stage transition matrix)
# taken from: https://compadre-db.org/Species/46451
trans_mat <- matrix(c(0, 0, 0, 0.219551, 0.219551, 0.219551, 0.219551, 0.219551, 0.304977, 0.827, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.866, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.919, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.911, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.911, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.911, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.911, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.911, 0.831 ), nrow = 9, byrow = T) 
# stages: yearling, subadult 2 years old, subadult 3 years old, adult 4 years old, adult 5 years old, adult 6 years old, adult 7 years old, adult 8 years old, adult 9+ years old

stg <- StageStructure(Stages = 9,
                      TransMatrix = trans_mat,
                      MaxAge = 35,
                      FecDensDep = F,
                      SurvDensDep = F, # mortality doesn't increase with density
                      SurvDensCoeff = 1 # Density dependence relative to 1/b
)

demo <- Demography(StageStruct = stg, ReproductionType = 0) # 0=all female


#3. Dispersal
# Define dispersal module / kernel
disp <-  Dispersal(
  # Emigration phase: stage 1 and 2 have constant emigration probability of 0.2
  Emigration = Emigration(StageDep=T, EmigProb = cbind(0:8,c(0, 0.2, 0.2, 0, 0, 0, 0, 0, 0)) ),
  # Transfer phase: negative exponential dispersal kernel with mean dispersal distance of 24km
  Transfer = DispersalKernel(Distances = 24000),
  # Settlement: if individual arrives in unsuitable cells, it can randomly chose a suitable neighbouring cell or will die
  Settlement = Settlement(Settle = 2)) 


#4. Initialisation
# single site-introduction because they will all be in patch 26 at the beginning
# put 2 female bears of age 4 (grown-ups)
# prepare dataframe for InitIndsFile
(init_df <- data.frame(Year=0,Species=0,PatchID=26,Ninds=c(2),Age=4,Stage=4))

# write InitIndsFile to file
write.table(init_df, file=paste0(dirpath,'Inputs/InitInds.txt'), sep='\t', row.names=F, quote=F)
# Set initialisation
init_pop <- Initialise(InitType = 2,       # Initialise from initial individuals list file
                       InitIndsFile = 'InitInds.txt')


#5. Simulation setting
# Set the number of replicate simulations to run 
RepNb <- 100
sim_years <- 50

sim <- Simulation(Simulation = 0, # ID
                  Replicates = RepNb, # number of replicate runs
                  Years = sim_years, # number of simulated years
                  OutIntPop = 1, # output interval
                  OutIntOcc = 1,
                  OutIntRange = 1)


#6. Parameter Master
s <- RSsim(batchnum = 1, simul = sim, land = land, demog = demo, dispersal = disp, init = init_pop, seed=324139)


#7. Simulation Run
RunRS(s,dirpath)


#Visualization and Interpretation

# plot the resulting abundances and occupancy
# general output of population size + occupancy
par(mfrow=c(1,2))
plotAbundance(s,dirpath,sd=T, rep=F)
plotOccupancy(s, dirpath, sd=T, rep=F)

#code for actual map
# Colonisation metrics
col_stats <- ColonisationStats(s, dirpath, years = 50, maps = T)
# mean occupancy probability in year 100
head(col_stats$occ_prob)

# map occupancy probability
par(mfrow=c(1,1))
mycol_occprob <- colorRampPalette(c('gold', 'darkgreen', 'black'))
plot(new.raster, axes=F, legend=F, col=c('lightgrey'))
plot(col_stats$map_occ_prob, add=TRUE, plg=list( title="Occupancy\nProbability", title.cex=0.9),
     pax=list(side=1:4, retro=TRUE), axes=F, range=c(0,1), col=mycol_occprob(10), type="continuous")
text(forest_patches_pol,labels=values(forest_patches_pol)[,1])

# map colonisation time on landscape background
par(mfrow=c(1,1))
mycol_coltime <- colorRampPalette(c('ivory2', 'gold', 'olivedrab3', 'black'))
plot(new.raster, axes=F, legend=F, col=c('lightgrey'))
plot(col_stats$map_col_time, add=TRUE, plg=list( title="Colonisation\nTime", title.cex=0.9), 
     axes=F, breaks=c(-9,seq(-9,50,length=11)), col=mycol_coltime(6), 
     type="continuous")
text(forest_patches_pol,labels=values(forest_patches_pol)[,1])

# read population output file into a data frame
pop <- readPop(s, dirpath)


# Calculate survival probability as number of replicate with surviving individuals per year
# Extinction probability is 1- survival probability:

# Calculate  extinction probability by hand:
pop %>%
  group_by(Rep,Year) %>%
  # Sum individuals over all cells per year and replicate
  summarise(sumPop = sum(NInd), .groups='keep') %>%
  group_by(Year) %>%
  # Average extinction probability (1 minus the proportion of replicates with surviving populations)
  summarise(extProb = 1-sum(sumPop>0, na.rm=T)/RepNb)

# Mean time to extinction:
pop %>%
  group_by(Rep,Year) %>%
  # Sum individuals over all cells per year and replicate    
  summarise(sumPop = sum(NInd), .groups='keep') %>% 
  # Identify in which year they go extinct
  filter(sumPop==0) %>% 
  pull(Year) %>% mean

# Define a function for calculating extinction probability
Calc_ExtProb <- function(pop_df,s) {
  require(dplyr)
  
  pop_df %>%
    group_by(Rep,Year) %>%
    # Sum individuals over all cells per year and replicate
    summarise(sumPop = sum(NInd), .groups='keep') %>%
    group_by(Year) %>%
    # Average extinction probability (1 minus the proportion of replicates with surviving populations)
    summarise(extProb = 1-sum(sumPop>0, na.rm=T)/RepNb) %>%
    # Make sure that data frame is filled until last year of simulation
    right_join(tibble(Year = seq_len(s@simul@Years)), by='Year') %>% replace_na(list(extProb=1))
}

# Define a function for calculating mean time to extinction
Calc_ExtTime <- function(pop_df) {
  require(dplyr)
  
  pop_df %>%
    group_by(Rep,Year) %>%
    # Sum individuals over all cells per year and replicate    
    summarise(sumPop = sum(NInd), .groups='keep') %>% 
    # Identify in which year they go extinct
    filter(sumPop==0) %>% 
    pull(Year) %>% mean
}

# extinction probability
extProb <- Calc_ExtProb(pop,s)

# Plot extinction probabilities
ggplot(data = extProb, mapping = aes(x = Year, y = extProb)) + 
  labs(y = "Extinction Probability") +
  geom_line(linewidth=1, color = "darkgreen") +
  ylim(0,1)

# mean time to extinction
Calc_ExtTime(pop)

#--------------------------------------------
#---Sensitivity Analysis---------------------
#--------------------------------------------

#-NUMBER OF BEARS----------------------------

# what happens if 5 or 10 female bears cross the border from Slovakia to Czechia? 
# single site-introduction 
# put 5 or 10 female bears of age 4 (grown-ups)
# prepare dataframe for InitIndsFile
(init_df_10bears <- data.frame(Year=0,Species=0,PatchID=26,Ninds=c(10),Age=4,Stage=4))
(init_df_5bears <- data.frame(Year=0,Species=0,PatchID=26,Ninds=c(5),Age=4,Stage=4))

# write InitIndsFile to file
write.table(init_df_10bears, file=paste0(dirpath,'Inputs/InitInds_10bears.txt'), sep='\t', row.names=F, quote=F)
write.table(init_df_5bears, file=paste0(dirpath,'Inputs/InitInds_5bears.txt'), sep='\t', row.names=F, quote=F)

# Set initialisation
init_pop_10bears <- Initialise(InitType = 2,       # Initialise from initial individuals list file
                               InitIndsFile = 'InitInds_10bears.txt')
init_pop_5bears <- Initialise(InitType = 2,       # Initialise from initial individuals list file
                              InitIndsFile = 'InitInds_5bears.txt')

#6. Parameter Master
s_10bears <- RSsim(batchnum = 4, simul = sim, land = land, demog = demo, dispersal = disp, init = init_pop_10bears, seed = 324139)
s_5bears <- RSsim(batchnum = 5, simul = sim, land = land, demog = demo, dispersal = disp, init = init_pop_5bears, seed = 324139)

#check for parameter conflicts
validateRSparams(s_10bears)
validateRSparams(s_5bears)

#7. Simulation Run
RunRS(s_10bears,dirpath)
RunRS(s_5bears,dirpath)

#Visualization and Interpretation

# plot the resulting abundances and occupancy
par(mfrow=c(1,3))
# Abundances
plotAbundance(s, dirpath,sd=T, rep=F, main="2 bears")
plotAbundance(s_5bears, dirpath,sd=T, rep=F, main="5 bears")
plotAbundance(s_10bears, dirpath,sd=T, rep=F, main="10 bears")

par(mfrow=c(1,3))
# Occupancy
plotOccupancy(s, dirpath,sd=T, rep=F, main="2 bears")
plotOccupancy(s_5bears, dirpath,sd=T, rep=F, main="5 bears")
plotOccupancy(s_10bears, dirpath,sd=T, rep=F, main="10 bears")

#code for actual map
# Colonisation metrics
col_stats_5bears <- ColonisationStats(s_5bears, dirpath, years = 50, maps = T)
# mean occupancy probability in year 100
head(col_stats_5bears$occ_prob)

# read population output file into a data frame
pop_5bears <- readPop(s_5bears, dirpath)

# Colonisation metrics
col_stats_10bears <- ColonisationStats(s_10bears, dirpath, years = 50, maps = T)
# mean occupancy probability in year 100
head(col_stats_10bears$occ_prob)

# read population output file into a data frame
pop_10bears <- readPop(s_10bears, dirpath)


# Calculate survival probability as number of replicate with surviving individuals per year
# Extinction probability is 1- survival probability:

# Calculate  extinction probability by hand:
pop_5bears %>%
  group_by(Rep,Year) %>%
  # Sum individuals over all cells per year and replicate
  summarise(sumPop = sum(NInd), .groups='keep') %>%
  group_by(Year) %>%
  # Average extinction probability (1 minus the proportion of replicates with surviving populations)
  summarise(extProb = 1-sum(sumPop>0, na.rm=T)/RepNb)

# Mean time to extinction:
pop_5bears %>%
  group_by(Rep,Year) %>%
  # Sum individuals over all cells per year and replicate    
  summarise(sumPop = sum(NInd), .groups='keep') %>% 
  # Identify in which year they go extinct
  filter(sumPop==0) %>% 
  pull(Year) %>% mean

# Calculate  extinction probability by hand:
pop_10bears %>%
  group_by(Rep,Year) %>%
  # Sum individuals over all cells per year and replicate
  summarise(sumPop = sum(NInd), .groups='keep') %>%
  group_by(Year) %>%
  # Average extinction probability (1 minus the proportion of replicates with surviving populations)
  summarise(extProb = 1-sum(sumPop>0, na.rm=T)/RepNb)

# Mean time to extinction:
pop_10bears %>%
  group_by(Rep,Year) %>%
  # Sum individuals over all cells per year and replicate    
  summarise(sumPop = sum(NInd), .groups='keep') %>% 
  # Identify in which year they go extinct
  filter(sumPop==0) %>% 
  pull(Year) %>% mean

# extinction probability
extProb_5bears <- Calc_ExtProb(pop_5bears,s_5bears)
extProb_10bears <- Calc_ExtProb(pop_10bears,s_10bears)

# mean time to extinction
Calc_ExtTime(pop_10bears)


#-SETTLEMENT PARAMETER-----------------------

disp_settle_3 <-  Dispersal(
  # Emigration phase: stage 1 and 2 have constant emigration probability of 0.2
  Emigration = Emigration(StageDep=T, EmigProb = cbind(0:8,c(0, 0.2, 0.2, 0, 0, 0, 0, 0, 0)) ),
  # Transfer phase: negative exponential dispersal kernel with mean dispersal distance of 24km
  Transfer = DispersalKernel(Distances = 24000),
  # Settlement: if individual arrives in unsuitable cells, it can survive a time step before further dispersing
  Settlement = Settlement(Settle = 3)) 
s_settle_3 <- RSsim(batchnum = 4, simul = sim, land = land, demog = demo, dispersal = disp_settle_3, init = init_pop, seed = 324139)
RunRS(s_settle_3,dirpath)

# Plots
plotAbundance(s_settle_3, dirpath,sd=T, rep=F, main="Settlement Parameter 3")
plotOccupancy(s_settle_3, dirpath,sd=T, rep=F, main="Settlement Parameter 3")


# Colonisation metrics
col_stats_set_par_3 <- ColonisationStats(s_settle_3, dirpath, years = 50, maps = T)
# mean occupancy probability in year 100
head(col_stats_set_par_3$occ_prob)

# read population output file into a data frame
pop_set_par_3 <- readPop(s_settle_3, dirpath)


# Calculate survival and extinction probability 
pop_set_par_3 %>%
  group_by(Rep,Year) %>%
  # Sum individuals over all cells per year and replicate
  summarise(sumPop = sum(NInd), .groups='keep') %>%
  group_by(Year) %>%
  # Average extinction probability (1 minus the proportion of replicates with surviving populations)
  summarise(extProb = 1-sum(sumPop>0, na.rm=T)/RepNb)

# Mean time to extinction:
pop_set_par_3 %>%
  group_by(Rep,Year) %>%
  # Sum individuals over all cells per year and replicate    
  summarise(sumPop = sum(NInd), .groups='keep') %>% 
  # Identify in which year they go extinct
  filter(sumPop==0) %>% 
  pull(Year) %>% mean

# extinction probability
extProb_set_par_3 <- Calc_ExtProb(pop_set_par_3,s_settle_3)

# mean time to extinction
Calc_ExtTime(pop_set_par_3)


#-DENSITY DEPENDENCE-------------------------

# We set the Survival Density dependence Parameter to TRUE, which means that mortality increases with density
stg_T <- StageStructure(Stages = 9,
                        TransMatrix = trans_mat,
                        MaxAge = 35,
                        FecDensDep = F,
                        SurvDensDep = T, # mortality does increase with density
                        SurvDensCoeff = 1 # Density dependence relative to 1/b
)
demo_stg_T <- Demography(StageStruct = stg_T, ReproductionType = 0) # 0=all female

s_stg_T <- RSsim(batchnum = 4, simul = sim, land = land, demog = demo_stg_T, dispersal = disp, init = init_pop, seed = 324139)
RunRS(s_stg_T,dirpath)

# Plots
plotAbundance(s_stg_T, dirpath,sd=T, rep=F, main="Survival Density Dependence")
plotOccupancy(s_stg_T, dirpath,sd=T, rep=F, main="Survival Density Dependence")

# Colonisation metrics
col_stats_stg_T <- ColonisationStats(s_stg_T, dirpath, years = 50, maps = T)

# read population output file into a data frame
pop_stg_T <- readPop(s_stg_T, dirpath)


# Calculate survival and extinction probability 
pop_stg_T %>%
  group_by(Rep,Year) %>%
  # Sum individuals over all cells per year and replicate
  summarise(sumPop = sum(NInd), .groups='keep') %>%
  group_by(Year) %>%
  # Average extinction probability (1 minus the proportion of replicates with surviving populations)
  summarise(extProb = 1-sum(sumPop>0, na.rm=T)/RepNb)

# Mean time to extinction:
pop_stg_T %>%
  group_by(Rep,Year) %>%
  # Sum individuals over all cells per year and replicate    
  summarise(sumPop = sum(NInd), .groups='keep') %>% 
  # Identify in which year they go extinct
  filter(sumPop==0) %>% 
  pull(Year) %>% mean

# extinction probability
extProb_stg_T <- Calc_ExtProb(pop_stg_T,s_stg_T)

# mean time to extinction
Calc_ExtTime(pop_stg_T)


#-SUITABLE LAND COVER TYPES------------------

#1. Landscape Module
# this time we use the raster where forest and grassland = 1
# Read in patch files
forest_and_grassland_patches <- terra::rast(paste0(dirpath,'Inputs/wood_and_grassland_patchIDs_2018_4km_Czechia.asc'))
# Makes spatial polygons
forest_and_grassland_patches_pol <- as.polygons(forest_and_grassland_patches, dissolve=T)

# Plot forest patches and add labels
par(mfrow=c(1,1))
plot(forest_and_grassland_patches,axes=F, legend=F, col = c('grey',rep(brewer.pal(n = 10, name = "Paired"),4)))
forest_and_grassland_patches_pol <- as.polygons(forest_and_grassland_patches, dissolve=T) # Makes spatial polygons
text(forest_and_grassland_patches_pol,labels=values(forest_and_grassland_patches_pol)[,1])

#define landscape parameter
land_with_grassland <- ImportedLandscape(LandscapeFile = "Corine_2018_4km_Czechia.asc",
                          PatchFile = "wood_and_grassland_patchIDs_2018_4km_Czechia.asc", 
                          Resolution = 4000,
                          Nhabitats = 7, # 7 types of habitats
                          K_or_DensDep = c(0, 0, 0, 0.0006, 0, 0,0)
)

#4. Initialisation
# different patches and IDs, so now our initial population is in patch 28
# single site-introduction because they will all be in patch 26 at the beginning
# put 2 female bears of age 4 (grown-ups)
# prepare dataframe for InitIndsFile
(init_df_with_grassland <- data.frame(Year=0,Species=0,PatchID=28,Ninds=c(2),Age=4,Stage=4))

# write InitIndsFile to file
write.table(init_df_with_grassland, file=paste0(dirpath,'Inputs/InitInds_with_grassland.txt'), sep='\t', row.names=F, quote=F)
# Set initialisation
init_pop_with_grassland <- Initialise(InitType = 2,       # Initialise from initial individuals list file
                       InitIndsFile = 'InitInds_with_grassland.txt')

#6. Parameter Master
s_with_grassland <- RSsim(batchnum = 1, simul = sim, land = land_with_grassland, demog = demo, dispersal = disp, init = init_pop_with_grassland, seed=324139)


#7. Simulation Run
RunRS(s_with_grassland,dirpath)

# Plots
plotAbundance(s_with_grassland, dirpath,sd=T, rep=F, main="Survival Density Dependence")
plotOccupancy(s_with_grassland, dirpath,sd=T, rep=F, main="Survival Density Dependence")

# Colonisation metrics
col_stats_with_grassland <- ColonisationStats(s_with_grassland, dirpath, years = 50, maps = T)

# read population output file into a data frame
pop_with_grassland <- readPop(s_with_grassland, dirpath)


# Calculate survival and extinction probability 
pop_with_grassland %>%
  group_by(Rep,Year) %>%
  # Sum individuals over all cells per year and replicate
  summarise(sumPop = sum(NInd), .groups='keep') %>%
  group_by(Year) %>%
  # Average extinction probability (1 minus the proportion of replicates with surviving populations)
  summarise(extProb = 1-sum(sumPop>0, na.rm=T)/RepNb)

# Mean time to extinction:
pop_with_grassland %>%
  group_by(Rep,Year) %>%
  # Sum individuals over all cells per year and replicate    
  summarise(sumPop = sum(NInd), .groups='keep') %>% 
  # Identify in which year they go extinct
  filter(sumPop==0) %>% 
  pull(Year) %>% mean

# extinction probability
extProb_with_grassland <- Calc_ExtProb(pop_with_grassland,s_with_grassland)

# mean time to extinction
Calc_ExtTime(pop_with_grassland)


#-COMBINING PLOTS----------------------------

# 1. Combining all data in one data frame for put all in 1 Plot such as the Extinction Prob
data <- data.frame(extProb, extProb_5bears, extProb_10bears, extProb_set_par_3, extProb_stg_T, extProb_with_grassland)

# Have 6 szenarios
# Extinction Probability for all Sensitivity Analysis Things.

pdf(file = "Output_Maps/ext_prob.pdf", width = 7, height = 4) 

par(mfrow=c(1,1))
ggplot(data = data) + 
  labs(y = "Extinction Probability", color = "Simulation run") +
  geom_line(aes(x = Year, y = extProb, color = "2 Bears"), linewidth=1) +
  geom_line(aes(x = Year, y = extProb.1, color = "5 Bears"), linewidth=1) +
  geom_line(aes(x = Year, y = extProb.2, color = "10 Bears"), linewidth=1) +
  geom_line(aes(x = Year, y = extProb.3, color = "Settlement par. 3"), linewidth=1)+
  geom_line(aes(x = Year, y = extProb.4, color = "Survival dens. dep."), linewidth=1)+
  geom_line(aes(x = Year, y = extProb.5, color = "Forest and grassland"), linewidth=1)+
  scale_color_manual(values = c("2 Bears" = "darkgreen", 
                                "5 Bears" = "green2", 
                                "10 Bears" = "yellow3", 
                                "Settlement par. 3" = "red2",
                                "Survival dens. dep." = "orchid3", 
                                "Forest and grassland" = "deepskyblue3")) +
  ylim(0, 1)

dev.off()

# Maps together for all Sensitivity analysis
# Occupancy Prob
pdf(file = "Output_Maps/occupancy_prob.pdf") 

par(mfrow=c(3,2))
mycol_coltime <- colorRampPalette(c('ivory2', 'gold', 'olivedrab3', 'black'))
#Occupancy Prob: 2 bears
plot(new.raster, axes=F, legend=F, col=c('lightgrey'), main = "2 Bears")
plot(col_stats$map_occ_prob, add=TRUE, plg=list( title="Occupancy\nProbability", title.cex=0.9), 
     pax=list(side=1:4, retro=TRUE), axes=F, range=c(0,1), col=mycol_occprob(7), type="continuous")
text(forest_patches_pol,labels=values(forest_patches_pol)[,1])
#Occupancy Prob: 5 bears
plot(new.raster, axes=F, legend=F, col=c('lightgrey'), main = "5 Bears")
plot(col_stats_5bears$map_occ_prob, add=TRUE, plg=list( title="Occupancy\nProbability", title.cex=0.9), 
     pax=list(side=1:4, retro=TRUE), axes=F, range=c(0,1), col=mycol_occprob(7), type="continuous")
text(forest_patches_pol,labels=values(forest_patches_pol)[,1])
#Occupancy Prob: 10 bears
plot(new.raster, axes=F, legend=F, col=c('lightgrey'), main = "10 Bears")
plot(col_stats_10bears$map_occ_prob, add=TRUE, plg=list( title="Occupancy\nProbability", title.cex=0.9), 
     pax=list(side=1:4, retro=TRUE), axes=F, range=c(0,1), col=mycol_occprob(7), type="continuous")
text(forest_patches_pol,labels=values(forest_patches_pol)[,1])
#Occupancy Prob: Settlement Parameter 3
plot(new.raster, axes=F, legend=F, col=c('lightgrey'), main = "Settlement par. 3")
plot(col_stats_set_par_3$map_occ_prob, add=TRUE, plg=list( title="Occupancy\nProbability", title.cex=0.9),
     pax=list(side=1:4, retro=TRUE), axes=F, range=c(0,1), col=mycol_occprob(7), type="continuous")
text(forest_patches_pol,labels=values(forest_patches_pol)[,1])
#Occupancy Prob: Stg T = Survival Dependency Dependence
plot(new.raster, axes=F, legend=F, col=c('lightgrey'), main = "Survival density dependent")
plot(col_stats_stg_T$map_occ_prob, add=TRUE, plg=list( title="Occupancy\nProbability", title.cex=0.9), 
     pax=list(side=1:4, retro=TRUE), axes=F, range=c(0,1), col=mycol_occprob(7), type="continuous")
text(forest_patches_pol,labels=values(forest_patches_pol)[,1])
#Occupancy Prob: forest and grassland
plot(new.raster, axes=F, legend=F, col=c('lightgrey'), main = "Forest and grassland")
plot(col_stats_with_grassland$map_occ_prob, add=TRUE, plg=list( title="Occupancy\nProbability", title.cex=0.9),
     pax=list(side=1:4, retro=TRUE), axes=F, range=c(0,1), col=mycol_occprob(7), type="continuous")
text(forest_and_grassland_patches_pol,labels=values(forest_and_grassland_patches_pol)[,1])

dev.off()


# Colonisation Time
pdf(file = "Output_Maps/col_time.pdf") 

par(mfrow=c(3,2))
mycol_coltime <- colorRampPalette(c('ivory2', 'gold', 'olivedrab3', 'black'))
#Colonization Time: 2 bears
plot(new.raster, axes=F, legend=F, col=c('lightgrey'), main = "2 Bears")
plot(col_stats$map_col_time, add=TRUE, plg=list( title="Colonization\nTime in years", title.cex=0.9), 
     axes=F, breaks=c(-9,seq(-9,50,length=11)), col=mycol_coltime(7), type="continuous")
text(forest_patches_pol,labels=values(forest_patches_pol)[,1])
#Colonization Time: 5 bears
plot(new.raster, axes=F, legend=F, col=c('lightgrey'), main = "5 Bears")
plot(col_stats_5bears$map_col_time, add=TRUE, plg=list( title="Colonization\nTime in years", title.cex=0.9), 
     axes=F, breaks=c(-9,seq(-9,50,length=11)), col=mycol_coltime(7), 
     type="continuous")
text(forest_patches_pol,labels=values(forest_patches_pol)[,1])
#Colonisation Time: 10 bears
plot(new.raster, axes=F, legend=F, col=c('lightgrey'), main = "10 Bears")
plot(col_stats_10bears$map_col_time, add=TRUE, plg=list( title="Colonization\nTime in years", title.cex=0.9), 
     axes=F, breaks=c(-9,seq(-9,50,length=11)), col=mycol_coltime(7), 
     type="continuous")
text(forest_patches_pol,labels=values(forest_patches_pol)[,1])
#Colonization Time: Settlement Parameter 3
plot(new.raster, axes=F, legend=F, col=c('lightgrey'), main = "Settlement Parameter 3")
plot(col_stats_set_par_3$map_col_time, add=TRUE, plg=list( title="Colonization\nTime in years", title.cex=0.9),
     pax=list(side=1:4, retro=TRUE), axes=F, range=c(0,1), col=mycol_coltime(7), 
     type="continuous")
text(forest_patches_pol,labels=values(forest_patches_pol)[,1])
#Colonisation Time: Stg T = Survival Dependency Dependence
plot(new.raster, axes=F, legend=F, col=c('lightgrey'), main = "Survival Density Dependence")
plot(col_stats_stg_T$map_col_time, add=TRUE, plg=list( title="Colonization\nTime in years", title.cex=0.9), 
     axes=F, breaks=c(-9,seq(-9,50,length=11)), col=mycol_coltime(7), 
     type="continuous")
text(forest_patches_pol,labels=values(forest_patches_pol)[,1])
#Colonisation Time: forest and grassland
plot(new.raster, axes=F, legend=F, col=c('lightgrey'), main = "Forest and grassland")
plot(col_stats_with_grassland$map_col_time, add=TRUE, plg=list( title="Colonization\nTime in years", title.cex=0.9),
     pax=list(side=1:4, retro=TRUE), axes=F, range=c(0,1), col=mycol_coltime(7), type="continuous")
text(forest_and_grassland_patches_pol,labels=values(forest_and_grassland_patches_pol)[,1])

dev.off()




# ---------
# plot the resulting abundances and occupancy

# Abundance
pdf(file = "Output_Maps/Abundance.pdf") 

par(mfrow=c(3,2)) 
RunRS(s, dirpath)
plotAbundance(s, dirpath,sd=T, rep=F, main="2 Bears")
RunRS(s_5bears, dirpath)
plotAbundance(s_5bears, dirpath,sd=T, rep=F, main="5 Bears")
RunRS(s_10bears, dirpath)
plotAbundance(s_10bears, dirpath,sd=T, rep=F, main="10 Bears")
RunRS(s_settle_3,dirpath)
plotAbundance(s_settle_3, dirpath,sd=T, rep=F, main="Settlement par. 3")
RunRS(s_stg_T,dirpath)
plotAbundance(s_stg_T, dirpath,sd=T, rep=F, main="Survival density dependent")
RunRS(s_with_grassland)
plotAbundance(s_with_grassland, dirpath,sd=T, rep=F, main="Forest and grassland")

dev.off()

# Occupancy
pdf(file = "Output_Maps/Occupancy.pdf") 

par(mfrow=c(3,2))
RunRS(s, dirpath)
plotOccupancy(s, dirpath,sd=T, rep=F, main="2 Bears")
RunRS(s_5bears, dirpath)
plotOccupancy(s_5bears, dirpath,sd=T, rep=F, main="5 Bears")
RunRS(s_10bears, dirpath)
plotOccupancy(s_10bears, dirpath,sd=T, rep=F, main="10 Bears")
RunRS(s_settle_3,dirpath)
plotOccupancy(s_settle_3, dirpath,sd=T, rep=F, main="Settlement par. 3")
RunRS(s_stg_T,dirpath)
plotOccupancy(s_stg_T, dirpath,sd=T, rep=F, main="Survival density dependent")
RunRS(s_with_grassland)
plotOccupancy(s_with_grassland, dirpath,sd=T, rep=F, main="Forest and grassland")

dev.off()

# mean time to extinction
print("mean time to extinction:")
print("2 Bears:")
Calc_ExtTime(pop)
print("5 Bears:")
Calc_ExtTime(pop_5bears)
print("10 Bears:")
Calc_ExtTime(pop_10bears)
print("Settlement par. 3:")
Calc_ExtTime(pop_set_par_3)
print("Survival density dependence:")
Calc_ExtTime(pop_stg_T)
print("Forest and grassland:")
Calc_ExtTime(pop_with_grassland)