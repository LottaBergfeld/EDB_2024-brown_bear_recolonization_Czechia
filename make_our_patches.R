setwd("C:/Users/lotta/OneDrive/Dokumente/CLEWS_Uni_Potsdam/4_Semester/Ecosystem_Dynamics_and_Biodiversity/Seminar_paper/R-Code")

source("make.patches.R")

library(RColorBrewer)
library(RangeShiftR)
library(terra)
library(viridis)
library(RColorBrewer)
library(paletteer)
library(grid)
library(gridExtra)
library(tidyverse)
library(ggplot2) 

# set directory for input data
data_dir <- "Inputs"

# read in spatial data and make patch map - we use the raster package as terra ascii format is problematic for RangeShiftR
corine <- terra::rast(file.path(data_dir, "Corine_2018_4km_Czechia.asc"))

# land cover classes: 
# 1 "Artificial surfaces (e.g. urban areas)"
# 2 "Arable land"
# 3 "Pastures and heterogeneous agricultural areas"
# 4 "Forests"
# 5 "Grasslands, moors and heathland"
# 6 "Open spaces with little or no vegetation (e.g. beaches, rocks, burnt areas)"
# 7 "Wetlands and water bodies"

sort(unique(values(corine)))

# create patch file with "Forests"----------------------------------------------

# Make new raster for woodland
woodland <- corine

# woodland has class 4. Set all other classes to zero and then afterwards set woodland to 1
values(woodland)[values(woodland)==1] <- 0 # set all non-woodland cells to 0
values(woodland)[values(woodland)==2] <- 0 # set all non-woodland cells to 0
values(woodland)[values(woodland)==3] <- 0 # set all non-woodland cells to 0
values(woodland)[values(woodland)==5] <- 0 # set all non-woodland cells to 0
values(woodland)[values(woodland)==6] <- 0 # set all non-woodland cells to 0
values(woodland)[values(woodland)==7] <- 0 # set all non-woodland cells to 0

values(woodland)[values(woodland)==4] <- 1 # set all woodland cells to 1

# plot woodland map
plot(woodland)

# write raster to file
terra::writeRaster(woodland, filename=file.path(data_dir, "woodland_2018_4km_Czechia.asc"),
                  datatype= 'INT2S',  NAflag = -999, overwrite=T)

##########################################

# discern woodland patches

# raster file to apply function on
file.in <- file.path(data_dir, "woodland_2018_4km_Czechia.asc")

# where to store the new file
file.out <- file.path(data_dir, "woodland_patchIDs_2018_4km_Czechia.asc")

# maximal patch size accepted
maxPatchSize <- 36 # in cell numbers

# minimal patch size accepted
minPatchSize <- 13 # in cell numbers

# size (in (largest) diagonal) of the hexagonal blocks in metres
block_size <- 29099 # in m be aware of your cell resolution and calculate the length of the largest diagonal in hexagon accordingly
#block_size <- 90000 # in m be aware of your cell resolution and calculate the length of the largest diagonal in hexagon accordingly

# integer indicating which cells are considered adjacent: Should be 8 (Queen's case) or 4 (Rook's case)
directions <- 4

# apply function
new.raster <- make.patches(file.in, maxPatchSize, minPatchSize, block_size, directions)
values(new.raster)<-values(new.raster)-1
plot(new.raster)

plot(new.raster,axes=F, legend=F, col = c('grey',rep(brewer.pal(n = 12, name = "Paired"),4)))
patches_pol <- as.polygons(new.raster, dissolve=T) # Makes spatial polygons
text(patches_pol,labels=values(patches_pol)[,1])

# write resulting patch map to file
terra::writeRaster(new.raster, filename = file.out,
                  datatype= 'INT2S',  NAflag = -999, overwrite=T)


# create patch file with "Forests" and "Grasslands, moors and heathland"--------

# Make new raster for wood and grassland
wood_and_grassland <- corine

# woodland has class 4. Set all other classes to zero and then afterwards set woodland to 1
values(wood_and_grassland)[values(wood_and_grassland)==1] <- 0  
values(wood_and_grassland)[values(wood_and_grassland)==2] <- 0 
values(wood_and_grassland)[values(wood_and_grassland)==3] <- 0 
values(wood_and_grassland)[values(wood_and_grassland)==6] <- 0 
values(wood_and_grassland)[values(wood_and_grassland)==7] <- 0 

values(wood_and_grassland)[values(wood_and_grassland)==4] <- 1 
values(wood_and_grassland)[values(wood_and_grassland)==5] <- 1

# plot wood and grassland map
plot(wood_and_grassland)

# write raster to file
terra::writeRaster(wood_and_grassland, filename=file.path(data_dir, "wood_and_grassland_2018_4km_Czechia.asc"),
                   datatype= 'INT2S',  NAflag = -999, overwrite=T)

##########################################

# discern woodland patches

# raster file to apply function on
file.in <- file.path(data_dir, "wood_and_grassland_2018_4km_Czechia.asc")

# where to store the new file
file.out <- file.path(data_dir, "wood_and_grassland_patchIDs_2018_4km_Czechia.asc")

# integer indicating which cells are considered adjacent: Should be 8 (Queen's case) or 4 (Rook's case)
directions <- 4

# apply function
new.raster <- make.patches(file.in, maxPatchSize, minPatchSize, block_size, directions)
values(new.raster)<-values(new.raster)-1
plot(new.raster)

plot(new.raster,axes=F, legend=F, col = c('grey',rep(brewer.pal(n = 12, name = "Paired"),4)))
patches_pol <- as.polygons(new.raster, dissolve=T) # Makes spatial polygons
text(patches_pol,labels=values(patches_pol)[,1])

# write resulting patch map to file
terra::writeRaster(new.raster, filename = file.out,
                   datatype= 'INT2S',  NAflag = -999, overwrite=T)
