## function -----------------------------------------------------------------
#'  Splitting patches into smaller ones
#'  
#'  Necessary R packages will be installed automatically.
#'  
#'  One limitation of RangeShifter is that during the reproduction phase,
#'  juveniles are assigned an initial location within the same patch 
#'  but not necessarily the same cell as the mother.
#'  This issue can potentially result in a juvenile commencing dispersal from
#'  the opposite end of the natal patch to the mother’s nominal location,
#'  giving a potentially false impression of dispersal, but this is only of 
#'  substantial concern when habitat patches are large in extent or 
#'  very elongated.
#'  
#'  This function takes the patch landscape file (file.in parameter) and 
#'  splits patches larger than the maxPatchSize parameter [nb of cells] into hexagonal
#'  patches with a given diameter of block_size parameter [m] given in meters.
#'  
#'  The function assigns new patch IDs to the divided patches. If a hexagonal
#'  block includes less cells than the minPatchSize parameter [nb of cells], it will assign
#'  these cells to the nearest neighbouring patch.

make.patches <- function(file.in, maxPatchSize, minPatchSize, block_size, directions){
  
  require(terra)
  require(blockCV)
  require(sf)
  require(data.table)
  
  ## Prerequesites ----------------------------------------------------------
  
  # load raster and set CRS
  orig <- terra::rast(file.in)
  crs(orig) <- "EPSG:7845"
  
  rasters <- terra::patches(orig, directions=4, allowGaps=FALSE, zeroAsNA=T)
  
  plot(rasters, main="Initial patches")
  print(paste("Initial number of patches detected:",length(unique(values(rasters)))))
  
  # set background to 0
  values(rasters)[is.na(values(rasters))] <- 0
  
  # help file to add new patch ids
  rasters.xyz <- as.data.table(rasters, xy = TRUE)
  colnames(rasters.xyz) <- c("x","y","z")
  
  # extract only large patches
  largePatchIDs <- as.numeric(names(table(values(rasters)[values(rasters) > 0]))[which(table(values(rasters)[values(rasters) > 0]) > maxPatchSize)])
  
  # maximum of patch ID currently assigned
  max_patch_id <- max(values(rasters))
  
  ## Split large patches ---------------------------------------------------
  
  print(paste("Number of large patches that need splitting:",length(largePatchIDs)))
  print("Note that splitting can take a while.")
  
  # circle through large patches, split them and add new ID to newly split patches
  for (i in seq_len(length(largePatchIDs))) {
    print(paste("Starting to split patch",i))
    max.test <- max_patch_id # for debugging
    # subset rasters to only current PatchID
    rasters.tmp <- rasters
    values(rasters.tmp)[values(rasters.tmp) != largePatchIDs[i]] <- NA
    
    # transform habitat map raster to xzy data
    points <- as.data.frame(rasters.tmp, xy = TRUE)
    colnames(points) <- c("x", "y", "occ")
    
    # transform to correct data format
    pa_data <- sf::st_as_sf(points, coords = c("x", "y"), crs = 7845)
    
    # create polygon blocks
    sb1 <- cv_spatial(x = pa_data,
                      column = "occ", # the response column (binary or multi-class)
                      k = 2, # number of folds
                      size = block_size, # size of the blocks in meters
                      selection = "random", # random blocks-to-fold
                      iteration = 50, # find evenly dispersed folds
                      biomod2 = F, # also create folds for biomod2
                      report = F,
                      plot = F) # avoid plotting
    
    # plot raster and polygon blocks
    cv_plot(cv = sb1,
            r = rasters.tmp,
            raster_colors = terrain.colors(10, alpha = 0.5),
            label_size = 4)
    
    # extract block ids
    blocks <- sb1$blocks[,-2]
    
    # update patch ID for each generated block
    for (block in blocks$block_id) {
      # subset
      block.tmp <- blocks[blocks$block_id == block,]
      v <- terra::vect(block.tmp)
      
      # crop and mask
      cm <- terra::crop(terra::rast(rasters.xyz, type = "xyz"), v, mask = T) 
      
      
      # plot(cm)
      values(cm)[values(cm) != largePatchIDs[i]] <- 0
      # only if patch is large enough
      # AND!!! cells are clumped, otherwise reject
      # freq of values
      (freq.cm <- freq(cm))
      if (length(freq.cm[freq.cm$value == largePatchIDs[i],3]) > 0) { # are there cells with the specific ID?
        if (nrow(freq(patches(cm, directions = directions, zeroAsNA = TRUE, allowGaps = TRUE))) == 1) { # is there only one connected patch?
          if (freq.cm[freq.cm$value == largePatchIDs[i],3] >= minPatchSize) { # is this patch large enough?
            
            # increase the max of patch id
            max_patch_id <- max_patch_id + 1
            
            # set new patch id value
            values(cm)[values(cm) != largePatchIDs[i]] <- 0 # just to make sure there are no other patch ID in this subset
            values(cm)[values(cm) == largePatchIDs[i]] <- max_patch_id
            
            # assign new patch ID in rasters object
            
            # transform cm to xyz
            cm.xyz <- as.data.table(cm, xy = TRUE)
            colnames(cm.xyz) <- c("x","y","z_new")
            cm.xyz <- cm.xyz[z_new == max_patch_id]
            # merge with rasters
            rasters.xyz <- merge(rasters.xyz,cm.xyz, all.x = T)
            # update value
            rasters.xyz[!is.na(z_new) & z == largePatchIDs[i],z := z_new] # test values
            # delete temporary column
            rasters.xyz[,z_new := NULL]
            # plot(terra::rast(rasters.xyz, type="xyz"))
            # freq(terra::rast(rasters.xyz, type="xyz"))
          } else{
            # print(paste0("patch too small in block ", block))
          }
        }else{
          # choose the biggest patch -> is it large enough? -> then select these cells for new patch
          if (max(freq(patches(cm, directions = directions, zeroAsNA = TRUE, allowGaps = TRUE))$count) >= minPatchSize) {
            # tmp includes all patches
            tmp <- patches(cm, directions = directions, zeroAsNA = TRUE, allowGaps = TRUE)
            tmp_sizes <- freq(tmp)
            tmp_sizes <- tmp_sizes[tmp_sizes$count >= minPatchSize,]
            for (j in tmp_sizes$value) {
              max_patch_id <- max_patch_id + 1
              values(tmp)[values(tmp) != j] <- 0 # just to make sure there are no other patch ID in this subset
              values(tmp)[values(tmp) == j] <- max_patch_id
              
              cm.xyz <- as.data.table(tmp, xy = TRUE)
              colnames(cm.xyz) <- c("x","y","z_new")
              cm.xyz <- cm.xyz[z_new == max_patch_id]
              # merge with rasters
              rasters.xyz <- merge(rasters.xyz,cm.xyz, all.x = T)
              # update value
              rasters.xyz[!is.na(z_new) & z == largePatchIDs[i],z := z_new] # test values
              # delete temporary column
              rasters.xyz[,z_new := NULL]
            }
            
          } else {
            # print(paste0("more than one patch in block ", block," and all too small"))
            }
        }
      }else{
        # print(paste0("no cells of patch in block ", block))
      }
    } # loop over blocks
    
    # current state after splitting one patch
    rasters.splitted <- terra::rast(rasters.xyz, type = "xyz")
    # plot(rasters.splitted)
    
    # check if there are still cells with the current patch ID of the large patch
    rasters.tmp <- rasters.splitted
    
    if (nrow(freq(rasters.tmp)[freq(rasters.tmp)$value == largePatchIDs[i],]) > 0) {
      
      # set values to NA
      values(rasters.tmp)[values(rasters.tmp) == largePatchIDs[i]] <- NA
      
      # function
      expand <- function(x){if (any(!is.na(x)) && max(x, na.rm = T) > 0) {return(max(x, na.rm = T))} else{return(NA)}}
      if (directions == 4) {w <- matrix(c(NA,1,NA,1,1,1,NA,1,NA),3,3)# directions = 4
      } else {w=matrix(c(1,1,1,1,1,1,1,1,1),3,3) # directions = 8
              }
      
      
      
      # set value to maximal value in 3x3 cells as long as there are still NAs
      iteration <- 1
      
      while (nrow(table(is.na(rasters.tmp[]))) > 1 & iteration < 500) {
        
        rasters.tmp <- terra::focal(rasters.tmp, w = w, fun = expand, #fill.na, 
                                    na.policy = "only" )
        iteration <- iteration + 1      
        (table(is.na(rasters.tmp[])))
        
      }
      
      # if (nrow(table(is.na(rasters.tmp[])))>1)stop("NAs still left in landscape")
      
      col.rasters <- colnames(rasters.xyz)
      
      rasters.xyz <- as.data.table(rasters.tmp, xy = TRUE)
      
      colnames(rasters.xyz) <- col.rasters
      
    } # if blocks were too small for a new patch
    
  } # loop over large patches
  
  # update patch IDs
  colnames(rasters.xyz) <- c("x", "y", "z")
  
  id_count <- 1
  
  for (nb in unique(rasters.xyz$z)) {
    
    if (nb > 0) {
      
      rasters.xyz[z == nb,z_new := id_count]
      
      id_count <- id_count + 1
      
    } else {
      rasters.xyz[z == nb,z_new := z]
      }
  }
  
  rasters.xyz[,z := z_new]
  
  rasters.xyz[,z_new := NULL]
  
  rasters <- terra::rast(rasters.xyz, type = "xyz")
  
  # mask again to original mask
  rasters.old <- terra::rast(file.in)
  rasters <- mask(rasters, rasters.old)
  
  patch.sizes <- values(rasters) %>% table
  values(rasters)[values(rasters) %in% as.numeric(names(patch.sizes)[which(patch.sizes < minPatchSize)])] <- 0
  
  id_count <- 1
  
  for (nb in sort(unique(values(rasters)))) {
      
    values(rasters)[values(rasters)==nb] <- id_count
      id_count <- id_count + 1
  }
  
  return(rasters)
  
}
