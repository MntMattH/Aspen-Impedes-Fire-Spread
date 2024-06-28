#Aspen impedes wildfire spread analysis methods and results: 

#Load packages: 
library(dplyr)
library(terra)
library(sf)
library(sp)
library(terra)
library(viridis)
library(ggplot2)
library(tidyterra)
library(landscapemetrics)
library(stringr)
library(glmmTMB)
library(effects)
library(MuMIn)
library(ggpubr)


#Create linear spread rasters based on Day of Burning (DOB) interpolations developed by Sean Parks. 
for (year in 2002:2020) {
  
  ## Use the DOB interpolation output directory to get a list of all successfully interpolated fires for current year.
  fire.list <- list.dirs(path = paste0('C:/Users/matth/Desktop/DOB/', year, '/'),
                         full.names = FALSE, recursive = TRUE)
  fire.list <- fire.list[-1] 
  
  ### Loop over fire.list to rasterize forward distances.
  for (fire in fire.list) { 
    
    # Enter fire's folder
    setwd(paste0('C:/Users/matth/Desktop/DOB/', year, '/', fire))
    
    # Load DOB interpolation and hotspots, confirm they are on same projection 
    if (file.exists(paste0(fire,"_dob.QC.tif"))) {
      dob.rast <- rast(paste0(fire,"_dob.QC.tif")) # if the fire was reinterpolated after QCAC, use that interpolation.
    } else {
      dob.rast <- rast(paste0(fire,"_dob.tif")) # Otherwise grab default interpolation if the fire was not pinged by QCAC.
    }
    
    
    hotspots <- vect(paste0(fire,"_hotspots.shp"))
    hotspots <- project(hotspots, dob.rast)
    
    ## Create DOBLOB raster of spatially discretized Day of Burn patches, or contiguous pixels burning on the same day.
    # To do this we'll use the get_patches() function from the landscapemetrics package.
    doblob.rast <- get_patches(dob.rast, class = "all", directions = 8, return_raster = TRUE) 
    # get_patches() actually returns a list of raster layers (1 layer for each DOB containing its DOBLOB patches). 
    # Let's rename these layers for what they are; DOBs. (Default is something nondescript like "layer_1")
    names(doblob.rast) <- "DOB"
    # Next let's squish those DOB layers containing their DOBLOB patches into a single layer with merge(). We'll also
    # use Reduce() to hit all layers at once. Finally, use rast() to make this object a spatraster again.
    doblob.rast <- rast(Reduce(merge, doblob.rast$DOB)) 
    # Let's rename this newly merged single layer for what it contains; DOBLOBS
    names(doblob.rast) <- "DOBLOB"
    
    
    ## Convert rasters to polygons for spatial analyses via as.polygons()
    doblob.poly <- as.polygons(doblob.rast)
    dob.poly <- as.polygons(dob.rast)
    # Convert spatvectors to sf objects. 
    doblob.poly <- st_as_sf(doblob.poly)
    dob.poly <- st_as_sf(dob.poly)
    
    
    #### Loop through each DOB within current fire. Iterator variable 'i' used to 
    #### process the fire's DOBs in chronological order.
    for (i in 1:dim(dob.poly)[1]) { 
      
      # subset DOBLOBs of current DOB. 
      doblobs <- doblob.poly[dob.poly[i,], ,op=st_within]
      
      ############## if statement to process Day/DOB 1 as special case (i.e., dob.poly i == 1)
      if (i == 1) {
        
        ### Prep objects to store rasterized distance. 
        # output.rast.forward = distance in the direction of fire progression,
        # ie. towards the sequential day
        output.rast.forward <- doblob.rast
        values(output.rast.forward) <- NA
        
        ### Process DOBLOBs belonging to DOB==1 by looping iterator variable
        # z over the doblobs row index
        for (z in 1:dim(doblobs)[1]) {
          # Day 1 case handling: For now I'll use the first hotspots as the fire
          # origin. 
          
          # Prepare a raster mask of DOBLOB == z to clip the distance rasters.
          doblob.mask <- mask(dob.rast, doblobs[(z),], touches = FALSE)
          
          # Prepare raster of DOBLOB == z perimeter (DOBLOB == z). First must
          # cast the polygon perimeter to a multilinestring (b/c we don't care
          # about polygon interior). Second, use vect() to convert this sf
          # object back to a spatvector for terra::rasterize(). Use doblob.mask
          # as the second argument to terra::rasterize() which defines spatial 
          # extent & resolution of the rasterization.
          doblob.perim <- st_cast(doblobs[(z),],"MULTILINESTRING")
          doblob.perim <- rasterize(vect(doblob.perim),
                                    doblob.mask)
          
          # Prepare earliest hotspots (minimum datetime) and use these as the 
          # defacto starting place for the fire. Even if the fire didn't actually
          # start here, they are where our DOB interpolation 'starts'.
          earliest.hotspots <- subset(hotspots, 
                                      hotspots$loc_JDT == min(hotspots$loc_JDT))
          
          # Check whether hotspots occur inside or outside of the DOBLOB; this 
          # determines how forward distances are estimated. The line 
          # below will count how many hotspots occur inside the DOBLOB.
          if (nrow(earliest.hotspots[vect(doblobs[z,]),]) == 0){           
            ## If none of the earliest hotspots occur inside the DOBLOB, find
            # the spot on the DOBLOB perimeter closest to a hotspot detection.
            
            # Cast the disjoint hotspots spatvector object to a sf object.
            disjoint.hotspots <- st_as_sf(earliest.hotspots)
            
            # sf::st_nearest_points() creates linestrings connecting points
            # (here, any disjoint hotspots) to a polygon perimeter (here, DOBLOB
            # == z). We'll use the location of this "connection" as the starting 
            # point for forward distance rasterization.
            disjoint.hotspots <- st_nearest_points(disjoint.hotspots, 
                                                   doblobs[z,])
            
            # Find the distance from the nearest hotspot to the DOBLOB perimeter
            # and prepare to select all hotspots within 100m of that distance
            # (from any side of the DOBLOB; this is not a buffer of the nearest
            # hotspot!)
            nearest.distance <- as.vector(min(st_length(disjoint.hotspots))) + 100
            
            # now grab any hotspots within the distance parameterized above.
            # which() will create a vector of indicies for elements closer than
            # the nearest.distance.
            nearest.distance <- which(as.vector(st_length(disjoint.hotspots)) <= nearest.distance)
            disjoint.hotspots <- disjoint.hotspots[nearest.distance]
            
            # Find points where the linestrings created above connect with the
            # DOBLOB. 
            disjoint.hotspots <- st_intersection(doblobs[z,], 
                                                 disjoint.hotspots)
            
            ### Add DOBLOB's "Forward Distance" to output.rast.forward
            # Calculate rasterized distance from disjoint hotspot detections to 
            # DOBLOB border. Call this object tmp.dist since we don't need to 
            # keep it and can reuse the name later to keep code & working 
            # environment a bit cleaner.
            tmp.dist <- distance(doblob.mask, disjoint.hotspots, unit = "m")
            tmp.dist <- mask(tmp.dist, doblob.mask) # Mask values beyond current DOB.
            # add to output.rast.forward
            output.rast.forward <- cover(tmp.dist, output.rast.forward)
            
          } else {
            # Else, the earliest hotspots occur inside the DOBLOB and we can use
            # them as-is. Calculate rasterized distance from hotspot detection to 
            # DOBLOB border. 
            earliest.hotspots <- earliest.hotspots[vect(doblobs[z,]),]
            tmp.dist <- distance(doblob.mask, earliest.hotspots, unit = "m")
            tmp.dist <- mask(tmp.dist, doblob.mask) # Mask values beyond current DOB.
            # add to output.rast.forward
            output.rast.forward <- cover(tmp.dist, output.rast.forward)
          }
        } # END DOB i = 1 PROCESSING
        
        ############### if...else, process Day/DOB i > 1 (i.e., dob.poly i > 1). There
        ############### will be three DOBLOB cases to check for and proc as appropriate: 
        ############### 1) temporal gaps, 2) non-gap contiguous and 3) non-gap disjoint.
      } else { 
        
        #*# Recall DOBLOBs for DOB = i were prepared above @ line 143, "doblobs"
        
        ############### 1) Check for temporal gap in fire progression. If we have a temporal
        # gap, we will use hotspots instead of the shared border to rasterize
        # forward distance. To check for a temporal gap, check whether the 
        # difference between DOB dates is greater than 1.
        if (dob.poly[i,]$DOB - dob.poly[i-1,]$DOB > 1) {
          
          # if this if statement is triggered, we want to rasterize spread for 
          # these DOBLOBs using their hotspot(s). Use iterator variable z to 
          # loop over all DOBLOBs for DOB i.
          for (z in 1:dim(doblobs)[1]) { 
            gc()
            # Prepare mask of current DOBLOB
            doblob.mask <- mask(dob.rast, 
                                doblobs[(z),], 
                                touches = FALSE)
            
            # prepare rasterized boundary of current DOBLOB
            doblob.perim <- st_cast(doblobs[(z),],
                                    "MULTILINESTRING")
            doblob.perim <- rasterize(doblob.perim, doblob.mask)
            
            # Because the starting point of the fire is unknown,
            # use the assumption of the hotspot detection location(s).
            dob.hotspots <- subset(hotspots, 
                                   hotspots$date == dob.poly[i,]$DOB | 
                                     (hotspots$date-1) == dob.poly[i,]$DOB) 
            
            # If none of the DOB's hotspots occur inside the DOBLOB, find
            # the spot on the DOBLOB perimeter closest to a hotspot detection.
            # The line below will count how many hotspots occur inside the DOBLOB.
            if (nrow(dob.hotspots[vect(doblobs[z,]),]) == 0){
              
              # st_nearest_points() fxn creates linestrings connecting points
              # (here, disjoint hotspots) to the current DOBLOB. We'll use the 
              # location of the "connection" as the starting point for distance
              # rasterization.
              dob.hotspots <- st_as_sf(dob.hotspots)
              dob.hotspots <- st_nearest_points(dob.hotspots, 
                                                doblobs[z,])
              
              # Find the distance from the nearest hotspot to the DOBLOB perimeter
              # and prepare to select all hotspots within 100m of that distance
              # (from any side of the DOBLOB; this is not a buffer of the nearest
              # hotspot!)
              nearest.distance <- as.vector(min(st_length(dob.hotspots))) + 100
              
              # now grab any hotspots within the distance parameterized above.
              # which() will create a vector of indicies for elements closer than
              # the nearest.distance.
              nearest.distance <- which(as.vector(st_length(dob.hotspots)) <= nearest.distance)
              dob.hotspots <- dob.hotspots[nearest.distance]
              
              # find the points where the hotspot linestring connects with the
              # DOBLOB. 
              dob.hotspots <- st_intersection(doblobs[z,], 
                                              dob.hotspots)
              
            } else {
              # else, use whichever hotspots occur inside the DOBLOB,
              
              # first grab only those inside the current DOBLOB.
              dob.hotspots <- dob.hotspots[vect(doblobs[z,]),]
              
              # second, grab the earliest hotspots
              dob.hotspots <- subset(dob.hotspots, 
                                     dob.hotspots$date == min(dob.hotspots$date))
              dob.hotspots <- subset(dob.hotspots, 
                                     dob.hotspots$time == min(dob.hotspots$time))
              dob.hotspots <- st_as_sf(dob.hotspots)
            }
            
            ### Add DOBLOB's "Forward Distance" to output.rast.forward
            # Here we'll rasterize distance from the shared boundary
            tmp.dist <- distance(doblob.mask, dob.hotspots, unit = "m")
            tmp.dist <- mask(tmp.dist, doblob.mask)
            # add this to the output raster
            output.rast.forward <- cover(output.rast.forward, tmp.dist)
            
          } 
          # END LOOP FOR TEMPORAL GAP DOBLOBS
          
        } else {
          # else, there was not a temporal gap, so rasterize distance using original
          # workflow (check whether contiguous or disjoint and proceed as appropriate)
          
          ############### 2) Next we need to know whether DOB i's DOBLOBs are contiguous or
          # disjoint. We'll start with contiguous. To identify them, subset DOBLOBs 
          # that grew from the previous day's DOB area (i.e., any of the previous 
          # DOB's DOBLOBs, fire could have 'started' from any of these and in any
          # case spread is modeled from the shared perimeter.
          contiguous.doblob.poly <- doblobs[dob.poly[i-1,], op=st_touches]
          
          ## Calculate Rasterized Spread for contiguous DOBLOB polygons; First 
          # check whether there are any contiguous DOBLOBs
          if (dim(contiguous.doblob.poly)[1] > 0) {
            
            # Use iterator variable z to loop over contiguous DOBLOBs for DOB i
            for (z in 1:dim(contiguous.doblob.poly)[1]) { 
              gc()
              # Prepare mask of current DOBLOB
              doblob.mask <- mask(dob.rast, 
                                  contiguous.doblob.poly[z,], 
                                  touches = FALSE)
              
              # Find the shared border between each DOBLOB and any of the prior 
              # DOB's DOBLOBs. This shared.border represents the flaming front 
              # along which (we presume) fire spread occurred or "began".
              shared.border <- st_intersection(contiguous.doblob.poly[z,], 
                                               dob.poly[(i-1),]) 
              
              # Now cast the DOBLOB polygon to a linestring and subtract the 
              # shared border to identify the nonshared border.
              doblob.perim <- st_cast(contiguous.doblob.poly[(z),],
                                      "MULTILINESTRING")
              
              # if shared.border is only 1 pixel, it requires special processing 
              # to allow us to calculate the nonshared border.
              if (st_geometry_type(shared.border) == "POINT" | 
                  st_geometry_type(shared.border) =="MULTIPOINT") {
                nonshared.border <- st_difference(doblob.perim, 
                                                  st_as_sf(buffer(vect(shared.border), 10)))
              } else {
                # else, it needs no special treatment, just cast to multilinestring
                # in case it is a geometry collection (occurs when shared border is
                # huge)
                shared.border <- st_cast(shared.border)[! st_is(st_cast(shared.border), c('POINT', 'MULTIPOINT')),]
                shared.border <- st_cast(shared.border, 'MULTILINESTRING')
                nonshared.border <- st_difference(doblob.perim, 
                                                  shared.border)
              }
              
              # Rasterize pieces of the doblob border to accelerate processing.
              # Define extent and resolution using the doblob.mask.
              shared.border <- rasterize(shared.border, doblob.mask)
              doblob.perim <- rasterize(doblob.perim, doblob.mask)
              
              ### Add DOBLOB's "Forward Distance" to output.rast.forward
              # Here we'll rasterize distance from the shared boundary.
              tmp.dist <- distance(shared.border, unit = "m")
              tmp.dist <- mask(tmp.dist, doblob.mask)
              # add this to the output raster
              output.rast.forward <- cover(output.rast.forward, tmp.dist)
            }
            # END LOOP FOR CONTIGUOUS DOBLOBS
          }
          
          ############### 3) Identify disjoint DOBLOBs by finding any DOBLOBs that aren't in
          # the contiguous.doblob.poly object created above. Spread in these disjoint polygons 
          # is modeled from the day's hotspot detections. 
          disjoint.doblob.poly <- doblobs[!doblobs[[1]] %in% contiguous.doblob.poly[[1]],]
          
          ## Calculate Rasterized Spread for disjoint DOBLOB polygons; First 
          # check whether there are any disjoint DOBLOBs
          if (dim(disjoint.doblob.poly)[1] > 0) {
            
            # Use iterator variable z to loop over disjoint DOBLOBs for DOB i
            for (z in 1:dim(disjoint.doblob.poly)[1]) { 
              gc()
              # Prepare mask of current DOBLOB
              doblob.mask <- mask(dob.rast, 
                                  disjoint.doblob.poly[(z),], 
                                  touches = FALSE)
              
              # prepare rasterized boundary of current DOBLOB
              doblob.perim <- st_cast(disjoint.doblob.poly[(z),],
                                      "MULTILINESTRING")
              doblob.perim <- rasterize(doblob.perim, doblob.mask)
              
              
              # Because the starting point of the fire is unknown,
              # use the assumption of the hotspot detection location(s).
              dob.hotspots <- subset(hotspots, 
                                     hotspots$date == dob.poly[i,]$DOB | # Selects hotspots from this DOB
                                       (hotspots$date-1) == dob.poly[i,]$DOB) #Selects hotspots from prior DOB
              
              # If none of the DOB's hotspots occur inside the DOBLOB, find
              # the spot on the DOBLOB perimeter closest to a hotspot detection.
              # The line below will count how many hotspots occur inside the DOBLOB.
              if (nrow(dob.hotspots[vect(disjoint.doblob.poly[z,]),]) == 0){
                
                # st_nearest_points() fxn creates linestrings connecting points
                # (here, disjoint hotspots) to the current DOBLOB. We'll use the 
                # location of the "connection" as the starting point for distance
                # rasterization.
                dob.hotspots <- st_as_sf(dob.hotspots)
                dob.hotspots <- st_nearest_points(dob.hotspots, 
                                                  disjoint.doblob.poly[z,])
                
                # Find the distance from the nearest hotspot to the DOBLOB perimeter
                # and prepare to select all hotspots within 100m of that distance
                # (from any side of the DOBLOB; this is not a buffer of the nearest
                # hotspot!)
                nearest.distance <- as.vector(min(st_length(dob.hotspots))) + 100
                
                # now grab any hotspots within the distance parameterized above.
                # which() will create a vector of indicies for elements closer than
                # the nearest.distance.
                nearest.distance <- which(as.vector(st_length(dob.hotspots)) <= nearest.distance)
                dob.hotspots <- dob.hotspots[nearest.distance]
                
                # find the points where the hotspot linestrings connect with
                # the DOBLOB. 
                dob.hotspots <- st_intersection(disjoint.doblob.poly[z,], 
                                                dob.hotspots)
                
              } else {
                # else, use whichever hotspots occur inside the DOBLOB,
                
                # first grab only those inside the current DOBLOB.
                dob.hotspots <- dob.hotspots[vect(disjoint.doblob.poly[z,]),]
                
                # second, grab the earliest hotspots
                dob.hotspots <- subset(dob.hotspots, 
                                       dob.hotspots$date == min(dob.hotspots$date))
                dob.hotspots <- subset(dob.hotspots, 
                                       dob.hotspots$time == min(dob.hotspots$time))
                dob.hotspots <- st_as_sf(dob.hotspots)
              }
              
              ### Add DOBLOB's "Forward Distance" to output.rast.forward
              # Here we'll rasterize distance from the shared boundary
              tmp.dist <- distance(doblob.mask, dob.hotspots, unit = "m")
              tmp.dist <- mask(tmp.dist, doblob.mask)
              # add this to the output raster
              output.rast.forward <- cover(output.rast.forward, tmp.dist)
              plot(output.rast.forward)
              
              
              
            }
            # END LOOP FOR DISJOINT DOBLOBS
          }
        }
      } # END DOB i > 1 PROCESSING
      
    } # END SPATIALIZED FIRE SPREAD PROCESSING
    
    
    
    # Lastly, save the forward spread raster
    # setwd(paste0('C:/Users/matth/Desktop/Thesis Dev/Daily_Fire_Progression'))
    # writeRaster(output.rast.forward, 
    #             filename = paste0("rast_dist_", fire, ".tiff"),
    #             overwrite = TRUE)
    # #Save max spread rate rasters: 
    # setwd(paste0('C:/Users/matth/Desktop/Thesis Dev/Max_Dist_Rasters'))
    # writeRaster(output.rast.max.dist,
    #             filename = paste0("Max_dist_", fire, ".tiff"),
    #             overwrite = TRUE)
  } # END FIRE PROCESSING, MOVE TO NEXT FIRE IN CPU'S FIRE LIST
}

###Create 0.1% random points for sampling of covariates per DOB:

# Create a random seed for reproducibility
set.seed(369)

for (year in 2002:2020) {
  
  fire.list <- list.dirs(path = paste0('C:/Users/matth/Desktop/DOB/', year, '/'),
                         full.names = FALSE, recursive = TRUE)
  fire.list <- fire.list[-1] 
  
  ### Loop over fire.list to create 0.1% random point per DOB
  for (fire in fire.list) { 
    # Enter fire's folder. Modify file pathway to where you've stored the data.
    setwd(paste0('C:/Users/matth/Desktop/DOB/', year, '/', fire))
    # Load DOB interpolation 
    if (file.exists(paste0(fire,"_dob.QC.tif"))) {
      dob.rast <- rast(paste0(fire,"_dob.QC.tif")) # if the fire was reinterpolated after QCAC, use that interpolation.
    } else {
      dob.rast <- rast(paste0(fire,"_dob.tif")) # Otherwise grab default interpolation if the fire was not pinged by QCAC.
    }
    
    # Convert raster to points
    points <- terra::as.points(dob.rast$DOB)
    
    # Get unique days from the "DOB" variable
    unique_days <- unique(points$DOB)
    
    sampled_points_per_day <- lapply(unique_days, function(day){
      
      # Subset the data to all points for each day
      day_data <- points[points$DOB == day, ]
      
      # Calculate the number of points to sample (0.1% of total points for the day)
      num_sample <- round(nrow(day_data) * 0.001)
      
      # I bet num_sample could round to 0 for small DOBs. This ensures a minimum of
      # 1 point per DOB. Optional whether or not you use this (a true random sample
      # could miss days entirely) but it's here if you want it.
      if (num_sample <= 0) { num_sample <- 1 }
      
      # Sample 0.1% of the points randomly for the day
      sampled_points <- day_data[sample(nrow(day_data), num_sample), ]
      
      return(sampled_points)
    })
    
    #Join all sampled_points_per_day
    sampled_points_per_day 
    sampled_points_per_day <- do.call(rbind, sampled_points_per_day) 
    sampled_points_per_day 
    
    output_dir<-"C:/Users/matth/Desktop/Thesis Dev/DOB_Random_Points_01P"
    
    # Define the output file path and name 
    output_file <- file.path(output_dir, paste0(basename(fire), "_points.gpkg"))
    
    # now save to disk with terra's writeVector
    writeVector(sampled_points_per_day,
                output_file)
    
    
  }
}  

###Create a loop to extract daily weather, climate, and topography for all DOBS: 
#Grab list of fires 
setwd("C:/Users/matth/Desktop/Thesis Dev/Master_DF")
Master_DF<-read.csv("Master_df.csv")
#Seperate by year:
YDF<-subset(Master_DF, Master_DF$FireYear==2020)
Fire_List<-unique(YDF$Fire_ID)
Fire_List <- paste0(Fire_List, "_points.gpkg")

#       Running Year: 2020
#Load stacks of each variable per year: 
setwd("C:/Users/matth/Desktop/Data Dungeon/ERA5_Daily_Weather/ERA5_NorthAmerica_tiffs/Outputs/build_up_index")
bui.stack<-rast("build_up_index_2020.tiff")
setwd("C:/Users/matth/Desktop/Data Dungeon/ERA5_Daily_Weather/ERA5_NorthAmerica_tiffs/Outputs/daily_severity_rating")
dsr.stack<-rast("daily_severity_rating_2020.tiff")
setwd("C:/Users/matth/Desktop/Data Dungeon/ERA5_Daily_Weather/ERA5_NorthAmerica_tiffs/Outputs/drought_code")
dc.stack<-rast("drought_code_2020.tiff")
setwd("C:/Users/matth/Desktop/Data Dungeon/ERA5_Daily_Weather/ERA5_NorthAmerica_tiffs/Outputs/duff_moisture_code")
dmc.stack<-rast("duff_moisture_code_2020.tiff")
setwd("C:/Users/matth/Desktop/Data Dungeon/ERA5_Daily_Weather/ERA5_NorthAmerica_tiffs/Outputs/fine_fuel_moisture_code")
ffm.stack<-rast("fine_fuel_moisture_code_2020.tiff")
setwd("C:/Users/matth/Desktop/Data Dungeon/ERA5_Daily_Weather/ERA5_NorthAmerica_tiffs/Outputs/fire_weather_index")
fwi.stack<-rast("fire_weather_index_2020.tiff")
setwd("C:/Users/matth/Desktop/Data Dungeon/ERA5_Daily_Weather/ERA5_NorthAmerica_tiffs/Outputs/initial_spread_index")
isi.stack<-rast("initial_spread_index_2020.tiff")
setwd("C:/Users/matth/Desktop/Data Dungeon/ERA5_Daily_Weather/ERA5_NorthAmerica_tiffs/Inputs/24hr_max_temperature")
tmax.stack<-rast("24hr_max_temperature_2020.tiff")
setwd("C:/Users/matth/Desktop/Data Dungeon/ERA5_Daily_Weather/ERA5_NorthAmerica_tiffs/Inputs/24hr_accumulated_precipitation")
prec.stack<-rast("24hr_accumulated_precipitation_2020.tiff")
setwd("C:/Users/matth/Desktop/Data Dungeon/ERA5_Daily_Weather/ERA5_NorthAmerica_tiffs/Inputs/global_noon_LST_wind_speed")
ws.stack<-rast("global_noon_LST_wind_speed_2020.tiff")
setwd("C:/Users/matth/Desktop/Data Dungeon/ERA5_Daily_Weather/ERA5_NorthAmerica_tiffs/Inputs/global_noon_LST_relative_humidity")
rh.stack<-rast("global_noon_LST_relative_humidity_2020.tiff")

#Load topography data (Not bound to days extraction)
setwd("C:/Users/matth/Desktop/Data Dungeon/Topography/LF2020_Asp_220_CONUS")
Asp.stack<-rast("LC20_Asp_220.tif")
setwd("C:/Users/matth/Desktop/Data Dungeon/Topography/LF2020_Elev_220_CONUS")
Elev.stack<-rast("LC20_Elev_220.tif")
setwd("C:/Users/matth/Desktop/Data Dungeon/Topography/LF2020_SlpD_220_CONUS")
Slp.stack<-rast("LC20_SlpD_220.tif")
setwd("C:/Users/matth/Desktop/Data Dungeon/Topography/Ruggedness")
Rugg.stack<-rast("Ruggedness.tif")


#Alter Aspect: 
northness <- cos(Asp.stack * pi / 180)
eastness <- sin(Asp.stack * pi / 180)

#Make covariates df
#Covariate_df<-data.frame()


for (fire in Fire_List) {
  #Load DOB_RP (0.1% random points)
  setwd("C:/Users/matth/Desktop/Thesis Dev/DOB_Random_Points_01P")
  DOB_RP<-st_read(fire)
  
  
  lapply(unique(DOB_RP$DOB), function(x){
    DOB_RP[DOB_RP$DOB == x, "bui_raw"]  <<- terra::extract(bui.stack[[x]], DOB_RP[DOB_RP$DOB == x,])[,2]
    DOB_RP[DOB_RP$DOB == x, "dsr_raw"]  <<- terra::extract(dsr.stack[[x]], DOB_RP[DOB_RP$DOB == x,])[,2]
    DOB_RP[DOB_RP$DOB == x, "dc_raw"]  <<- terra::extract(dc.stack[[x]], DOB_RP[DOB_RP$DOB == x,])[,2]
    DOB_RP[DOB_RP$DOB == x, "dmc_raw"]  <<- terra::extract(dmc.stack[[x]], DOB_RP[DOB_RP$DOB == x,])[,2]
    DOB_RP[DOB_RP$DOB == x, "ffm_raw"]  <<- terra::extract(ffm.stack[[x]], DOB_RP[DOB_RP$DOB == x,])[,2]
    DOB_RP[DOB_RP$DOB == x, "fwi_raw"]  <<- terra::extract(fwi.stack[[x]], DOB_RP[DOB_RP$DOB == x,])[,2]
    DOB_RP[DOB_RP$DOB == x, "isi_raw"]  <<- terra::extract(isi.stack[[x]], DOB_RP[DOB_RP$DOB == x,])[,2]
    DOB_RP[DOB_RP$DOB == x, "tmax_raw"]  <<- terra::extract(tmax.stack[[x]], DOB_RP[DOB_RP$DOB == x,])[,2]
    DOB_RP[DOB_RP$DOB == x, "prec_raw"]  <<- terra::extract(prec.stack[[x]], DOB_RP[DOB_RP$DOB == x,])[,2]
    DOB_RP[DOB_RP$DOB == x, "ws_raw"]  <<- terra::extract(ws.stack[[x]], DOB_RP[DOB_RP$DOB == x,])[,2]
    DOB_RP[DOB_RP$DOB == x, "rh_raw"]  <<- terra::extract(rh.stack[[x]], DOB_RP[DOB_RP$DOB == x,])[,2]
    #Now extract topography variables: 
    DOB_RP[DOB_RP$DOB == x, "Northness"]  <<- terra::extract(northness, DOB_RP[DOB_RP$DOB == x,])[,2]
    DOB_RP[DOB_RP$DOB == x, "Eastness"]  <<- terra::extract(eastness, DOB_RP[DOB_RP$DOB == x,])[,2]
    DOB_RP[DOB_RP$DOB == x, "Elev_raw"]  <<- terra::extract(Elev.stack, DOB_RP[DOB_RP$DOB == x,])[,2]
    DOB_RP[DOB_RP$DOB == x, "Slp_raw"]  <<- terra::extract(Slp.stack, DOB_RP[DOB_RP$DOB == x,])[,2]
    DOB_RP[DOB_RP$DOB == x, "Rugg_raw"]  <<- terra::extract(Rugg.stack, DOB_RP[DOB_RP$DOB == x,])[,2]
  })
  
  
  #Remove geom:
  DOB_RP_df<-as.data.frame(DOB_RP)
  DOB_RP_df<-DOB_RP_df[, -2]
  #Create means df: 
  DOB_RP_Mean<- DOB_RP_df %>%
    group_by(DOB) %>%
    summarise_all(mean)
  
  #Make each 3 decimal places:
  DOB_RP_Mean <- round(DOB_RP_Mean, digits = 3)
  
  #Add DOB_ID:
  DOB_ID<-as.character(str_sub_all(fire,0 ,-12))
  DOB_RP_Mean$DOB_ID<-paste0(DOB_ID,DOB_RP_Mean$DOB)
  
  
  
  Covariate_df<-rbind(Covariate_df, DOB_RP_Mean)
  
}  

#2020 Done

#Save Final Covariate_df
write.csv(Covariate_FDF, "C:\\Users\\matth\\Desktop\\Data Dungeon\\Daily_Covariate_DF_01P.csv")

###Create a loop to extract monthly climate for all DOBS: 

#Grab list of fires 
setwd("C:/Users/matth/Desktop/Thesis Dev/Master_DF")
Master_DF<-read.csv("Master_df.csv")
#Seperate by year:
YDF<-subset(Master_DF, Master_DF$FireYear==2020)
Fire_List<-unique(YDF$Fire_ID)
Fire_List <- paste0(Fire_List, "_points.gpkg")

#       Running Year: 2020
#Load stacks of each variable per year: 

#Stack monthly covariates into full years: 
#CWD:
CWD_file_path_base <- "C:/Users/matth/Desktop/Data Dungeon/TerraClimate_Raw/CWDef/def.2020_"
# List of months (assuming raster files are named accordingly)
months <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12")
# Initialize an empty list to store the monthly rasters
CWD_raster_list <- list()
# Loop through each month to load rasters and add them to the list
for (month in months) {
  CWD_file_path <- paste0(CWD_file_path_base, month, ".tiff")  
  raster <- rast(CWD_file_path)
  CWD_raster_list[[month]] <- raster
}
# Create a raster stack from the list of monthly rasters
CWD_raster_stack <- rast(CWD_raster_list)
#TMAX:
TMAX_file_path_base <- "C:/Users/matth/Desktop/Data Dungeon/TerraClimate_Raw/TMax/tmax.2020_"
# Initialize an empty list to store the monthly rasters
TMAX_raster_list <- list()
# Loop through each month to load rasters and add them to the list
for (month in months) {
  TMAX_file_path <- paste0(TMAX_file_path_base, month, ".tiff")  
  raster <- rast(TMAX_file_path)
  TMAX_raster_list[[month]] <- raster
}
# Create a raster stack from the list of monthly rasters
TMAX_raster_stack <- rast(TMAX_raster_list)
#VPD:
VPD_file_path_base <- "C:/Users/matth/Desktop/Data Dungeon/TerraClimate_Raw/VPD/vpd.2020_"
# Initialize an empty list to store the monthly rasters
VPD_raster_list <- list()
# Loop through each month to load rasters and add them to the list
for (month in months) {
  VPD_file_path <- paste0(VPD_file_path_base, month, ".tiff")  
  raster <- rast(VPD_file_path)
  VPD_raster_list[[month]] <- raster
}
# Create a raster stack from the list of monthly rasters
VPD_raster_stack <- rast(VPD_raster_list)

#Add zeros: 
new_colnames <- matrix(1:12, ncol = 12)
new_colnames[1:9] <- paste0("0", 1:9)
names(VPD_raster_stack) <- new_colnames
names(TMAX_raster_stack) <- new_colnames
names(CWD_raster_stack) <- new_colnames

#Make covariates df
#Covariate_df<-data.frame()

for (fire in Fire_List) {
  #Load DOB_RP (0.1% random points)
  setwd("C:/Users/matth/Desktop/Thesis Dev/DOB_Random_Points_01P")
  DOB_RP<-st_read(fire)
  DOB_RP$Month<-as.numeric((str_sub_all(fire,18 ,-15)))
  
  #Use lapply to iterate over every month:
  lapply(unique(DOB_RP$Month), function(x){
    DOB_RP[DOB_RP$Month == x, "cwd_raw"]  <<- terra::extract(CWD_raster_stack[[x]], DOB_RP[DOB_RP$Month == x,])[,2]
    DOB_RP[DOB_RP$Month == x, "tmax_raw"]  <<- terra::extract(TMAX_raster_stack[[x]], DOB_RP[DOB_RP$Month == x,])[,2]
    DOB_RP[DOB_RP$Month == x, "vpd_raw"]  <<- terra::extract(VPD_raster_stack[[x]], DOB_RP[DOB_RP$Month == x,])[,2]
  })
  
  #Remove geom and DOB:
  DOB_RP_df<-as.data.frame(DOB_RP)
  DOB_RP_df<-DOB_RP_df[, -1]
  DOB_RP_df<-DOB_RP_df[, -1]
  #Create means df: 
  DOB_RP_Mean<- DOB_RP_df %>%
    group_by(Month) %>%
    summarise_all(mean)
  
  #Make each 3 decimal places:
  DOB_RP_Mean <- round(DOB_RP_Mean, digits = 3)
  
  #Add DOB_ID:
  Fire_ID<-as.character(str_sub_all(fire,0 ,-13))
  DOB_RP_Mean$Fire_ID<-Fire_ID
  
  
  
  Covariate_df<-rbind(Covariate_df, DOB_RP_Mean)
  
}  

#2020 Done

#Save Final Covariate_df
write.csv(Covariate_df, "C:\\Users\\matth\\Desktop\\Data Dungeon\\Monthly_Covariate_DF_01P.csv")


###Join climate, weather, and topographic covariates to master dataframe: 

#Load master dataframe: 
setwd("C:/Users/matth/Desktop/Thesis Dev/Master_DF")
Master_DF<-read.csv("Master_df.csv")

#Load Grouped_Percent_Landcover dataframe:
setwd("C:/Users/matth/Desktop/Data Dungeon/DOB_Products")
GPL_DF<-read.csv("GPL_by_DOB_PB10.csv")

#Load daily covariates: 
setwd("C:/Users/matth/Desktop/Data Dungeon")
Daily_Covariates<-read.csv("Daily_Covariate_DF_01P.csv")

#Load Monthly covariates: 
setwd("C:/Users/matth/Desktop/Data Dungeon")
Monthly_Covariates<-read.csv("Monthly_Covariate_DF_01P.csv")

#Master first: 
#Spread
GPL_DF$DOB_Area_m2<-Master_DF$Area_m2[match(GPL_DF$DOB_ID, Master_DF$DOB_ID)]
GPL_DF$DOB_Avg_DS<-Master_DF$Avg_DS[match(GPL_DF$DOB_ID, Master_DF$DOB_ID)]
GPL_DF$DOB_Max_DS<-Master_DF$Max_DS[match(GPL_DF$DOB_ID, Master_DF$DOB_ID)]
#Daily
GPL_DF$bui_raw<-Daily_Covariates$bui_raw[match(GPL_DF$DOB_ID, Daily_Covariates$DOB_ID)]
GPL_DF$dsr_raw<-Daily_Covariates$dsr_raw[match(GPL_DF$DOB_ID, Daily_Covariates$DOB_ID)]
GPL_DF$dc_raw<-Daily_Covariates$dc_raw[match(GPL_DF$DOB_ID, Daily_Covariates$DOB_ID)]
GPL_DF$dmc_raw<-Daily_Covariates$dmc_raw[match(GPL_DF$DOB_ID, Daily_Covariates$DOB_ID)]
GPL_DF$ffm_raw<-Daily_Covariates$ffm_raw[match(GPL_DF$DOB_ID, Daily_Covariates$DOB_ID)]
GPL_DF$fwi_raw<-Daily_Covariates$fwi_raw[match(GPL_DF$DOB_ID, Daily_Covariates$DOB_ID)]
GPL_DF$isi_raw<-Daily_Covariates$isi_raw[match(GPL_DF$DOB_ID, Daily_Covariates$DOB_ID)]
GPL_DF$tmax_raw<-Daily_Covariates$tmax_raw[match(GPL_DF$DOB_ID, Daily_Covariates$DOB_ID)]
GPL_DF$prec_raw<-Daily_Covariates$prec_raw[match(GPL_DF$DOB_ID, Daily_Covariates$DOB_ID)]
GPL_DF$ws_raw<-Daily_Covariates$ws_raw[match(GPL_DF$DOB_ID, Daily_Covariates$DOB_ID)]
GPL_DF$rh_raw<-Daily_Covariates$rh_raw[match(GPL_DF$DOB_ID, Daily_Covariates$DOB_ID)]
GPL_DF$Asp_raw<-Daily_Covariates$Asp_raw[match(GPL_DF$DOB_ID, Daily_Covariates$DOB_ID)]
GPL_DF$Elev_raw<-Daily_Covariates$Elev_raw[match(GPL_DF$DOB_ID, Daily_Covariates$DOB_ID)]
GPL_DF$Slp_raw<-Daily_Covariates$Slp_raw[match(GPL_DF$DOB_ID, Daily_Covariates$DOB_ID)]
GPL_DF$Rugg_raw<-Daily_Covariates$Rugg_raw[match(GPL_DF$DOB_ID, Daily_Covariates$DOB_ID)]
GPL_DF$Northness<-Daily_Covariates$Northness[match(GPL_DF$DOB_ID, Daily_Covariates$DOB_ID)]
GPL_DF$Eastness<-Daily_Covariates$Eastness[match(GPL_DF$DOB_ID, Daily_Covariates$DOB_ID)]
#Monthly
GPL_DF$M_cwd_raw<-Monthly_Covariates$cwd_raw[match(GPL_DF$Fire_ID, Monthly_Covariates$Fire_ID)]
GPL_DF$M_tmax_raw<-Monthly_Covariates$tmax_raw[match(GPL_DF$Fire_ID, Monthly_Covariates$Fire_ID)]
GPL_DF$M_vpd_raw<-Monthly_Covariates$vpd_raw[match(GPL_DF$Fire_ID, Monthly_Covariates$Fire_ID)]

#Add DOY:
GPL_DF$DOY<-as.numeric(str_sub_all(GPL_DF$DOB_ID,23 ,-1))


#Save:
write.csv(Master_DF, "C:\\Users\\matth\\Desktop\\Data Dungeon\\GPLXCovariate_df_01P.csv")

#Begin Analysis (Results)

#Load data:
setwd("C:/Users/matth/Desktop/Thesis Dev/Final Data For Publication")
ALL_LC_df<-read.csv("Harris_Aspen_Fire_Spread_Data.csv")#Formerly "GPLXCovariate_df_01P.csv"

unique(ALL_LC_df$Group_Class)

#Create aspen only df: 
Aspen_DOB<-ALL_LC_df[ALL_LC_df$Group_Class=="aspen", ]

#Circularize DOY:
# convert doy to radians
Aspen_DOB$day.of.year.radians <- (Aspen_DOB$DOY - 1) * (2 * pi / 365)
# sin() transform
Aspen_DOB$sin.doy <- sin(Aspen_DOB$day.of.year.radians)


#Base model:
model1<-glmmTMB(log10(DOB_Area_m2)~GC_Percent +(1|Event_ID), family=gaussian, data=Aspen_DOB)
summary(model1)#AIC 2397.8 

#Additive term model:
model2<-glmmTMB(log10(DOB_Area_m2)~GC_Percent + I(sin.doy) + (1|Event_ID), family=gaussian, data=Aspen_DOB)
summary(model2)#AIC 2370.7

#Two dimension model:
model3<-glmmTMB(log10(DOB_Area_m2)~GC_Percent * M_cwd_raw +(1|Event_ID), family=gaussian, data=Aspen_DOB)
summary(model3)#AIC 2391.7

#Two dimension model with additional term:
model4<-glmmTMB(log10(DOB_Area_m2)~GC_Percent * M_cwd_raw + rh_raw +(1|Event_ID), family=gaussian, data=Aspen_DOB)
summary(model4)#AIC 2325.3

#Two dimension model with two interactions: potential bui_raw,dmc_raw,prec_raw, M_cwd_raw(?)
model5<-glmmTMB(log10(DOB_Area_m2)~GC_Percent * Slp_raw + I(sin.doy) * M_cwd_raw  +(1|Event_ID), family=gaussian, data=Aspen_DOB)
summary(model5)#AIC 2355.4  

#Two dimension model with two interactions: potential bui_raw,dmc_raw,prec_raw, M_cwd_raw(?), 
model6<-glmmTMB(log10(DOB_Area_m2)~GC_Percent * Slp_raw +  prec_raw * I(sin.doy)  +(1|Event_ID), family=gaussian, data=Aspen_DOB)
summary(model6)#AIC 2346.2


#All additive: bui_raw + dsr_raw + dc_raw	+ dmc_raw	+ ffm_raw	+ fwi_raw	+ isi_raw	+ tmax_raw + prec_raw	+ ws_raw + rh_raw	+ Elev_raw + Slp_raw	+ Rugg_raw	+ Northness +	Eastness	+ M_cwd_raw	+ M_tmax_raw + sin.doy

#Additive term model:
modelA<-glmmTMB(log10(DOB_Area_m2)~GC_Percent +  bui_raw + dsr_raw + dc_raw	+ dmc_raw	+ ffm_raw	+ fwi_raw	+ isi_raw	+ tmax_raw + prec_raw	+ ws_raw + rh_raw	 + Slp_raw	+ Rugg_raw + Elev_raw	+ Northness	+ Eastness	+ M_cwd_raw	+ M_tmax_raw + sin.doy + (1|Event_ID), family=gaussian, data=Aspen_DOB)
summary(modelA)

#All significant additive:
#Additive term model:
modelAS<-glmmTMB(log10(DOB_Area_m2)~GC_Percent + fwi_raw +  sin.doy + (1|Event_ID), family=gaussian, data=Aspen_DOB)
summary(modelAS)#AIC 2172.8

#2207 AIC for fwi individual effect

####################
#Max Spread Based:

#Base model:
modelB<-glmmTMB(log10(DOB_Max_DS)~GC_Percent +(1|Event_ID), family=gaussian, data=Aspen_DOB)
summary(modelB)#AIC 807.8

#All additive:
#Additive term model:
modelA<-glmmTMB(log10(DOB_Max_DS)~GC_Percent + bui_raw + dsr_raw + dc_raw	+ dmc_raw	+ ffm_raw	+ fwi_raw	+ isi_raw	+ tmax_raw + prec_raw	+ ws_raw + rh_raw	+ Elev_raw + Slp_raw	+ Rugg_raw	+ Northness +	Eastness	+ M_cwd_raw	+ M_tmax_raw + sin.doy + (1|Event_ID), family=gaussian, data=Aspen_DOB)
summary(modelA)

#All significant additive:
#Additive term model:
modelAS<-glmmTMB(log10(DOB_Max_DS)~GC_Percent + fwi_raw  + Slp_raw + sin.doy + (1|Event_ID), family=gaussian, data=Aspen_DOB)
summary(modelAS)#AIC 669.3

#Two dimension model:
model2D<-glmmTMB(log10(DOB_Max_DS)~GC_Percent * Slp_raw + sin.doy +(1|Event_ID), family=gaussian, data=Aspen_DOB)
summary(model2D)#AIC 764.3

#Find r^2
R2 <- r.squaredGLMM(modelAS)

#Make correlation matrix:

# Subset the data to select only columns 8 through 26
subset_data <- Aspen_DOB[, 8:26]
subset_data$sin.doy<-Aspen_DOB$sin.doy
subset_data <- na.omit(subset_data)
# Calculate the correlation matrix
cor_matrix <- cor(subset_data)

# Extract the R-squared values
r_squared <- cor_matrix^2

rounded_r_squared <- round(r_squared, digits = 3)

# Print the R-squared matrix
print(rounded_r_squared)

###Make figures:

#Create aspen only df: 
Aspen_DOB<-ALL_LC_df[ALL_LC_df$Group_Class=="aspen", ]

hist(Aspen_DOB$GC_Percent)

#Make test subset: 25% = 135DOBs, 20% = 186DOBs, 15% = 269DOBs
A10<-Aspen_DOB[Aspen_DOB$GC_Percent>10, ]
A25<-Aspen_DOB[Aspen_DOB$GC_Percent>25, ]

#Base model:
model1<-glmmTMB(log10(DOB_Area_m2)~GC_Percent +(1|Event_ID), family=gaussian, data=Aspen_DOB)
summary(model1)

#Subset model:
model10<-glmmTMB(log10(DOB_Area_m2)~GC_Percent +(1|Event_ID), family=gaussian, data=A10)
summary(model10)

#Subset model:
model25<-glmmTMB(log10(DOB_Area_m2)~GC_Percent +(1|Event_ID), family=gaussian, data=A25)
summary(model25)

#Find r^2 
R2 <- r.squaredGLMM(model25)


eff_data <- as.data.frame(Effect(mod = model1,
                                 term = "GC_Percent",
                                 focal.predictors = "GC_Percent",
                                 xlevels = 100))
# Back-transform DOB_Area
eff_data$fit <- 10^eff_data$fit/10000
eff_data$lower <- 10^eff_data$lower/10000
eff_data$upper <- 10^eff_data$upper/10000


E10 <- as.data.frame(Effect(mod = model10,
                            term = "GC_Percent",
                            focal.predictors = "GC_Percent",
                            xlevels = 100))
# Back-transform DOB_Area
E10$fit <- 10^E10$fit/10000
E10$lower <- 10^E10$lower/10000
E10$upper <- 10^E10$upper/10000

E25 <- as.data.frame(Effect(mod = model25,
                            term = "GC_Percent",
                            focal.predictors = "GC_Percent",
                            xlevels = 100))
# Back-transform DOB_Area
E25$fit <- 10^E25$fit/10000
E25$lower <- 10^E25$lower/10000
E25$upper <- 10^E25$upper/10000



# Daily Area Burned 
DABW <- ggplot() +
  geom_line(data = eff_data, aes(x = GC_Percent, y = fit, color = "Aspen Present"), linewidth = 1.2, lineend = "round") +
  geom_ribbon(data = eff_data, aes(x = GC_Percent, ymin = lower, ymax = upper, fill = "Aspen Present 95% CI"), alpha = 0.4) +
  geom_line(data = E10, aes(x = GC_Percent, y = fit, color = ">10% Aspen"), linewidth = 1.2, lineend = "round", linetype = "longdash") +
  geom_ribbon(data = E10, aes(x = GC_Percent, ymin = lower, ymax = upper, fill = ">10% Aspen 95% CI"), alpha = 0.2) +
  geom_line(data = E25, aes(x = GC_Percent, y = fit, color = ">25% Aspen"), linewidth = 1.2, lineend = "round", linetype = "dotted") +
  geom_ribbon(data = E25, aes(x = GC_Percent, ymin = lower, ymax = upper, fill = ">25% Aspen 95% CI"), alpha = 0.2) +
  labs(x = "Aspen Cover (%)",
       y = "Daily Area Burned (ha)",
       color = "",
       fill = "") +
  scale_x_continuous(limits = c(0, 100)) +
  scale_color_manual(values = c("Aspen Present" = "#00A6B6", ">10% Aspen" = "#7B02A8", ">25% Aspen" = "#F1731D")) +
  scale_fill_manual(values = c("Aspen Present 95% CI" = "#00A6B6", ">10% Aspen 95% CI" = "#7B02A8", ">25% Aspen 95% CI" = "#F1731D")) +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 0.5), 
        axis.text.y = element_text(hjust = 0.5),
        legend.position = c(0.95, 0.95),  # Position the legend in the top right
        legend.justification = c("right", "top"),  # Justify the legend to the top right
        legend.box.just = "right",  # Align legend box to the right
        legend.margin = margin(t = 0.1, r = 0, b = 0, l = 0), # Adjust bottom margin to move legend downwards
        legend.box.margin = margin(t = 1, r = 1, b = 1, l = 1),
        legend.text = element_text(size = 20, face="bold"),
        axis.text = element_text(face="bold", size = 18, color = "black"), 
        axis.title = element_text(face="bold", size = 22, color = "black")
  )
DABW

#Make linear spread figure:

#Base model:
model1<-glmmTMB(log10(DOB_Max_DS)~GC_Percent +(1|Event_ID), family=gaussian, data=Aspen_DOB)
summary(model1)

#Subset model:
model10<-glmmTMB(log10(DOB_Max_DS)~GC_Percent +(1|Event_ID), family=gaussian, data=A10)
summary(model10)

#Subset model:
model25<-glmmTMB(log10(DOB_Max_DS)~GC_Percent +(1|Event_ID), family=gaussian, data=A25)
summary(model25)

#Find r^2 
R2 <- r.squaredGLMM(model25)


eff_data <- as.data.frame(Effect(mod = model1,
                                 term = "GC_Percent",
                                 focal.predictors = "GC_Percent",
                                 xlevels = 100))
# Back-transform DOB_Area
eff_data$fit <- 10^eff_data$fit
eff_data$lower <- 10^eff_data$lower
eff_data$upper <- 10^eff_data$upper


E10 <- as.data.frame(Effect(mod = model10,
                            term = "GC_Percent",
                            focal.predictors = "GC_Percent",
                            xlevels = 100))
# Back-transform DOB_Area
E10$fit <- 10^E10$fit
E10$lower <- 10^E10$lower
E10$upper <- 10^E10$upper

E25 <- as.data.frame(Effect(mod = model25,
                            term = "GC_Percent",
                            focal.predictors = "GC_Percent",
                            xlevels = 100))
# Back-transform DOB_Area
E25$fit <- 10^E25$fit
E25$lower <- 10^E25$lower
E25$upper <- 10^E25$upper



# Daily Area Burned 
DLSW <- ggplot() +
  geom_line(data = eff_data, aes(x = GC_Percent, y = fit, color = "Aspen Present"), linewidth = 1.2, lineend = "round") +
  geom_ribbon(data = eff_data, aes(x = GC_Percent, ymin = lower, ymax = upper, fill = "Aspen Present 95% CI"), alpha = 0.4) +
  geom_line(data = E10, aes(x = GC_Percent, y = fit, color = ">10% Aspen"), linewidth = 1.2, lineend = "round", linetype = "longdash") +
  geom_ribbon(data = E10, aes(x = GC_Percent, ymin = lower, ymax = upper, fill = ">10% Aspen 95% CI"), alpha = 0.2) +
  geom_line(data = E25, aes(x = GC_Percent, y = fit, color = ">25% Aspen"), linewidth = 1.2, lineend = "round", linetype = "dotted") +
  geom_ribbon(data = E25, aes(x = GC_Percent, ymin = lower, ymax = upper, fill = ">25% Aspen 95% CI"), alpha = 0.2) +
  labs(x = "Aspen Cover (%)",
       y = "Daily Maximum Linear Spread (m)",
       color = "",
       fill = "") +
  scale_x_continuous(limits = c(0, 100)) +
  scale_color_manual(values = c("Aspen Present" = "#00A6B6", ">10% Aspen" = "#7B02A8", ">25% Aspen" = "#F1731D")) +
  scale_fill_manual(values = c("Aspen Present 95% CI" = "#00A6B6", ">10% Aspen 95% CI" = "#7B02A8", ">25% Aspen 95% CI" = "#F1731D")) +
  theme_classic() +
  theme(axis.text.x = element_text(hjust = 0.5), 
        axis.text.y = element_text(hjust = 0.5),
        legend.position = c(0.95, 0.95),  # Position the legend in the top right
        legend.justification = c("right", "top"),  # Justify the legend to the top right
        legend.box.just = "right",  # Align legend box to the right
        legend.margin = margin(t = 0.1, r = 0, b = 0, l = 0), # Adjust bottom margin to move legend downwards
        legend.box.margin = margin(t = 1, r = 1, b = 1, l = 1),
        legend.text = element_text(size = 20, face="bold"),
        axis.text = element_text(face="bold", size = 18, color = "black"), 
        axis.title = element_text(face="bold", size = 22, color = "black")
  )
DLSW

#Just DAB:
ggsave("C:/Users/matth/Desktop/Thesis Dev/Figures/DAB_Multi.tiff",
       plot = DAB,
       width = 8,
       height = 8,
       units = "in",
       dpi = 600)

#Make multipanel of battlement
ROS.multipanel <- ggarrange(DABW,
                            DLSW,
                            nrow = 1, ncol = 2,
                            labels = c("A)", "B)"),
                            common.legend = F) 
ROS.multipanel


ggsave("C:/Users/matth/Desktop/Thesis Dev/Figures/ROS.multipanel.tiff",
       plot = ROS.multipanel,
       width = 16,
       height = 8,
       units = "in",
       dpi = 600)


