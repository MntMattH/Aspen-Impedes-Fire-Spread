#Aspen impedes wildfire spread perimeter analysis methods and results:


#Load Packages:
library(landscapemetrics)
library(sf)
library(sp)
library(stringr)
library(terra)
library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)
library(mgcv)
library(MuMIn)

###Methods:

#1) Load Landfire Datasets:
#2001
setwd("C:/Users/matth/Desktop/Data Dungeon/Prior Burn")
LF_2001<-rast("LF_2001_PB10.tif")
#2012
setwd("C:/Users/matth/Desktop/Data Dungeon/Prior Burn")
LF_2012<-rast("LF_2012_PB10.tif")
#2014
setwd("C:/Users/matth/Desktop/Data Dungeon/Prior Burn")
LF_2014<-rast("LF_2014_PB10.tif")
#2016
setwd("C:/Users/matth/Desktop/Data Dungeon/Prior Burn")
LF_2016<-rast("LF_2016_PB10.tif")
#2020
setwd("C:/Users/matth/Desktop/Data Dungeon/Prior Burn")
LF_2020<-rast("LF_2020_PB10.tif")


#2) Begin loop to load in fire perimeters:
setwd("C:/Users/matth/Desktop/Data Dungeon/Perimeters")
Perimeters<-vect("Final_Perimeters.gpkg")

#3 Load landcover classification:  
setwd("C:/Users/matth/Desktop/Data Dungeon/Landfire")
LF_Classes_V3<- read.csv("LF_Classes_V3_PB10.csv")

# Iterate over unique Event_ID values
unique_ids <- unique(Perimeters$Event_ID)
#Make a starting ID
start_id<-"CO3926710806720010417"
start_index <- which(unique_ids == start_id)


if (length(start_index) > 0) {
  # Start the loop from the specified fire ID
  for (i in start_index:length(unique_ids)) {
    event_id <- unique_ids[i]
    
    # Subset the perimeter for the current Event_ID
    subset_Perimeter <- Perimeters[Perimeters$Event_ID == event_id, ]
    plot(subset_Perimeter)
    
    #Create df to assign correct years LC to by fire perimeter year: 
    dfLF<- data.frame(year=seq(2000,2020,1), range= c(rep("2001", 12),rep("2012", 2), rep("2014", 2), rep("2016", 4), rep("2020", 1)))
    #Get Fire Year
    FireYear<-str_sub_all(subset_Perimeter$Event_ID,14,-5)
    FireYear<-as.numeric(FireYear)
    #Define Landcover by range
    range= dfLF$range[which(dfLF$year==FireYear)]
    Landcover<-paste0("LF_", range)
    #select landcover:
    selected_landcover <- switch(Landcover,
                                 "LF_2001" = LF_2001,
                                 "LF_2012" = LF_2012,
                                 "LF_2014" = LF_2014,
                                 "LF_2016" = LF_2016,
                                 "LF_2020" = LF_2020,
                                 default = NULL
    )
    
    #NOTE: Originally we made separate buffers that were 60 meters interior and 60 meters exterior, these were later joined to be a single 120 meter buffer around fire perimeters. The 120 meter buffer was used for all analysize depicted in the manuscript. 
    
    #Create perimeter buffers: 
    SP_Line<-as.lines(subset_Perimeter)
    buffer_distance <- 60
    exterior_buffered <- terra::buffer(subset_Perimeter, width = buffer_distance)
    interior_buffer<-terra::buffer(subset_Perimeter, width = -buffer_distance)
    plot(subset_Perimeter, add=T)
    plot(exterior_buffered, add=T)
    plot(interior_buffer, add=T)
    
    #Clip correct LC to exterior perimeter buffer: 
    clipped_landcover <- crop(selected_landcover, exterior_buffered)
    plot(clipped_landcover)
    
    #Make interior raster vegetation ratio: 
    Interior_landcover <- crop(selected_landcover, interior_buffer)
    Interior_landcover <- mask(clipped_landcover, interior_buffer)
    plot(Interior_landcover)
    Int_Df<-as.data.frame(Interior_landcover)
    
    #Count then Convert: 
    Count<-table(Int_Df$EVT_NAME)
    Interior<-as.data.frame(Count)
    Full_Count<-sum(Interior$Freq)
    Interior$Percent<-Interior$Freq/Full_Count*100
    Interior$Group_Class<-NA
    Interior$Group_Class <-  LF_Classes_V3$Group_Classes[match(Interior$Var1, LF_Classes_V3$VALUE)]
    Interior$EVT<-NA
    Interior$EVT <-  LF_Classes_V3$SAF_SRM[match(Interior$Var1, LF_Classes_V3$VALUE)]
    Interior<- Interior[Interior[, "Freq"]!=0, ]
    
    
    #Make perimeters LC Counts: 
    #Exterior
    ext_clipped <- mask(clipped_landcover, exterior_buffered)
    ext_buff_rast<-mask(ext_clipped, SP_Line)
    clipped_landcover <- crop(selected_landcover, SP_Line)
    #Interior
    int_clipped <- mask(clipped_landcover, subset_Perimeter)
    int_buff_rast<-mask(int_clipped, interior_buffer, inverse=T)
    plot(ext_buff_rast)
    plot(int_buff_rast, add=T)
    plot(subset_Perimeter, add=T)
    
    #Make exterior dataframe to count pixels:
    Ext_Buff<-as.data.frame(ext_buff_rast)
    Count<-table(Ext_Buff$EVT_NAME)
    exterior_buffer<-as.data.frame(Count)
    Full_Count<-sum(exterior_buffer$Freq)
    exterior_buffer$Percent<-exterior_buffer$Freq/Full_Count*100
    exterior_buffer$Group_Class<-NA
    exterior_buffer$Group_Class <-  LF_Classes_V3$Group_Classes[match(exterior_buffer$Var1, LF_Classes_V3$VALUE)]
    exterior_buffer$EVT<-NA
    exterior_buffer$EVT <-  LF_Classes_V3$SAF_SRM[match(exterior_buffer$Var1, LF_Classes_V3$VALUE)]
    exterior_buffer<- exterior_buffer[exterior_buffer[, "Freq"]!=0, ]
    
    #Make interior dataframe to count pixels:
    Int_Buff<-as.data.frame(int_buff_rast)
    Count<-table(Int_Buff$EVT_NAME)
    interior_buffer<-as.data.frame(Count)
    Full_Count<-sum(interior_buffer$Freq)
    interior_buffer$Percent<-interior_buffer$Freq/Full_Count*100
    interior_buffer$Group_Class<-NA
    interior_buffer$Group_Class <- LF_Classes_V3$Group_Classes[match(interior_buffer$Var1, LF_Classes_V3$VALUE)]
    interior_buffer$EVT<-NA
    interior_buffer$EVT <- LF_Classes_V3$SAF_SRM[match(interior_buffer$Var1, LF_Classes_V3$VALUE)]
    interior_buffer<- interior_buffer[interior_buffer[, "Freq"]!=0, ]
    
    # Finalize:
    #Add fire ID
    Interior$Event_id <- event_id
    exterior_buffer$Event_id <- event_id
    interior_buffer$Event_id <- event_id
    
    #Add fire size to all patches:
    subset_Perimeter$Area<-expanse(subset_Perimeter, unit = "ha")
    
    Interior$Fire.Area <- subset_Perimeter$Area
    exterior_buffer$Fire.Area <- subset_Perimeter$Area
    interior_buffer$Fire.Area <- subset_Perimeter$Area
    
    #Prepare Vectors: 
    Interior_landcover<-as.polygons(Interior_landcover)
    ext_buff_vect<-as.polygons(ext_buff_rast)
    int_buff_vect<-as.polygons(int_buff_rast)
    
    #Add all variables: 
    Interior_landcover <- merge(Interior_landcover, 
                                Interior, 
                                by.x = "EVT_NAME",
                                by.y = "Var1",
                                all = TRUE) 
    
    #Ext: 
    ext_buff_vect <- merge(ext_buff_vect, 
                           exterior_buffer, 
                           by.x = "EVT_NAME",
                           by.y = "Var1",
                           all = TRUE)
    #Int: 
    int_buff_vect <- merge(int_buff_vect, 
                           interior_buffer,
                           by.x = "EVT_NAME",
                           by.y = "Var1",
                           all = TRUE)
    
    #9) Save the dataframes and rasters per fire:
    output_dir<-"C:/Users/matth/Desktop/Data Dungeon/Perimeters/PB_Interior_LCR"
    output_dir2<-"C:/Users/matth/Desktop/Data Dungeon/Perimeters/PB_Exterior_Buffer_LCR"
    output_dir3<-"C:/Users/matth/Desktop/Data Dungeon/Perimeters/PB_Interior_Buffer_LCR"
    #Interior:
    output_file <- file.path(output_dir, paste0(basename(event_id), "_Interior_LCR_PB10.gpkg"))
    writeVector(Interior_landcover, output_file, overwrite=T)
    #Ext_Buff:
    output_file2 <- file.path(output_dir2, paste0(basename(event_id), "_Exterior_Buffer_LCR_PB10.gpkg"))
    writeVector(ext_buff_vect, output_file2, overwrite=T)
    #Int_Buff
    output_file3 <- file.path(output_dir3, paste0(basename(event_id), "_Interior_Buffer_LCR_PB10.gpkg"))
    writeVector(int_buff_vect, output_file3, overwrite=T)
    
  }
} else {
  cat("Start ID not found in the list.")
}

###Use save perimeters and interiors to create a dataframe comprising all fires: 

#First create dataframe with all interior LC:
folder_path <- "C:/Users/matth/Desktop/Data Dungeon/Perimeters/PB10_Interior_LCR" 

spatial_vectors <- lapply(list.files(folder_path, full.names = TRUE), st_read)

Interior_LCR_df <- bind_rows(spatial_vectors)

rm(spatial_vectors)

Interior_LCR_DF<-as.data.frame(Interior_LCR_df)
Interior_LCR_DF<-subset(Interior_LCR_DF, select = -geom)

Interior_LCR_DF$CLASSNAME <- ifelse(is.na(Interior_LCR_DF$CLASSNAME), "", Interior_LCR_DF$CLASSNAME)
Interior_LCR_DF$EVT_NAME <- ifelse(is.na(Interior_LCR_DF$EVT_NAME), "", Interior_LCR_DF$EVT_NAME)
Interior_LCR_DF$OID_ <- ifelse(is.na(Interior_LCR_DF$OID_), "", Interior_LCR_DF$OID_)

Interior_LCR_DF$EVT_Code<-paste0(Interior_LCR_DF$CLASSNAME, Interior_LCR_DF$EVT_NAME, Interior_LCR_DF$OID_)
Interior_LCR_Full<-subset(Interior_LCR_DF, select = -c(CLASSNAME, EVT_NAME, OID_))

#Create full exterior buffer DF:
folder_path <- "C:/Users/matth/Desktop/Data Dungeon/Perimeters/PB10_Exterior_Buffer_LCR" 

spatial_vectors <- lapply(list.files(folder_path, full.names = TRUE), st_read)

Exterior_Buffer_LCR_df <- bind_rows(spatial_vectors)

rm(spatial_vectors)

Exterior_Buffer_LCR_DF<-as.data.frame(Exterior_Buffer_LCR_df)
Exterior_Buffer_LCR_DF<-subset(Exterior_Buffer_LCR_DF, select = -geom)

Exterior_Buffer_LCR_DF$CLASSNAME <- ifelse(is.na(Exterior_Buffer_LCR_DF$CLASSNAME), "", Exterior_Buffer_LCR_DF$CLASSNAME)
Exterior_Buffer_LCR_DF$EVT_NAME <- ifelse(is.na(Exterior_Buffer_LCR_DF$EVT_NAME), "", Exterior_Buffer_LCR_DF$EVT_NAME)
Exterior_Buffer_LCR_DF$OID_ <- ifelse(is.na(Exterior_Buffer_LCR_DF$OID_), "", Exterior_Buffer_LCR_DF$OID_)

Exterior_Buffer_LCR_DF$EVT_Code<-paste0(Exterior_Buffer_LCR_DF$CLASSNAME, Exterior_Buffer_LCR_DF$EVT_NAME, Exterior_Buffer_LCR_DF$OID_)
Exterior_Buffer_LCR_Full<-subset(Exterior_Buffer_LCR_DF, select = -c(CLASSNAME, EVT_NAME, OID_))
#Create full interior buffer DF:
folder_path <- "C:/Users/matth/Desktop/Data Dungeon/Perimeters/PB10_Interior_Buffer_LCR" 

spatial_vectors <- lapply(list.files(folder_path, full.names = TRUE), st_read)

Interior_Buffer_LCR_df <- bind_rows(spatial_vectors)

rm(spatial_vectors)

Interior_Buffer_LCR_DF<-as.data.frame(Interior_Buffer_LCR_df)
Interior_Buffer_LCR_DF<-subset(Interior_Buffer_LCR_DF, select = -geom)

Interior_Buffer_LCR_DF$CLASSNAME <- ifelse(is.na(Interior_Buffer_LCR_DF$CLASSNAME), "", Interior_Buffer_LCR_DF$CLASSNAME)
Interior_Buffer_LCR_DF$EVT_NAME <- ifelse(is.na(Interior_Buffer_LCR_DF$EVT_NAME), "", Interior_Buffer_LCR_DF$EVT_NAME)
Interior_Buffer_LCR_DF$OID_ <- ifelse(is.na(Interior_Buffer_LCR_DF$OID_), "", Interior_Buffer_LCR_DF$OID_)

Interior_Buffer_LCR_DF$EVT_Code<-paste0(Interior_Buffer_LCR_DF$CLASSNAME, Interior_Buffer_LCR_DF$EVT_NAME, Interior_Buffer_LCR_DF$OID_)
Interior_Buffer_LCR_Full<-subset(Interior_Buffer_LCR_DF, select = -c(CLASSNAME, EVT_NAME, OID_))


#Save all compiled dataframes: 
#Interior
fwrite(Interior_LCR_Full, "C:\\Users\\matth\\Desktop\\Data Dungeon\\Perimeters\\PB10_Interior_LCR_Full.csv", append= FALSE)
#Ext_Buff
fwrite(Exterior_Buffer_LCR_Full, "C:\\Users\\matth\\Desktop\\Data Dungeon\\Perimeters\\PB10_Exterior_Buffer_LCR_Full.csv", append= FALSE)
#Int_Buff
fwrite(Interior_Buffer_LCR_Full, "C:\\Users\\matth\\Desktop\\Data Dungeon\\Perimeters\\PB10_Interior_Buffer_LCR_Full.csv", append= FALSE)


#Begin Analysis (Results)

#Load data
setwd("C:/Users/matth/Desktop/Data Dungeon/Perimeters")
Interior<-read.csv("PB10_Interior_LCR_Full.csv")
Ext_Buffer<-read.csv("PB10_Exterior_Buffer_LCR_Full.csv")
Int_Buffer<-read.csv("PB10_Interior_Buffer_LCR_Full.csv")

###The above dataframes have not been provided as the complete dataframe ("Harris_Aspen_Perimeter_Analysis_Data") has all of this information already compiled.  

#Run once the LC_Ratio_DF has been assembled to subset for 0.1% aspen in ther perimeter
Interior <- Interior[Interior$Event_id %in% A01_Firelist, ]
Ext_Buffer <- Ext_Buffer[Ext_Buffer$Event_id %in% A01_Firelist, ]
Int_Buffer <- Int_Buffer[Int_Buffer$Event_id %in% A01_Firelist, ]

#Prepare dataframe for sign test (make proportion and ratio for each LC type across all fires)
#Combine to one dataframe: 
Ext_Buffer$Position<-"Perimeter"
Int_Buffer$Position<-"Perimeter"
Interior$Position<-"Interior"

#Make interior Count
Interior_Count <- Interior %>%
  group_by(Group_Class, Event_id, Position) %>%
  summarise(Total_Frequency = sum(Freq))

#Create perimeter (make into one continuous buffer): 
Perim_Df<-rbind(Ext_Buffer, Int_Buffer)
#Make perimeter Counts
Perimeter_Count <- Perim_Df %>%
  group_by(Group_Class, Event_id, Position) %>%
  summarise(Total_Frequency = sum(Freq))

#Create per fire pixel sums:
Perim_Total_Freq<-Perimeter_Count %>%
  group_by(Event_id) %>%
  summarise(PF_Sum = sum(Total_Frequency))

Interior_Total_Freq<-Interior_Count %>%
  group_by(Event_id) %>%
  summarise(PF_Sum = sum(Total_Frequency))


#Make proportions: 
#Add sums to total DF
#Int
Interior_Count$Int_PF_Sum <-  Interior_Total_Freq$PF_Sum[match(Interior_Count$Event_id, Interior_Total_Freq$Event_id)]
Interior_Count$Int_Prop<-Interior_Count$Total_Frequency/Interior_Count$Int_PF_Sum
#Perim
Perimeter_Count$Perim_PF_Sum <-  Perim_Total_Freq$PF_Sum[match(Perimeter_Count$Event_id, Perim_Total_Freq$Event_id)]
Perimeter_Count$Perim_Prop<-Perimeter_Count$Total_Frequency/Perimeter_Count$Perim_PF_Sum

#Now combine into one dataframe to make proportions:

Ratio_Df<-subset(Interior_Count, select = -Position)
Ratio_Df<-subset(Ratio_Df, select = -Total_Frequency)
Ratio_Df<-subset(Ratio_Df, select = -Int_PF_Sum)

#Make LC_ID for matching purposes: 
Ratio_Df$LC_ID<-paste(Ratio_Df$Event_id, Ratio_Df$Group_Class, sep = "_")
Perimeter_Count$LC_ID<-paste(Perimeter_Count$Event_id, Perimeter_Count$Group_Class, sep = "_")

#Now match ratios FPE:
matching_df <- inner_join(Ratio_Df, Perimeter_Count, by = "LC_ID")
matching_df$FPE<-((matching_df$Perim_Prop-matching_df$Int_Prop)/(matching_df$Perim_Prop+matching_df$Int_Prop))
matching_df$FPE<-round(matching_df$FPE, 2)


LC_Ratio_DF <- matching_df[-1, c(1, 2, 3, 10, 11)]#Got it all together!!! 

###For reviewers please load in ("Harris_Aspen_Perimeter_Analysis_Data.csv") at this time:
setwd("C:/Users/matth/Desktop/Thesis Dev/Final Data For Publication")
LC_Ratio_DF<-read.csv("Harris_Aspen_Perimeter_Analysis_Data.csv")

#Subsample for 0.1% aspen in perimeter:
Aspen_LC_Ratio<- subset(LC_Ratio_DF, LC_Ratio_DF$Group_Class=="aspen")
PA01<-Aspen_LC_Ratio[Aspen_LC_Ratio$Perim_Prop>0.001,]
A01_Firelist<-unique(PA01$Event_id)
filtered_df <- LC_Ratio_DF[LC_Ratio_DF$Event_id %in% A01_Firelist, ]

#Go back to the top to use the A01_Firelist to subset this dataframe to 0.1% aspen in the perimeter

#################
#Now run sign test 

#Running: Aspen

#Create aspen only DF: 
Aspen_LC_Ratio<- subset(filtered_df, filtered_df$Group_Class=="aspen")

#Sign test: 
sign_test <- function(id, Aspen_LC_Ratio) {
  # Subset the data for the given ID
  subset_data <- Aspen_LC_Ratio[Aspen_LC_Ratio$Event_id == id, "FPE"]
  
  # Count the number of positive and negative differences from 0
  positive_diff <- sum(subset_data > 0)
  negative_diff <- sum(subset_data < 0)
  
  # Perform the sign test
  sign_test_result <- sign(positive_diff - negative_diff)
  
  return(sign_test_result)
}

# Get unique IDs from the data
unique_ids <- unique(Aspen_LC_Ratio$Event_id)

# Perform sign test for each unique ID
sign_test_results <- sapply(unique_ids, function(id) sign_test(id, Aspen_LC_Ratio))

# Print the results
result_df <- data.frame(ID = unique_ids, SignTestResult = sign_test_results)
print(result_df)

#Sum positives and negatives:
total_positive <- sum(result_df$SignTestResult == 1)
total_negative <- sum(result_df$SignTestResult == -1)

# Perform a sign test on the aggregated data
All_fire_sign_test_result <- sign(total_positive - total_negative)

# Calculate a p-value for the aggregated sign test
p_value <- binom.test(min(total_positive, total_negative), n = total_positive + total_negative, alternative = "two.sided")$p.value
p_value 


#Now inspect the influence of climate (CWD_Zscores) and weather (FWI)

#Load covariate df: Already joined in "Harris_Aspen_Perimeter_Analysis_Data" skip to line 432. 
setwd("C:/Users/matth/Desktop/Data Dungeon/C&W_Perimeter_Analysis")
FWI_raw<-read.csv("Final_PA_Covariate_df_01P.csv")
FWI_ZC<-read.csv("FWI_Z_Score_df.csv")

#Load CWD Z-scores:
setwd("C:/Users/matth/Desktop/Data Dungeon/C&W_Perimeter_Analysis")
CWD_ZC<-read.csv("PA_FireSum_CWD_ZScores_01P.csv")

#Load ISI Z-scores:
setwd("C:/Users/matth/Desktop/Data Dungeon/C&W_Perimeter_Analysis")
ISI_ZC<-read.csv("ISI_Z_Score_df_01P.csv")

#Add day of the year:
day_df<-filtered_df
day_df$date<-as.character(str_sub_all(day_df$Event_id,14 ,-1))
# Convert to a regular date object
day_df$regular_date <- as.Date(day_df$date, format = "%Y%m%d")

day_df$DOY <- as.POSIXlt(day_df$regular_date)$yday + 1


#Match
filtered_df$fwi_fire_mean<-FWI_raw$fwi_raw[match(filtered_df$Event_id, FWI_raw$Fire_ID)]
filtered_df$fwi_ZScore<-FWI_ZC$FWI_ZScore[match(filtered_df$Event_id, FWI_ZC$Fire_ID)]
filtered_df$cwd_ZScore<-CWD_ZC$CWD_ZScore[match(filtered_df$Event_id, CWD_ZC$Fire_ID)]
filtered_df$isi_ZScore<-ISI_ZC$ISI_ZScore[match(filtered_df$Event_id, ISI_ZC$Fire_ID)]
filtered_df$isi_fire_mean<-ISI_ZC$FireMean[match(filtered_df$Event_id, ISI_ZC$Fire_ID)]
filtered_df$DOY<-day_df$DOY[match(filtered_df$Event_id, day_df$Event_id)]

#Subset to aspen: 
Aspen_LC_Ratio<- subset(filtered_df, filtered_df$Group_Class=="aspen")

#Aspen ratio by ISI Z-scores
model1<-lm(FPE~isi_ZScore, data=Aspen_LC_Ratio)
summary(model1)

#Aspen ratio by FWI z-score (fires lifetime)
model2<-lm(FPE~fwi_ZScore, data=Aspen_LC_Ratio)
summary(model2)

#Other covariates tested in model above: fwi_fire_mean, CWD_ZScore, isi_fire_mean, DOY

#Run a GAM to compare LCRs to Day of Year

# convert doy to radians
Aspen_LC_Ratio$day.of.year.radians <- (Aspen_LC_Ratio$DOY - 1) * (2 * pi / 365)
# sin() transform
Aspen_LC_Ratio$sin.doy <- sin(Aspen_LC_Ratio$day.of.year.radians)
# cos() transform 
Aspen_LC_Ratio$cos.doy <- cos(Aspen_LC_Ratio$day.of.year.radians)

# Fit a GAM 
model.sin <- gam(FPE ~ sin.doy, data = Aspen_LC_Ratio)

anova(model.sin) 
summary(model.sin)
R2 <- r.squaredGLMM(model.sin)

# extract model fit and 95% CI
tmp <- predict(model.sin, 
               se.fit = TRUE, # this argument will return the SE at each fit value
               type = "response") 
head(tmp) 
tmp <- data.frame(tmp)
# calculate 95% CI by hand from the SE using critical value from the normal distribution
tmp$lwr <- tmp$fit - tmp$se.fit*qnorm(0.975)
tmp$upr <- tmp$fit + tmp$se.fit*qnorm(0.975)
head(tmp) 

# add DOY back to this tmp dataframe
tmp <- cbind(tmp, 
             Aspen_LC_Ratio[,5|11]) # there are other fit and lwr columns here so I'll use the index to grab the FPE and DOY columns


#Aspen FPE by DOY
DOY<-ggplot(tmp, aes(x = DOY, y = fit)) +
  geom_point(aes(x = DOY, y = FPE), col="#E65D2FFF") +
  geom_ribbon(aes(ymin = lwr, 
                  ymax = upr),
              fill = "lightblue",
              alpha = 0.5) +
  geom_line(col="black", size=1) + # this is the fitted trendline
  labs(x = "Day of Year",
       y = "Aspen Fire Perimeter Effect") +
  scale_x_continuous(breaks = c(0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 330), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec")) +
  theme_classic()+
  geom_hline(yintercept = 0, color = "black", size=1, linetype = "dashed")+
  theme(axis.text.x = element_text(hjust = 0.5), 
        axis.text.y = element_text(hjust = 0.5),
        legend.position = "top",  # Set legend position to the top
        legend.justification = "right",  # Justify legend to the right
        legend.box.just = "right",  # Align legend box to the right
        legend.margin = margin(t = 0.1, r = 0, b = 0, l = 0), # Adjust bottom margin to move legend downwards
        legend.box.margin = margin(t = 1, r = 1, b = 1, l = 1),
        legend.text = element_text(size = 20, face="bold"),
        axis.text = element_text(face="bold", size = 18, color = "black"), 
        axis.title = element_text(face="bold", size = 22, color = "black")
        
  )
DOY

#Make plots: FPE boxplot

#load data: 
setwd("C:/Users/matth/Desktop/Data Dungeon/Perimeters")
FPE<-read.csv("IPR_FPE.csv")


custom_colors3 <- c("#FFCD00","#B3E0A6", "#98D688", "#86CA78",  "#76BD6A", "#68B05D","grey",
                    "#5AA354", "#4B974F", "#3B8A4A", "#2C7E41")


#Reorder
filtered_df <- filtered_df[order(filtered_df$FPE), ]
FPE$Landcover.Class <- factor(FPE$Landcover.Class, levels = FPE$Landcover.Class)

#Fix labels:
format_label <- function(label) {
  str_to_title(gsub("_", " ", label))
}

filtered_df <- filtered_df %>%
  mutate(formatted_labels = sapply(Group_Class, format_label))

#Add Stars and means: 
# Calculate median values for each category
df <- filtered_df %>%
  group_by(Group_Class) %>%
  mutate(median_value = round(median(FPE), 2))

medians <- aggregate(FPE ~ formatted_labels, filtered_df, median)
medians$Stars<-FPE$Stars[match(medians$formatted_labels, FPE$Landcover.Class)]
medians$Means<-FPE$Mean_FPE[match(medians$formatted_labels, FPE$Landcover.Class)]
filtered_df$median_value<-df$median_value[match(filtered_df$Group_Class, df$Group_Class)]

#Make boxplot

BP<-ggplot(filtered_df, aes(x = reorder(formatted_labels, FPE, FUN = median), y = FPE, fill = formatted_labels)) +
  geom_boxplot()+
  labs(x = "",
       y = "Fire Perimeter Effect") +
  theme_classic()+
  theme(axis.text.y = element_text(family = "TT Arial", face="bold", size = 20, color = "black"),
        axis.text.x = element_text(family = "TT Arial", face="bold", angle = 45, hjust = 1, size = 20, color = "black"),
        axis.title = element_text(family = "TT Arial",face="bold", size = 26, color = "black"))+
  scale_fill_manual(values = custom_colors3)+
  geom_hline(yintercept = 0, color = "black")+
  guides(fill = FALSE)+#Remove legend
  geom_text(data = medians, aes(x = formatted_labels, y = ifelse(FPE >= 0, FPE - 0.08, FPE - 0.2), label = Stars), vjust = -1, size = 12)+
  scale_x_discrete(expand = c(0.055, 0.055, 0.18, 0.18))
BP

ggsave("C:/Users/matth/Desktop/Thesis Dev/Figures/BP_5.tiff",
       plot = BP,
       width = 18,
       height = 10,
       units = "in",
       dpi = 600)





