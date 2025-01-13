# load libraries
library("neotoma2")
library("Bchron")
library("tidyr")
library("dplyr")

# create directories
dir.create("TempFiles") # for temporary files
dir.create("IndividualSummaries") # for time calibrated summary tables
dir.create("Results") # for results to be used in paper or supplement

###################################################################################
############################## CALL RADIOMETRIC DATA ############################## 
###################################################################################

# create a function that calls the geochronologic controls for each of the sites in the Neotoma database
neotomaGeochron <- function(site_ids) {
  
  # Create an empty list to store results
  geochronologic_data <- list()
  
  # Loop through the site IDs
  for (site_id in site_ids) {
    # Try to fetch site data
    site_data <- tryCatch({
      site <- neotoma2::get_sites(siteid = site_id)
    }, error = function(e) {
      message(paste("Error fetching data for site ID:", site_id))
      return(NULL)
    })
    
    # If site_data is valid, fetch geochronologic controls
    if (!is.null(site_data)) {
      geo_controls <- tryCatch({
        
        site_geochron <- neotoma2::get_datasets(site, all_data = TRUE) %>%
          neotoma2::filter(datasettype == "geochronologic")
        
        site_geochrondl <- site_geochron %>%
          neotoma2::get_downloads()
        
        site_geochrondl@sites[["site"]]@collunits@collunits[[1]]@chronologies@chronologies[[1]]@chroncontrols
        
      }, error = function(e) {
        message(paste("Error fetching geochronologic controls for site ID:", site_id))
        return(NULL)
      })
      
      # Store the results
      if (!is.null(geo_controls)) {
        geo_controls <- as.data.frame(geo_controls)
        geo_controls <- cbind(siteid = site_id, geo_controls)  # Add siteid as the first column
        geochronologic_data[[length(geochronologic_data) + 1]] <- geo_controls
      }
    }
    
    # pause before looping to avoid overwhelming Neotoma server
    Sys.sleep(2)
  }
  
  # Combine all geochronologic controls into a single data frame if needed
  combined_geo_data <- do.call(rbind, lapply(geochronologic_data, as.data.frame))
}

loc_controls <- neotomaGeochron(site_ids = c(10537, 10539, 513, 2271, 10538, 1396, 1748, 790, 992,  1974, 
                                         2270, 1973, 2245, 1977, 10102, 1955, 207, 2232, 1699, 1503,
                                         1355, 2551, 13690, 11575, 11579, 11583, 846))

# manually call geochronologic controls for Site 11583 (different formatting)
geochron11583 <- neotoma2::get_datasets(siteid = 11583, all_data = TRUE) %>%
  neotoma2::get_downloads()
 geo_controls11583 <- geochron11583@sites[["site"]]@collunits@collunits[[1]]@chronologies@chronologies[[1]]@chroncontrols
geo_controls11583
# since Site 11583 is dated using stratigraphy, we will exclude it from the calibration

# filter dataset to only include cores that are dated using radiocarbon dates
radiocarbon_sites <- unique(loc_controls$siteid[loc_controls$chroncontroltype == "Radiocarbon"])
filtered_controls <- loc_controls %>%
  dplyr::filter(siteid %in% radiocarbon_sites & (chroncontroltype == "Radiocarbon" | chroncontroltype == "Core top"))

# save dataframe as a .csv file for easy recall
write.csv(filtered_controls, "TempFiles/radiocarbonControl.csv", row.names = FALSE)

###################################################################################
###################################################################################
###################################################################################



##############################################################################
############################## CALL POLLEN DATA ##############################
##############################################################################

# function to pull pollen data from Neotoma database
neotomaPollen <- function(site_ids, taxa) {
  
  # create a list to store the results
  pollen_data <- list()
  index <- 1
  
  # Loop through the site IDs
  for (site_id in site_ids) {
  
      site <- neotoma2::get_sites(siteid = site_id)
  
        # find pollen data for site
        site_pollen <- neotoma2::get_datasets(site, all_data = TRUE) %>%
          neotoma2::filter(datasettype == "pollen") %>%
          get_downloads()
        allSamp = samples(site_pollen)
        
        for (taxon in taxa) {
          # harmonize taxa based on user input
          allSamp = allSamp %>% 
            dplyr::filter(ecologicalgroup %in% c("TRSH")) %>% 
            mutate(variablename = replace(variablename, 
                                          stringr::str_detect(variablename, taxon), 
                                          taxon))
          
          # create a function to check and add the specific taxon if not present
          ensure_taxon_present <- function(df, taxon) {
            if (!(taxon %in% df$variablename)) {
              df <- bind_rows(df, data.frame(sitename = df$sitename[1], lat = df$lat[1], long = df$long[1], siteid = df$siteid[1], datasetid = df$datasetid[1], age = df$age[1], variablename = taxon, value = 0, depth = df$depth[1]))
            }
            return(df)
          }
          
          # Apply the function to each group
          allSamp0 = allSamp %>%
            group_by(sitename, lat, long, siteid, datasetid, age, variablename, depth) %>%
            summarize(value = sum(value), .groups = "keep") %>%
            group_by(sitename, lat, long, siteid, datasetid, age, depth) %>%
            do(ensure_taxon_present(., taxon)) %>%
            ungroup() %>%
            dplyr::filter(variablename == taxon) %>%
            select(sitename, lat, long, siteid, datasetid, value, age, variablename, depth)
        
          # Append the result to the list
          if (!is.null(allSamp0)) {
            pollen_data[[index]] <- allSamp0
            index <- index + 1
          }
        }
  }

  # Combine all results after both loops
  combined_pollen_data <- bind_rows(pollen_data)
}

loc_pollen <- neotomaPollen(site_ids = c(10537, 10539, 513, 2271, 10538, 1396, 1748, 790, 992,  1974, 
                                         2270, 1973, 2245, 1977, 10102, 1955, 207, 2232, 1699, 1503,
                                         1355, 2551, 13690, 11575, 11579, 11583, 846),
                            taxa = c("Salix", "Populus", "Picea"))
  
# pivot table to display each taxon as a separate column
pollen_wide <- loc_pollen %>%
  pivot_wider(
    names_from = variablename, 
    values_from = value         
  )

# save dataframe as a .csv file for easy recall
write.csv(pollen_wide, "TempFiles/pollen_wide.csv", row.names= FALSE)

##############################################################################
##############################################################################
##############################################################################



#############################################################################
############################## CALIBRATE DATES ############################## 
#############################################################################

# if needed, read radiometric controls and pollen data from temporary .csv files
filtered_controls <- read.csv("TempFiles/radiocarbonControl.csv")
pollen_wide <- read.csv("TempFiles/pollen_wide.csv")

# create function for creating age-depth model for each site and calibrating dates
calibrateDates <- function(filtered_controls, pollen_wide){
  # create a new column in filtered_controls called age_sd
  filtered_controls$age_sd <- (filtered_controls$agelimitolder - filtered_controls$agelimityounger)/2
  
  # Handle NA age_sd values for core tops (set to 10 years or another reasonable value)
  filtered_controls$age_sd[is.na(filtered_controls$age_sd)] <- 10  # For example, 10 years uncertainty for core tops
  
  # Automatically assign calibration curves for radiocarbon dates
  filtered_controls$cal_curve <- ifelse(filtered_controls$chroncontroltype == "Radiocarbon", "intcal20", NA)
  
  # Initialize output dataframe with the same structure as pollen_wide
  output_df <- pollen_wide[0, ] # Empty dataframe with the same columns as pollen_wide
  
  # Loop through each siteid
  for (site in unique(pollen_wide$siteid)) {
    
    site_controls <- filtered_controls[filtered_controls$siteid == site, ]
    site_depths <- pollen_wide[pollen_wide$siteid == site, ]
    
    # Ensure `calibrated_age` column exists in site_depths
    site_depths$calibrated_age <- NA
    
    # Separate core tops and radiocarbon controls
    radiocarbon_controls <- site_controls[!is.na(site_controls$cal_curve), ]
    core_tops <- site_controls[is.na(site_controls$cal_curve), ]
    
    # Create the age-depth model only if there are radiocarbon controls
    if (nrow(radiocarbon_controls) > 0) {
      age_depth_model <- Bchronology(
        ages = radiocarbon_controls$chroncontrolage,
        ageSds = radiocarbon_controls$age_sd,
        positions = radiocarbon_controls$depth,
        positionThicknesses = radiocarbon_controls$thickness,
        calCurves = radiocarbon_controls$cal_curve
      )
      
      # Predict ages for all depths
      calibrated_dates <- data.frame(predict(age_depth_model, newPositions = site_depths$depth))
      median_predicted_ages <- apply(calibrated_dates, 2, median, na.rm = TRUE)
      site_depths$calibrated_age <- median_predicted_ages
      
    } else {
      # Use default ages if no radiocarbon controls exist
      site_depths$calibrated_age <- site_depths$age
    }
    
    # Append core tops back into the result with fixed ages
    if (nrow(core_tops) > 0) {
      for (i in 1:nrow(core_tops)) {
        depth <- core_tops$depth[i]
        age <- core_tops$chroncontrolage[i]
        site_depths$calibrated_age[site_depths$depth == depth] <- age
      }
    }
    
    # Bind the updated site_depths to the output dataframe
    output_df <- rbind(output_df, site_depths)
  }
  
  return(output_df)
}

# call the calibration function
output_df <- calibrateDates(filtered_controls, pollen_wide)

# move radiometric age ("age") after calibrated age for easier comparison
output_df <- output_df %>%
  relocate(age, .after = calibrated_age)

# # save dataframe as a .csv file for easy recall
write.csv(output_df, "TempFiles/calibratedDates.csv", row.names = FALSE)

#############################################################################
#############################################################################
#############################################################################



###################################################################################
############################## BIN AND ORGANIZE DATA ##############################
###################################################################################

# load libraries
library("neotoma2")
library("Bchron")
library("tidyr")
library("dplyr")

# if needed, read calibrated radiometric data from temporary .csv files
output_df <- read.csv("TempFiles/calibratedDates.csv")

# create function for organizing data into time bins
organizeData <- function(timeBin, taxon, samplingProtocol, yearMin, yearMax) {
  # create 500-yr time bins as a separate column
  timeCorrected = output_df %>%
    dplyr::filter(age >= 0) %>%
    mutate(Year_Bin = floor(calibrated_age / timeBin) * timeBin)
  
  # filter to look at one taxon
  timeCorrected = timeCorrected %>%
    select(sitename, lat, long, siteid, datasetid, depth, all_of(taxon), calibrated_age, Year_Bin) %>%
    rename(value = taxon)
  
  data_filtered = timeCorrected %>%
    group_by(sitename, Year_Bin) %>%
    slice_min(order_by = value, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::filter(Year_Bin >= yearMin) %>%
    dplyr::filter(Year_Bin <= yearMax)
  
  # creates pivot table with correctly ordered time bins
  ordered_years = sort(unique(data_filtered$Year_Bin))
  pivot_table = data_filtered %>%
    select(sitename, siteid, datasetid, lat, long, Year_Bin, value) %>%
    pivot_wider(names_from = Year_Bin, values_from = value, values_fill = list(Taxon_Abundance = NA)) %>%
    select(sitename, siteid, datasetid, lat, long, all_of(as.character(ordered_years)))
  
  # return pivot table
  return(pivot_table)
}

PiceaMin <- organizeData(timeBin = 500,
                         taxon = "Picea",
                         samplingProtocol = "Minimum",
                         yearMin = 0,
                         yearMax = 20000)

write.csv(PiceaMin, "IndividualSummaries/PiceaMin.csv", row.names = FALSE)

###################################################################################
###################################################################################
###################################################################################




















all_depths <- as.numeric(pollen_Samp$depth)
known_depths <- geochron$depth
ages = c(2540, 6620, 9350, 14140)
#errors = c(60, 90, 90, 440)

radiocarbon_data <- data.frame(
  depth = all_depths,
  age = ifelse(all_depths %in% known_depths, ages[match(all_depths, known_depths)], NA),
  error = ifelse(all_depths %in% known_depths, errors[match(all_depths, known_depths)], NA)
)

# Filter for known radiocarbon dates
known_dates <- radiocarbon_data[!is.na(radiocarbon_data$age), ]

# Calibrate radiocarbon dates
calibrated_dates <- BchronCalibrate(
  ages = known_dates$age,
  ageSds = known_dates$error,
  calCurves = rep("intcal20", nrow(known_dates))
)

median_calibrated_vector <- numeric(length(calibrated_dates))

# Access the calibration data for each date
for (i in 1:length(calibrated_dates)) {
  # Access the ageGrid and densities for the current date
  age_grid <- calibrated_dates[[i]]$ageGrid
  densities <- calibrated_dates[[i]]$densities
  
  # Calculate the median calibrated age
  # The median can be found by locating the ageGrid value where the cumulative density reaches 50%
  cumulative_density <- cumsum(densities) / sum(densities)
  median_index <- which.min(abs(cumulative_density - 0.5))
  
  # The median calibrated age is the corresponding value from ageGrid
  median_calibrated_age <- age_grid[median_index]
  median_calibrated_vector[i] <- median_calibrated_age
  
  print(paste("Median calibrated age for Date", i, ":", median_calibrated_age))
}

# Fit the age-depth age_depth_model
age_depth_model <- Bchronology(
  ages = known_dates$age,
  ageSds = known_dates$error,
  calCurves = rep("intcal20", nrow(known_dates)),
  positions = known_dates$depth,  # Depths with radiocarbon dates
  positionThicknesses = rep(1, nrow(known_dates))
)

# Predict ages for all depths
predicted_ages <- data.frame(predict(age_depth_model, newPositions = radiocarbon_data$depth))
median_predicted_ages <- apply(predicted_ages, 2, median, na.rm = TRUE)

loc1503_Samp$calibrated_age <- median_predicted_ages
write.csv(loc1503_Samp, file = "Locality20.csv", row.names=FALSE)


# Add predicted ages to the data frame
radiocarbon_data$calibrated_age <- ifelse(
  calibrate,
  median_calibrated_vector, # Use calibrated age for known depths
  predicted_ages # Use model-predicted age for estimated depths
)

age = ifelse(all_depths %in% known_depths, ages[match(all_depths, known_depths)], NA)





loc1503_geochron <-  neotoma2::get_sites(datasetid = 8257)

View(loc1503_geochron$samples)
geochron_dl <- neotoma2::get_downloads(loc1503_geochron)




loc513 <- neotoma2::get_sites(siteid = 513)
loc513_geochron <- neotoma2::get_datasets(loc513, all_data = TRUE) %>%
  neotoma2::filter(datasettype == "geochronologic")

loc513_geochron_dl <- loc513_geochron %>%
  neotoma2::get_downloads()

loc513_geochron<-loc513_geochron_dl@sites[["site"]]@collunits@collunits[[1]]@chronologies@chronologies[[1]]@chroncontrols











# Calibrate known dates
calibrated_dates <- BchronCalibrate(
  ages = known_dates$age,
  ageSds = known_dates$error,
  calCurves = known_dates$curve
)

# Fit the age-depth model
age_depth_model <- Bchronology(
  ages = known_dates$age,
  ageSds = known_dates$error,
  calCurves = known_dates$curve,
  positions = known_dates$depth, # Depths of known dates
  positionThicknesses = known_dates$thickness # Thickness of layers
)

predicted_ages <- predict(age_depth_model, newPositions = loc1503_Samp$depth)

calibrated_age <- ifelse(
  known_dates,
  summary(calibrated_dates)$CalBP_Median,
  predicted_ages
)
