# load libraries
library("neotoma2")
library("Bchron")
library("dplyr")

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

loc_test <- neotomaGeochron(site_ids = c(10537, 10539, 513, 2271, 10538, 1396, 1748, 790, 992,  1974, 
                                         2270, 1973, 2245, 1977, 10102, 1955, 207, 2232, 1699, 1503,
                                         1355, 2551, 13690, 11575, 11579, 11583, 846))

# manually call geochronologic controls for Site 11583 (different formatting)
neotomaManual <- function(site_id){

  geochron <- neotoma2::get_datasets(siteid = site_id, all_data = TRUE) %>%
    neotoma2::get_downloads()
  geo_controls<-geochron_dl@sites[["site"]]@collunits@collunits[[1]]@chronologies@chronologies[[1]]@chroncontrols
  
  # Store the results
  if (!is.null(geo_controls)) {
    geo_controls <- as.data.frame(geo_controls)
    geo_controls <- cbind(siteid = site_id, geo_controls)  # Add siteid as the first column
    geochronologic_data[[length(geochronologic_data) + 1]] <- geo_controls
  }
  
}

loc11583_test <- neotomaManual(site_id = 11583)




loc1503_geochron <-  neotoma2::get_sites(sitetid = 9701)

View(loc1503_geochron$samples)
geochron_dl <- neotoma2::get_downloads(loc1503_geochron)


loc11583 <- neotoma2::get_sites(siteid = 11583) %>% 
  get_downloads()


loc11583_geochron <- neotoma2::get_datasets(loc11583, all_data = TRUE) %>%
  neotoma2::filter(datasettype == "geochronologic")
geochron_dl <- loc11583_geochron %>%
  get_downloads()
geochron<-geochron_dl@sites[["site"]]@collunits@collunits[[1]]@chronologies@chronologies[[1]]@chroncontrols




loc1503_pollen <- neotoma2::get_datasets(loc1503, all_data = TRUE) %>%
  neotoma2::filter(datasettype == "pollen")


pollen_dl <- loc1503_pollen %>%
  get_downloads()



pollen_Samp <- samples(pollen_dl) %>%
  group_by(sitename, lat, long, siteid, datasetid, depth, age, variablename) %>%
  summarize(value = sum(value), .groups = "keep") %>%
  group_by(sitename, lat, long, siteid, datasetid, depth, age) %>%
  #do(ensure_taxon_present(., taxon)) %>%
  ungroup() %>%
  dplyr::filter(variablename == "Picea") %>%
  dplyr::select(sitename, lat, long, siteid, datasetid, value, depth, age)


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

# Fit the age-depth model
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
