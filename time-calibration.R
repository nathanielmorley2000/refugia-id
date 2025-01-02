library("neotoma2")
library("Bchron")
library("dplyr")

loc1503 <- neotoma2::get_sites(siteid = 1503)
loc1503_pollen <- neotoma2::get_datasets(loc1503, all_data = TRUE) %>%
  neotoma2::filter(datasettype == "pollen")
summary(loc1503)

loc1503_dl <- loc1503_pollen %>%
  get_downloads()

loc1503_Samp <- samples(loc1503_dl) %>%
  group_by(sitename, lat, long, siteid, datasetid, depth, age, variablename) %>%
  summarize(value = sum(value), .groups = "keep") %>%
  group_by(sitename, lat, long, siteid, datasetid, depth, age) %>%
  #do(ensure_taxon_present(., taxon)) %>%
  ungroup() %>%
  dplyr::filter(variablename == "Picea") %>%
  dplyr::select(sitename, lat, long, siteid, datasetid, value, depth, age)


all_depths <- as.numeric(loc1503_Samp$depth)
known_depths <- c(30, 280, 405, 445)
ages = c(2540, 6620, 9350, 14140)
errors = c(60, 90, 90, 440)

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

age = ifelse(all_depths %in% known_depths, ages[match(all_depths, known_depths)], NA),



























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
