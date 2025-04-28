# load libraries
library("neotoma2")
library("Bchron")
library("tidyr")
library("dplyr")

# create directories to store outputs
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
                                         2270, 2245, 1977, 10102, 1955, 207, 2232, 1699, 1503,
                                         1355, 2551, 13690, 11575, 11579, 11583))

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
                                         2270, 2245, 1977, 10102, 1955, 207, 2232, 1699, 1503,
                                         1355, 2551, 13690, 11575, 11579, 11583),
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

# save dataframe as a .csv file for easy recall
write.csv(output_df, "TempFiles/calibratedDates.csv", row.names = FALSE)

#############################################################################
#############################################################################
#############################################################################



###################################################################################
############################## BIN AND ORGANIZE DATA ##############################
###################################################################################

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
  
  # selects sample with smallest value (for specified taxon) in 500-year time bin
  if (samplingProtocol == "Minimum") {
  data_filtered = timeCorrected %>%
    group_by(sitename, Year_Bin) %>%
    slice_min(order_by = value, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::filter(Year_Bin >= yearMin) %>%
    dplyr::filter(Year_Bin <= yearMax)
  
  # selects sample with largest value (for specified taxon) in 500-year time bin
  } else if (samplingProtocol == "Maximum") {
    data_filtered = timeCorrected %>%
      group_by(sitename, Year_Bin) %>%
      slice_max(order_by = value, with_ties = FALSE) %>%
      ungroup() %>%
      dplyr::filter(Year_Bin >= yearMin) %>%
      dplyr::filter(Year_Bin <= yearMax)
  }
  
  # creates pivot table with correctly ordered time bins
  ordered_years = sort(unique(data_filtered$Year_Bin))
  pivot_table = data_filtered %>%
    select(sitename, siteid, datasetid, lat, long, Year_Bin, value) %>%
    pivot_wider(names_from = Year_Bin, values_from = value, values_fill = list(Taxon_Abundance = NA)) %>%
    select(sitename, siteid, datasetid, lat, long, all_of(as.character(ordered_years)))
  
  # return pivot table
  return(pivot_table)
}

# call function and export as .csv for Salix using minimum protocol
SalixMin <- organizeData(timeBin = 500,
                         taxon = "Salix",
                         samplingProtocol = "Minimum",
                         yearMin = 0,
                         yearMax = 20000)
write.csv(SalixMin, "IndividualSummaries/SalixMin.csv", row.names = FALSE)

# call function and export as .csv for Salix using maximum protocol
SalixMax <- organizeData(timeBin = 500,
                         taxon = "Salix",
                         samplingProtocol = "Maximum",
                         yearMin = 0,
                         yearMax = 20000)
write.csv(SalixMax, "IndividualSummaries/SalixMax.csv", row.names = FALSE)

# call function and export as .csv for Populus using minimum protocol
PopulusMin <- organizeData(timeBin = 500,
                         taxon = "Populus",
                         samplingProtocol = "Minimum",
                         yearMin = 0,
                         yearMax = 20000)
write.csv(PopulusMin, "IndividualSummaries/PopulusMin.csv", row.names = FALSE)

# call function and export as .csv for Populus using maximum protocol
PopulusMax <- organizeData(timeBin = 500,
                         taxon = "Populus",
                         samplingProtocol = "Maximum",
                         yearMin = 0,
                         yearMax = 20000)
write.csv(PopulusMax, "IndividualSummaries/PopulusMax.csv", row.names = FALSE)

# call function and export as .csv for Picea using minimum protocol
PiceaMin <- organizeData(timeBin = 500,
                         taxon = "Picea",
                         samplingProtocol = "Minimum",
                         yearMin = 0,
                         yearMax = 20000)
write.csv(PiceaMin, "IndividualSummaries/PiceaMin.csv", row.names = FALSE)

# call function and export as .csv for Picea using maximum protocol
PiceaMax <- organizeData(timeBin = 500,
                         taxon = "Picea",
                         samplingProtocol = "Maximum",
                         yearMin = 0,
                         yearMax = 20000)
write.csv(PiceaMax, "IndividualSummaries/PiceaMax.csv", row.names = FALSE)

# create function to produce regional summary for all taxa for presence through time
regionalSummary <- function(dataframes) {
  
  # initialize empty dataframe to store output
  all_summaries <- NULL
  
  # loop through each dataset
  for(df_name in dataframes){
    
    # Dynamically get the dataframe from the global environment using `get()`
    df <- get(df_name)
    
    # convert to long format for analysis
    long_df <- df %>%
      pivot_longer(
        cols = "0":last_col(),
        names_to = ("time"),
        values_to = ("abundance")
      )
    
    # summarize number of sites with values >0 per time bin
    summary = long_df %>%
      group_by(time) %>%
      summarize(
        localities_with_data = sum(!is.na(abundance)),
        localities_with_pollen = sum(abundance > 0, na.rm = TRUE)
      ) %>%
      mutate(
        time = as.numeric(time),
        localities_with_data = as.numeric(localities_with_data),
        localities_with_pollen = as.numeric(localities_with_pollen)
      ) %>%
      arrange(time)
    
    # If it's the first dataframe, initialize 'all_summaries'
    if (is.null(all_summaries)) {
      all_summaries <- data.frame(summary$localities_with_pollen)
      colnames(all_summaries) <- df_name  # Name the column with the dataframe name
    } else {
      # Add the new column with the dataframe's name
      all_summaries[[df_name]] <- summary$localities_with_pollen
    }
    
  }
  
  # add columns for time and data availability to front of output dataframe
  time = summary$time
  AvailableData = summary$localities_with_data
  all_summaries = cbind(time, AvailableData, all_summaries)
  colnames(all_summaries)[1:2] = c("Time", "AvailableData")
  
  return(all_summaries)
}

# call function and save as .csv file in TempFiles
presencethroughtime <- regionalSummary(dataframes = c("SalixMin", 
                                                      "SalixMax", 
                                                      "PopulusMin", 
                                                      "PopulusMax", 
                                                      "PiceaMin", 
                                                      "PiceaMax"))
write.csv(presencethroughtime, "TempFiles/PresenceThroughTime.csv", row.names = FALSE)

###################################################################################
###################################################################################
###################################################################################