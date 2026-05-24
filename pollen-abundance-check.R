# load libraries
library("neotoma2")
library("tidyr")
library("dplyr")
library("readr")
library("ggplot2")
library("ggtext")
library("ggpubr")

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
        dplyr::filter(elementtype == "pollen") %>%
        group_by(age) %>%
        mutate(value = (value / sum(value))) %>%
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

loc_pollen <- neotomaPollen(site_ids = c(2232, 1503, 1355),
                            taxa = c("Populus", "Salix", "Alnus", "Betula", "Picea"))

# pivot table to display each taxon as a separate column
pollen_wide <- loc_pollen %>%
  pivot_wider(
    names_from = variablename, 
    values_from = value         
  )

# import previously calibrated data
time_calibrations <- read.csv("TempFiles/calibratedDates.csv") %>%
  filter(siteid %in% c(2232, 1503, 1355)) %>%
  select(calibrated_age)

# append to test taxon
calibrated_test <- cbind(pollen_wide, time_calibrations)

# remove datapoints greater than 20 ka
calibrated_20max <- calibrated_test %>%
  filter(calibrated_age <= 20000) %>%
  mutate(sitename = recode(sitename, 
                           "Ruppert Lake" = "Locality 17",
                           "Kollioksak Lake" = "Locality 19",
                           "Joe Lake" = "Locality 20"))
    
# collect individual spectra between Kollioksak Lake and adjacent lakes with stronger geochronologies 
alnus <- ggplot(data = calibrated_20max, aes(x = calibrated_age, y = Alnus, colour = sitename, shape = sitename)) +
  geom_point() +
  ylab("Percent of sample") +
  xlab("Age (kcal yr BP)") + 
  ggtitle("*Alnus*") +
  labs(colour = "Locality", shape = "Locality") +
  theme_test() +
  theme(plot.title = ggtext::element_markdown())

betula <- ggplot(data = calibrated_20max, aes(x = calibrated_age, y = Betula, colour = sitename, shape = sitename)) +
  geom_point() +
  ylab("Percent of sample") +
  xlab("Age (kcal yr BP)") + 
  ggtitle("*Betula*") +
  labs(colour = "Locality", shape = "Locality") +
  theme_test() +
  theme(plot.title = ggtext::element_markdown())

picea <- ggplot(data = calibrated_20max, aes(x = calibrated_age, y = Picea, colour = sitename, shape = sitename)) +
  geom_point() +
  ylab("Percent of sample") +
  xlab("Age (kcal yr BP)") + 
  ggtitle("*Picea*") +
  labs(colour = "Locality", shape = "Locality") +
  theme_test() +
  theme(plot.title = ggtext::element_markdown())

populus <- ggplot(data = calibrated_20max, aes(x = calibrated_age, y = Populus, colour = sitename, shape = sitename)) +
  geom_point() +
  ylab("Percent of sample") +
  xlab("Age (kcal yr BP)") + 
  ggtitle("*Populus*") +
  labs(colour = "Locality", shape = "Locality") +
  theme_test() +
  theme(plot.title = ggtext::element_markdown())

salix <- ggplot(data = calibrated_20max, aes(x = calibrated_age, y = Salix, colour = sitename, shape = sitename)) +
  geom_point() +
  ylab("Percent of sample") +
  xlab("Age (kcal yr BP)") + 
  ggtitle("*Salix*") +
  labs(colour = "Locality", shape = "Locality") +
  theme_test() +
  theme(plot.title = ggtext::element_markdown())


# compile spectra and export as SVG
composite <- ggarrange(alnus, betula, populus, salix, picea,
                       ncol = 2, nrow = 3)

svg("Results/pollenSpikes.svg", width = 10, height = 7)
composite
dev.off()
