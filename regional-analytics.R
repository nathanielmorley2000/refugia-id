# import dataset and prepare for graphing
pollen <- read.csv("TempFiles/PresenceThroughTime.csv")
pollen$Time <- pollen$Time/1000


#####################################################################
############################## FIGURES ##############################
#####################################################################

# create function for summary plot
plotPollen <- function(taxonMin, taxonMax, taxon, color) {
  # create plot and plot data availability
  plot(pollen$AvailableData, pollen$Time, 
       ylim = rev(range(pollen$Time)), 
       xlim = c(0,30),
       type = "l", 
       lty = 1,
       lwd = 5,
       xlab = "Presence",
       ylab = "Time (ka)",
       main = taxon)
  
  # plot maximum taxon
  lines(taxonMax, pollen$Time,
        ylim =rev(range(pollen$Time)),
        type = "l",
        lty = 1,
        lwd = 4,
        col = color)
  
  # plot maximum taxon
  lines(taxonMin, pollen$Time,
        ylim =rev(range(pollen$Time)),
        type = "l",
        lty = 2,
        lwd = 4,
        col = color)
  abline(h=11.7, lty = 3, lwd = 5, col="navy")
}

# direct software to save as .svg to local directory
svg('Results/Figure3.svg', 
    width = 24,
    height = 24,
    pointsize = 40)

# set up graphical parameters
par(mfrow = c(2,2),
    mar = c(4.1, 4.4, 4.1, 1.9))

# call individual plots
plotPollen(pollen$SalixMin, pollen$SalixMax, "Salix", "red3")
plotPollen(pollen$PopulusMin, pollen$PopulusMax, "Populus", "deepskyblue2")
plotPollen(pollen$PiceaMin, pollen$PiceaMax, "Picea", "springgreen3")

# create legend
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("center", legend=c("Localities with Data", "Onset of Holocene",
                          "Salix (Maximum)", "Salix (Minimum)", "Populus (Maximum)", "Populus (Minimum)",
                          "Picea (Maximum)", "Picea (Minimum)"),
       col=c("black", "navy", "red3", "red3", "deepskyblue2", "deepskyblue2",  "springgreen3", "springgreen3"), 
       lty=c(1,3,1,2,1,2,1,2),
       lwd=c(5,5,4,4,4,4,4,4),
       cex=0.8,
       pt.cex=1,
       x.intersp = 0.5,
       xpd=T,
       ncol=2,
       text.width = 0.5,
       title = "Legend")

# stop saving as .svg file
dev.off()

# create function individual plots for Figs. 4, S1-S2
indiv.plot <- function(taxon, color){
  plot(pollen$AvailableData, pollen$Time, 
       ylim = rev(range(pollen$Time)),
       xlim = c(0,30),
       type = "l", 
       lty = 1,
       lwd = 5,
       xlab = "Presence",
       ylab = "Time (ka)",
       pch = 14)
  lines(taxon, pollen$Time,
        ylim =rev(range(pollen$Time)),
        type = "l",
        lty = 1,
        lwd = 5,
        col = color)
}

# call functions
# Picea
svg('Results/PiceaIndividual.svg', 
    width = 10,
    height = 12,
    pointsize = 30)
indiv.plot(pollen$PiceaMin, "springgreen3")
dev.off()

# Populus
svg('Results/PopulusIndividual.svg', 
    width = 10,
    height = 12,
    pointsize = 30)
indiv.plot(pollen$PopulusMin, "deepskyblue2")
dev.off()

#####################################################################
#####################################################################
#####################################################################



#########################################################################
############################## CORRELATION ##############################
#########################################################################

# Salix maximum and minimum
cor.Salix <- cor.test(pollen$SalixMax, pollen$SalixMin, 
                method = "spearman")
cor.Salix # rho = 0.99, p < 2.2e-16

# Populus maximum and minimum
cor.Populus <- cor.test(pollen$PopulusMax, pollen$PopulusMin, 
                method = "spearman")
cor.Populus # rho = 0.89, p = 4.239e-15

# Picea maximum and minimum
cor.Picea <- cor.test(pollen$PiceaMax, pollen$PiceaMin, 
                       method = "spearman")
cor.Picea # rho = 0.97, p < 2.2e-16

#########################################################################
#########################################################################
#########################################################################



############################################################################################
############################## SIGNIFICANCE OF Picea PATTERN ##############################
############################################################################################

# create function for Monte Carlo simultation (with replacement)
temporalMonteCarlo <- function(entry, decline, duration, nit) {
  # set counters for simulation
  it = 1
  score = 0
  count = 0
  
  # enter loop
  for(it in 1:nit){
    
    # randomly resample incidence
    vec = floor(runif(nrow(pollen), min = 0, max = pollen$AvailableData + 1))
    indices_greater_than_entry = which(vec >= entry)
    
    # check if conditions are satisfied based on inputs
    if(length(indices_greater_than_entry) > 1){
      for(i in 1:length(indices_greater_than_entry)){
        count = 0
        if(i + 1 <= length(indices_greater_than_entry)){
          index_lower = indices_greater_than_entry[i]
          index_higher = indices_greater_than_entry[i + 1]
          for(j in index_lower:index_higher){
            if(vec[j] <= decline){
              count <- count + 1
            }
          }
          if(count >= duration) {
            
            # if conditions are satisfied, count as "success." Otherwise move on
            score = score + 1
            break
          }
        }
      }
    }
  }
  return(score)
}


# significance of going from >= 8 down to <= 5 for 9 time bins, then back up to >= 8
success <- temporalMonteCarlo(9, 6, 10, 100000)
p <- success/100000
p # p = 0.00148


# Load libraries
library(raster)
library(dplyr)
library(tidyr)
library(stringr)

# Load data
spatialData <- read.csv("IndividualSummaries/PiceaMin.csv")

# Create custom function for creating a corrected presence dataframe
correctPresence <- function (data) {
  
  # Find presence/absence for spatial data
  presence <- data.frame(spatialData[,1:5])
  for (i in 6:46) {
    presence[,i] <- ifelse(spatialData[,i] > 0, 1, spatialData[,i])
  }
  colnames(presence) <- colnames(spatialData)
  
  # Create NA-adjusted P/A data
  corrected.Presence <- presence
  for (i in 45:6) {
    j = i + 1
    corrected.Presence[,i] <- ifelse(is.na(corrected.Presence[,i]), corrected.Presence[,j], corrected.Presence[,i])
  }
  colnames(corrected.Presence) <- colnames(spatialData) 
  
  return(corrected.Presence)
}

# Call function
corrected.Presence <- correctPresence(data = spatialData)

# Isolate coordinates for distances
coordinates <- as.matrix(cbind(spatialData$long, spatialData$lat))

# Find distances between each locality
distances <- pointDistance(p1 = coordinates, 
                           p2 = coordinates,
                           lonlat = TRUE, allpairs = TRUE) %>%
  as.data.frame() %>%
  rename_with(~ spatialData$sitename) %>%
  mutate(spatialData$sitename, .before = "Angal Lake") %>%
  pivot_longer(cols = 2:26) %>%
  rename(Site_1 = "spatialData$sitename",
         Site_2 = name,
         Distance = value)


# Create custom function to compare 
findRankDistance <- function (start, stop, refugium) {
  
  # Initialize data frame to store correlation data
  distance_travelled <- data.frame(Age = numeric(0),
                                   Distance = numeric(0))
  
  # Isolate initial and expanded time bins
  t0 <- dplyr::select(corrected.Presence, start)
  t1 <- corrected.Presence %>% dplyr::select(which(names(corrected.Presence) == start) - 1)
  t.final <- dplyr::select(corrected.Presence, stop)
  
  # Isolate distances to only include those with the refugium and rank them
  rankDistances <- distances %>%
    filter(Site_1 == refugium) %>%
    mutate(Rank_Distance = rank(Distance) - 1)
  
  while (names(t1) != names(t.final)) {
    
    # Are there any differences between t1 and t0?
    sites <- (t1 > t0) %>%
      ifelse(is.na(.), FALSE, .) # If there are any NAs that turn into a P/A, mark it FALSE
    
    # Create loop to calculate total distance from refugium for a given time bin.
    totalRank <- 0
    for (i in 1:as.numeric(length(sites))) {
      if (sites[i] == TRUE) {
        
        # Find site name for given difference
        sitename <- spatialData$sitename[i]
        
        # Find all distances between site and previous sites
        individualDistance <- rankDistances %>%
          filter(Site_2 == sitename)
        
        # Pick smallest distance
        totalRank <- totalRank + individualDistance$Rank_Distance
      }
    }
    
    # Calculate average ranked distance for a given time bin
    averageRank <- totalRank / length(which(sites == TRUE)) #length(which(!is.na(t1)))
    
    # Store results in data frame
    distance_travelled[nrow(distance_travelled) + 1,] = c(names(t1), averageRank)
    
    # Advance conditions
    t0 <- t1
    t1 <- corrected.Presence %>% dplyr::select(which(names(corrected.Presence) == names(t0)) - 1)
  }
  
  # Correct transcription artefacts so correlation can be performed
  distance_travelled$Age <- -1 * as.numeric(str_remove(distance_travelled$Age, "X"))
  distance_travelled$Distance <- as.numeric(distance_travelled$Distance)
  
  return(distance_travelled)
}
distance_travelled <- findRankDistance(start = "X12000", stop = "X6500", refugium = "Kollioksak Lake")


# Calculate correlation for observed data
result <- cor.test(distance_travelled$Age, distance_travelled$Distance, method = "spearman")

# Create custom Monte Carlo to identify spatial pattern
spatialMonteCarlo <- function(distance_travelled, result, nit) {
  
  # Reset counter to 0
  it <- 0
  for (i in 1:nit) {
    test <- runif(n = nrow(distance_travelled), 1, 24)
    
    testResult <- cor.test(distance_travelled$Age, test, method = "spearman")
    
    # Check that test correlation is significant
    if (testResult$p.value < 0.05) {
      
      # Check whether test correlation is better than actual correlation
      if (testResult$estimate > result$estimate) {
        
        # Add to counter
        it <- it + 1
      }
    }
  }
  
  # Calculate success ratio
  success <- it / nit
  success
}
spatialMonteCarlo(distance_travelled = distance_travelled, result = result, nit = 10000)

############################################################################################
############################################################################################
############################################################################################