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
monteCarlo <- function(entry, decline, duration, nit) {
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
success <- monteCarlo(8, 5, 9, 100000)
p <- success/100000
p # p = 0.00144


# Load libraries
library(raster)
library(dplyr)
library(tidyr)
library(stringr)

# Load data
spatialData <- read.csv("IndividualSummaries/PiceaMin.csv")

# Find presence/absence for spatial data
presence <- data.frame(spatialData[,1:5])
for (i in 6:46) {
  presence[,i] <- ifelse(spatialData[,i] > 0, 1, spatialData[,i])
}
colnames(presence) <- colnames(spatialData)

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

# Create custom function to calculate average distance between refugium and colonized localities in time interval of interest
findDistance <- function (start, end, refugium) {
  
  # Initialize empty data frame to store results
  results <- data.frame(time = numeric(0), 
                        average.distance = numeric(0))
  
  # Isolate first time bin
  t0 <- dplyr::select(presence, start)
  
  # Create while loop to calculate average distance for each time bin while loop is before end
  while (names(t0) != end) {
    
    # Create dataframe of which sites are colonized for a given time bin
    sites <- (t0 == 1) %>%
      ifelse(is.na(.), FALSE, .) # If there are any NAs that turn into a P/A, mark it FALSE
    
    # Create loop to calculate total distance from refugium for a given time bin.
    totalDistance <- 0
    for (i in 1:as.numeric(length(sites))) {
      if (sites[i] == TRUE) {
        
        # Find site name for given difference
        sitename = spatialData$sitename[i]
        
        # Find all distances between site and previous sites
        individualDistance <- distances %>%
          filter(Site_1 == refugium) %>%
          filter(Site_2 == sitename)
        
        # Pick smallest distance
        totalDistance <- totalDistance + individualDistance$Distance
      }
    }
    
    # Calculate average distance for a given time bin
    averageDistance <- totalDistance / sum(!is.na(t0))
    
    # Store results in data frame
    results[nrow(results) + 1,] = c(names(t0), averageDistance)
    
    # Advance conditions
    t0 <- presence %>% dplyr::select(which(names(presence) == names(t0)) - 1)
  }
  
  # Correct transcription artefacts so correlation can be performed
  results$time <- -1 * as.numeric(str_remove(results$time, "X"))
  results$average.distance <- as.numeric(results$average.distance)
  
  return(results)
}

# Run function
results <- findDistance(start = "X13000",
                        end = "X8000",
                        refugium = "Kollioksak Lake")

# Perform correlation on results
cor.test(results$time, results$average.distance, method = "spearman")

############################################################################################
############################################################################################
############################################################################################