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

# Load data
spatialData <- read.csv("IndividualSummaries/PiceaMin.csv")

# Find presence/absence for spatial data
presence <- data.frame(spatialData[,1:5])
for (i in 6:46) {
  presence[,i] <- ifelse(spatialData[,i] > 1, 1, spatialData[,i])
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
         Distance = value) %>%
  filter(Distance != 0)

# Rank distances between each locality
rankDistances <- distances %>%
  group_by(Site_1) %>%
  mutate(Rank_Distance = rank(Distance))

# t0 = 13 ka
t0 = presence$X13000 # do with previous one
t1 = presence$X12500

# If you have an NA in t1, replace with value from t0
t01 <- ifelse(is.na(t1), t0, t1)

# Are there any differences between corrected t1 and t0?
diff <- (t01 > t0) # make greater than 1?

presentLocalities <- spatialData$sitename[which(t0 == 1)]

distance <- 0
for (i in length(diff)) {
  
  
  if (diff[i] == TRUE) {
    
    # Find site name for given difference
    sitename = spatialData$sitename[i]
    
    # Find all distances between site and previous sites
    sites <- Distances %>%
      filter(Site_1 == sitename,
             Site_2 == presentLocalities)
    
    # Pick smallest distance
    distance <- distance + min(sites$Distance)
  }
  
  # replace initial t with t01!
}

# ASSUMPTIONS:
## Once a site it marked "present," it remains present through NAs until we see an absence

############################################################################################
############################################################################################
############################################################################################