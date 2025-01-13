# load libraries
library("dplyr")
library("tidyr")
library("ggplot2")
library("gganimate")
library("gifski")
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("viridis")


# load dataset
mapdata <- read.csv("IndividualSummaries/PiceaMin.csv", check.names = FALSE)

# Load country boundaries
countries <- ne_countries(scale = "medium", returnclass = "sf")

# Validate required columns
required_cols <- c("sitename", "lat", "long")

# remove siteid and datasetid columns, if present
unwanted_columns = c("siteid", "datasetid")
existing_columns = colnames(mapdata)
columns_to_remove = intersect(existing_columns, unwanted_columns)
mapdata <- mapdata %>% 
  dplyr::select(-all_of(columns_to_remove))


# Identify time columns (all except the required ones)
time_cols <- setdiff(names(mapdata), required_cols)
if(length(time_cols) == 0) {
  showNotification("No numerical columns found for time categories. Please include at least one time-based numerical column.", type = "error")
  return(NULL)
}

# transform data into a method that can 
mapdata <- mapdata %>%
  pivot_longer(cols = all_of(time_cols), 
               names_to = "time", values_to = 
                 "value") %>%
  drop_na(value) %>%
  mutate(time = as.numeric(as.character(time))) %>%
  mutate(time = factor(time, levels = sort(unique(time)))) %>%
  mutate(value = as.numeric(value))

find_bbox <- function(map_data) {
  # find original bounding box
  coordinates = st_as_sf(map_data, coords = c("long", "lat"), crs = 4326)
  bbox = st_bbox(coordinates)
  
  # set margin to 0.1 and find original height and width
  margin = 0.1
  width = bbox$xmax - bbox$xmin
  height = bbox$ymax - bbox$ymin
  
  # adjust coordinates to accommodate margin
  xmin = bbox$xmin - width * margin
  xmax = bbox$xmax + width * margin
  ymin = bbox$ymin - height * margin
  ymax = bbox$ymax + height * margin
  
  # create new bounding box and convert so it can be recognized by ggplot2
  bbox_coords = matrix(c(xmin, ymin,  # lower-left
                         xmax, ymin,  # lower-right
                         xmax, ymax,  # upper-right
                         xmin, ymax,  # upper-left
                         xmin, ymin), # close the polygon
                       ncol = 2, byrow = TRUE)
  bbox_polygon = st_polygon(list(bbox_coords))
  bbox_sf = st_sfc(bbox_polygon, crs = 4326)
  expanded_bbox <- st_bbox(bbox_sf)
  
  # finish function and return adjusted bounding box values  
  return(expanded_bbox)}


# custom function in functions.R that automatically finds bounding box from coordinates on datasheet
expanded_bbox <- find_bbox(mapdata)

# further modify the bounding box so it can be placed on map
expanded_bbox_sfc <- st_as_sfc(expanded_bbox)


map_data <- mapdata %>%
  dplyr::mutate(time = as.numeric(as.character(time)))

# create base map
base_map <- ggplot() +
  geom_sf(data = countries, fill = "lightgrey", color = "black") +  # Background countries
  geom_sf(data = expanded_bbox_sfc, fill = NA, color = "black", lwd = 2) +  # Bounding box
  coord_sf(xlim = c(expanded_bbox$xmin, expanded_bbox$xmax),
           ylim = c(expanded_bbox$ymin, expanded_bbox$ymax),
           expand = FALSE) +  # Zoom into bounding box
  xlab("Longitude") +
  ylab("Latitude")+
  theme_minimal(base_size = 20)

# plot data onto base map with black points being n=0 and coloured points reflecting number of grains >1
map_with_data <- base_map +
  geom_point(data = map_data, aes(x = long, y = lat, color = value, group = time), size = 10) +
  scale_size_continuous(guide = 'none') +
  
  # Set up a dual color scale: 0 values as white, others with viridis gradient
  scale_color_gradientn(
    colors = c("white", viridis(256)),  # White for 0, viridis for others
    values = scales::rescale(c(0, 1)),  # Ensure 0 is mapped to white
    limits = c(0, max(map_data$value, na.rm = TRUE)),  # Set limits starting from 0
    na.value = NA  # Ensure NA values are not plotted
  ) +
  
  labs(color = "Abundance")

# animate the map through time
map_with_animation <- map_with_data +
  transition_time(-time) +
  ggtitle('Year: {frame_time}',
          subtitle = 'Frame {frame} of {nframes}')
num_years <- length(time_cols)

# save animation to Results folder
anim_save("Results/PiceaMin.gif", animation = animate(map_with_animation, 
                                    nframes = num_years, 
                                    fps = 1.5, 
                                    width = 1600, 
                                    height = 1200,
                                    res = 150))

