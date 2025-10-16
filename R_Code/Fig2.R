################################################################################
# This code is written by Dazhi Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: yangdazhi.nus@gmail.com
################################################################################

#Clear all workspace
rm(list = ls(all = TRUE))
libs <- c("dplyr", "raster", "ggplot2", "RColorBrewer")
invisible(lapply(libs, library, character.only = TRUE))

################################################################################
# global input
################################################################################
# main working directory
dir0 <- "/Users/seryangd/Library/CloudStorage/Dropbox/Working papers/Angles"
# source the functions
source(file.path(dir0, "Code/functions.R"))
# some plotting parameters
plot.size = 9; line.size = 0.15; point.size = 2; legend.size = 0.4; text.size = plot.size*5/14;
################################################################################

# get nominal pixel coordinates
loc <- get_lon_lat(res = "4000M", directory = file.path(dir0, "Data"))
lat <- loc$lat
lon <- loc$lon

# get auxiliary inputs for solar positioning, which include altitude, surface pressure, air temperature, and delta T
alt <- get_dem(lat, lon, directory = file.path(dir0, "Data"))
sp <- get_sp(lat, lon, directory = file.path(dir0, "Data"))
t2m <- get_t2m(lat, lon, directory = file.path(dir0, "Data"))

# perform the most rigorous positioning
sat_angles <- satpos(lat = lat, lon = lon, alt = alt, 
                     ar_cor = TRUE, sp = sp/100, t2m = t2m)
theta <- sat_angles$SatZen
phi <- sat_angles$SatAzi

# performing differet versions of satellite positioning over the full disk area
sat_angles1 <- satpos(lat = lat, lon = lon)
sat_angles2 <- satpos(lat = lat, lon = lon, ar_cor = TRUE, sp = sp/100, t2m = t2m)

# to show that altitude has no effect on azimuth
phi1 <- sat_angles1$SatAzi
plot(raster(phi-phi1))

# prepare data
data.plot <- NULL
for(i in 1:2)
{
  # get the other versions of zenith and azimuth
  VACx <- get(paste0("sat_angles", i))
  theta_VACx <- VACx$SatZen
  
  # compute the differences
  Delta_theta <- theta - theta_VACx
  
  # store them for plotting
  ras_Delta_theta <- raster(Delta_theta, xmn = 1, xmx = 2748, ymn = 1, ymx = 2748)
  ras_Delta_theta <- aggregate(ras_Delta_theta, fact = 4)
  ras_Delta_theta <- rasterToPoints(ras_Delta_theta)
  ras_Delta_theta <- as_tibble(ras_Delta_theta) %>%
    mutate(y = rev(y)) 
  
  # arrange in tibble for plotting
  data.plot <- data.plot %>%
    bind_rows(., tibble(x = ras_Delta_theta$x, y = ras_Delta_theta$y, 
                        diff = ras_Delta_theta$layer, version = i, 
                        angle = "Satellite zenith"))
}

################################################################################
# plot zenith angle
################################################################################
min_index <- which.min(theta) 
position <- arrayInd(min_index, dim(theta))
# make a point to plot the position of the sun
data.point <- tibble(x = position[,2], y = position[,1])

ras_theta <- raster(theta)
extent(ras_theta) <- c(1, 2748, 1, 2748)  
ras_theta <- aggregate(ras_theta, fact=4)
ras_theta <- rasterToPoints(ras_theta)
ras_theta <- as_tibble(ras_theta) %>%
  mutate(y = rev(y)) 

ras_phi <- raster(phi)
extent(ras_phi) <- c(1, 2748, 1, 2748)  
#ras_phi <- aggregate(ras_phi, fact=4)
ras_phi <- rasterToPoints(ras_phi)
ras_phi <- as_tibble(ras_phi) %>%
  mutate(y = rev(y)) 

p1 <- ggplot(ras_theta) + 
  geom_raster(aes(x = x, y = y, fill = layer)) +
  geom_point(data = data.point, aes(x = x, y = y), shape = 1, size = point.size)+
  coord_fixed() +
  viridis::scale_fill_viridis(name = expression(paste(italic(theta), " [", degree, "]")), direction = -1, option = "B", limits = c(0,90), breaks = c(0, 30,60,90)) +
  scale_x_continuous("Column index", breaks = c(1, 1000, 2000, 2748)) +
  scale_y_reverse("Row index", breaks = c(1, 1000, 2000, 2748)) + #
  theme_bw() +
  theme(plot.margin = unit(c(0.8,0.4,0,0.4), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(family = "Times", size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size/2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1,0,0.1,0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0,0.1,0,0.1), "lines")), panel.spacing = unit(0 , "lines"), legend.position = "bottom", legend.key.width = unit(1.2, "lines"), legend.key.height = unit(0.4, "lines"), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7,0,0,0.0), "lines"), legend.direction = "horizontal")

pA <- ggpubr::ggarrange(p1, ncol = 1, labels = c("(a)"), font.label = list(size = plot.size, color = "black", face = "plain", family = "Times"))

################################################################################
# plot azimuth angle (not displayed in paper)
################################################################################

p2 <- ggplot(ras_phi) + 
  geom_raster(aes(x = x, y = y, fill = layer)) +
  geom_point(data = data.point, aes(x = x, y = y), shape = 1, size = point.size)+
  coord_fixed() +
  viridis::scale_fill_viridis(name = expression(paste(phi, " [", degree, "]")), direction = -1, option = "D", limits = c(0,360), breaks = c(0, 100, 200, 300)) +
  scale_x_continuous("Column index", breaks = c(1, 1000, 2000, 2748)) +
  scale_y_reverse("Row index", breaks = c(1, 1000, 2000, 2748)) + #
  theme_bw() +
  theme(plot.margin = unit(c(0.8,0.4,0,0.4), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(family = "Times", size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size/2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1,0,0.1,0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0,0.1,0,0.1), "lines")), panel.spacing = unit(0 , "lines"), legend.position = "bottom", legend.key.width = unit(1.2, "lines"), legend.key.height = unit(0.4, "lines"), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7,0,0,0.0), "lines"), legend.direction = "horizontal")

################################################################################
# plot zenith angle, using other inputs, i.e., VAC1 and VAC2
################################################################################

tmp <- data.plot %>% 
  filter(version == 1, angle == "Satellite zenith") 
# Calculate quantile breaks
no_classes <- 10
quantiles <- quantile(tmp$diff, probs = seq(0, 1, length.out = no_classes + 1), na.rm = TRUE)

p3 <- ggplot(data.plot %>% 
               filter(version %in% c(1:2), angle == "Satellite zenith") %>% 
               mutate(diff > quantiles[length(quantiles)], quantiles[length(quantiles)], diff) %>%
               mutate(diff < quantiles[1], quantiles[1], diff)) + 
  geom_raster(aes(x = x, y = y, fill = diff)) +
  geom_point(data = data.point, aes(x = x, y = y), shape = 1, size = point.size)+
  facet_wrap(~version, labeller = labeller(version = c("1" = "VAC 1", "2" = "VAC 2"))) +
  coord_fixed() +
  #scale_fill_gradientn(name = expression(paste("(W ", m^-2, ")")), colors = colorRampPalette(rev(brewer.pal(7, "PuBu")))(7), values = scales::rescale(quantiles, to = c(0,1))) +
  #scale_fill_gradient2(name = expression(paste(Delta, theta, " [", degree, "]")), low = brewer.pal(11, "RdBu")[11:8], mid = "white", high = brewer.pal(11, "RdBu")[4:1], midpoint = 0, transform = ggallin::ssqrt_trans, breaks = c(-0.5, -0.2, -0.05, 0), limits = c(-0.5, 0)) + 
  scale_fill_gradient2(name = expression(paste(Delta, theta, " [", degree, "]")), low = brewer.pal(11, "RdBu")[11:8], mid = "white", high = brewer.pal(11, "RdBu")[4:1], midpoint = 0, limits = c(-0.3, 0.01), transform = ggallin::ssqrt_trans, breaks = c(-0.2, -0.05, 0)) + 
  scale_x_continuous("Column index", breaks = c(1, 1000, 2000, 2748)) +
  scale_y_reverse("Row index", breaks = c(1, 1000, 2000, 2748)) + #
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0,0), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(family = "Times", size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size/2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1,0,0.1,0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0,0.1,0,0.1), "lines")), panel.spacing = unit(0.7, "lines"), legend.position = "bottom", legend.key.width = unit(1.2, "lines"), legend.key.height = unit(0.4, "lines"), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7,0,0,0), "lines"), legend.direction = "horizontal", legend.title = element_text(margin = margin(r = 15)))

pB <- ggpubr::ggarrange(p3, ncol = 1, labels = c("(b)"), font.label = list(size = plot.size, color = "black", face = "plain", family = "Times"))

multiplot <- list(pA, pB)

lay <- rbind(c(1, 2, 2))

p <- gridExtra::grid.arrange(grobs = multiplot, layout_matrix = lay)

setwd(file.path(dir0, "tex"))
ggsave(filename = "zenithSat.pdf", plot = p, width = 160, height = 60, unit = "mm")

