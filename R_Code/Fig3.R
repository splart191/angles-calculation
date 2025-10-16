################################################################################
# This code is written by Dazhi Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: yangdazhi.nus@gmail.com
################################################################################

#Clear all workspace
rm(list = ls(all = TRUE))
libs <- c("dplyr", "lubridate", "raster", "ggplot2", "RColorBrewer")
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

datetime <- lubridate::ymd_hm("2024-09-23 04:00")

# get nominal pixel coordinates
loc <- get_lon_lat(res = "4000M", directory = file.path(dir0, "Data"))
lat <- loc$lat
lon <- loc$lon

# get AGRI scanning time from a HDF file
row_scan_time <- get_row_scan_time(datetime, directory = file.path(dir0, "Data"))
scan_time <- allocate_scan_time(row_scan_time$time_start, row_scan_time$time_end, lat = lat, lon = lon)
scan_time <- as.POSIXct(scan_time, origin = "1970-01-01", tz = "UTC")

# get auxiliary inputs for solar positioning, which include altitude, surface pressure, air temperature, and delta T
alt <- get_dem(lat, lon, directory = file.path(dir0, "Data"))
sp <- get_sp(lat, lon, directory = file.path(dir0, "Data"))
t2m <- get_t2m(lat, lon, directory = file.path(dir0, "Data"))
deltaT <- get_deltaT(datetime, directory = file.path(dir0, "Data"))

# perform the most rigorous solar positioning
sun_angles <- solpos(lat = lat, lon = lon, alt = alt, 
                     sp = sp/100, t2m = t2m, deltaT = deltaT, 
                     datetime = scan_time, tz = 0)
theta0 <- sun_angles$SolZen
phi0 <- sun_angles$SolAzi
min_index <- which.min(theta0) 
position <- arrayInd(min_index, dim(theta0))
data.point1 <- tibble(x = position[,2], y = position[,1])

# perform the most rigorous satellite positioning
sat_angles <- satpos(lat = lat, lon = lon, alt = alt, 
                     ar_cor = TRUE, sp = sp/100, t2m = t2m)
theta <- sat_angles$SatZen
phi <- sat_angles$SatAzi
min_index <- which.min(theta) 
position <- arrayInd(min_index, dim(theta))
data.point2 <- tibble(x = position[,2], y = position[,1])

# compute the sun glint angles
sun_glint <- sunglint(SolZen = theta0, SolAzi = phi0,
                      SatZen = theta, SatAzi = phi)
plot(raster(sun_glint))

################################################################################
# quantify the difference between our calculation with L1_GEO
################################################################################

# read the Fengyun-4B official angles
sun_angles_FY <- hdf_read(datetime = datetime, name = "GEO")
theta0_FY <- as.matrix(sun_angles_FY$NOMSunZenith)
phi0_FY <- as.matrix(sun_angles_FY$NOMSunAzimuth)
theta_FY <- as.matrix(sun_angles_FY$NOMSatelliteZenith)
phi_FY <- as.matrix(sun_angles_FY$NOMSatelliteAzimuth)

# find out whether the NA locations in the L1_GEO files are identical with the lat-lon look-up table
identical(which(is.na(loc$lat)), which(is.na(theta0_FY)))
# and they are not!
# therefore, we need to mask the zen_FY and azi_FY with lat
remove <- which(is.na(loc$lat), arr.ind = TRUE)
theta0_FY[remove] <- NA
phi0_FY[remove] <- NA
theta_FY[remove] <- NA
phi_FY[remove] <- NA

# difference in theta_0
Delta_theta0 <- raster(theta0 - theta0_FY)
extent(Delta_theta0) <- c(1, 2748, 1, 2748) 
Delta_theta0 <- aggregate(Delta_theta0, fact=4)
Delta_theta0 <- rasterToPoints(Delta_theta0)
Delta_theta0 <- as_tibble(Delta_theta0) %>%
  mutate(y = rev(y)) 

# difference in phi_0
Delta_phi0 <- raster(phi0 - phi0_FY)
extent(Delta_phi0) <- c(1, 2748, 1, 2748) 
Delta_phi0 <- aggregate(Delta_phi0, fact=4)
Delta_phi0 <- rasterToPoints(Delta_phi0)
Delta_phi0 <- as_tibble(Delta_phi0) %>%
  mutate(y = rev(y)) 

# difference in theta
Delta_theta <- raster(theta - theta_FY)
extent(Delta_theta) <- c(1, 2748, 1, 2748) 
Delta_theta <- aggregate(Delta_theta, fact=4)
Delta_theta <- rasterToPoints(Delta_theta)
Delta_theta <- as_tibble(Delta_theta) %>%
  mutate(y = rev(y)) 

# difference in phi
Delta_phi <- raster(phi - phi_FY)
extent(Delta_phi) <- c(1, 2748, 1, 2748) 
Delta_phi <- aggregate(Delta_phi, fact=4)
Delta_phi <- rasterToPoints(Delta_phi)
Delta_phi <- as_tibble(Delta_phi) %>%
  mutate(y = rev(y)) 

################################################################################
# plot angle differences
################################################################################

p1 <- ggplot(Delta_theta0) + 
  geom_raster(aes(x = x, y = y, fill = layer)) +
  geom_point(data = data.point1, aes(x = x, y = y), shape = 1, size = point.size)+
  coord_fixed() +
  scale_fill_gradient2(name = expression(paste(Delta, theta[0], " [", degree, "]")), low = brewer.pal(11, "RdBu")[11:10], mid = "white", high = brewer.pal(11, "RdBu")[2:1], midpoint = 0, transform = ggallin::ssqrt_trans, breaks = c(-0.6, -0.2, 0)) +
  scale_x_continuous("Column index", breaks = c(1, 1000, 2000, 2748)) +
  scale_y_reverse("Row index", breaks = c(1, 1000, 2000, 2748)) + #
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(family = "Times", size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size/2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1,0,0.1,0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0,0.1,0,0.1), "lines")), panel.spacing = unit(0 , "lines"), legend.position = "bottom", legend.key.width = unit(1.2, "lines"), legend.key.height = unit(0.4, "lines"), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7,1,0,0), "lines"), legend.direction = "horizontal")

p2 <- ggplot(Delta_phi0 %>% filter(layer < 1)) + 
  geom_raster(aes(x = x, y = y, fill = layer)) +
  geom_point(data = data.point1, aes(x = x, y = y), shape = 1, size = point.size)+
  coord_fixed() +
  scale_fill_gradient2(name = expression(paste(Delta, phi[0], " [", degree, "]")), low = brewer.pal(11, "RdBu")[11:10], mid = "white", high = brewer.pal(11, "RdBu")[2:1], midpoint = 0, transform = ggallin::ssqrt_trans) +
  scale_x_continuous("Column index", breaks = c(1, 1000, 2000, 2748)) +
  scale_y_reverse("Row index", breaks = c(1, 1000, 2000, 2748)) + #
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(family = "Times", size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size/2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1,0,0.1,0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0,0.1,0,0.1), "lines")), panel.spacing = unit(0 , "lines"), legend.position = "bottom", legend.key.width = unit(1.2, "lines"), legend.key.height = unit(0.4, "lines"), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7,1,0,0), "lines"), legend.direction = "horizontal")

p3 <- ggplot(Delta_theta) + 
  geom_raster(aes(x = x, y = y, fill = layer)) +
  geom_point(data = data.point2, aes(x = x, y = y), shape = 1, size = point.size)+
  coord_fixed() +
  scale_fill_gradient2(name = expression(paste(Delta, theta, " [", degree, "]")), low = brewer.pal(11, "RdBu")[11:10], mid = "white", high = brewer.pal(11, "RdBu")[2:1], midpoint = 0, transform = ggallin::ssqrt_trans, breaks = c(-0.6, -0.2, 0)) +
  scale_x_continuous("Column index", breaks = c(1, 1000, 2000, 2748)) +
  scale_y_reverse("Row index", breaks = c(1, 1000, 2000, 2748)) + #
  theme_bw() +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(family = "Times", size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size/2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1,0,0.1,0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0,0.1,0,0.1), "lines")), panel.spacing = unit(0 , "lines"), legend.position = "bottom", legend.key.width = unit(1.2, "lines"), legend.key.height = unit(0.4, "lines"), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7,1,0,0), "lines"), legend.direction = "horizontal")

pA <- ggpubr::ggarrange(p1, p2, p3, nrow = 1, labels = c("(a)", "(b)", "(c)"), widths = c(1, 1, 1), font.label = list(size = plot.size, color = "black", face = "plain", family = "Times"))

# also check phi
ggplot(Delta_phi) + 
  geom_raster(aes(x = x, y = y, fill = layer)) +
  geom_point(data = data.point2, aes(x = x, y = y), shape = 1, size = point.size)+
  coord_fixed() +
  scale_fill_gradient2(name = expression(paste(Delta, phi, " [", degree, "]")), low = brewer.pal(11, "RdBu")[11:10], mid = "white", high = brewer.pal(11, "RdBu")[2:1], midpoint = 0, transform = ggallin::ssqrt_trans, limits  = c(-0.05, 0.05))


################################################################################
# plot angle distributions
################################################################################

p4 <- ggplot(Delta_theta0) + 
  geom_histogram(aes(x = layer, y = ..density..), bins = 100, fill = "grey80", color = "gray30", linewidth = line.size) +
  scale_x_continuous(name = expression(paste(Delta, theta[0], " [", degree, "]")), limits = c(-0.2, 0.02)) +
  scale_y_continuous(name = expression(paste("Prob. density"))) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.5,0.5), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(family = "Times", size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size/2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1,0,0.1,0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0,0.1,0,0.1), "lines")), panel.spacing = unit(0 , "lines"), legend.position = "bottom", legend.key.width = unit(1.5, "lines"), legend.key.height = unit(0.4, "lines"), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7,0,0,0), "lines"), legend.direction = "horizontal")

p5 <- ggplot(Delta_phi0 %>% filter(layer < 1)) + 
  geom_histogram(aes(x = layer, y = ..density..), bins = 100, fill = "grey80", color = "gray30", linewidth = line.size) +
  scale_x_continuous(name = expression(paste(Delta, phi[0], " [", degree, "]")), limits = c(-0.1, 0.1)) +
  scale_y_continuous(name = expression(paste("Prob. density"))) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.5,0.5), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(family = "Times", size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size/2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1,0,0.1,0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0,0.1,0,0.1), "lines")), panel.spacing = unit(0 , "lines"), legend.position = "bottom", legend.key.width = unit(1.5, "lines"), legend.key.height = unit(0.4, "lines"), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7,0,0,0), "lines"), legend.direction = "horizontal")

p6 <- ggplot(Delta_theta) + 
  geom_histogram(aes(x = layer, y = ..density..), bins = 100, fill = "grey80", color = "gray30", linewidth = line.size) +
  scale_x_continuous(name = expression(paste(Delta, theta, " [", degree, "]")), limits = c(-0.2, 0.02)) +
  scale_y_continuous(name = expression(paste("Prob. density"))) +
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0.5,0.5), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(family = "Times", size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size/2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1,0,0.1,0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0,0.1,0,0.1), "lines")), panel.spacing = unit(0 , "lines"), legend.position = "bottom", legend.key.width = unit(1.5, "lines"), legend.key.height = unit(0.4, "lines"), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7,0,0,0), "lines"), legend.direction = "horizontal")


pB <- ggpubr::ggarrange(p4, p5, p6, nrow = 1, labels = c("(d)", "(e)", "(f)"), widths = c(1, 1, 1), font.label = list(size = plot.size, color = "black", face = "plain", family = "Times"))


multiplot <- list(pA, pB)

lay <- rbind(c(1),
             c(1),
             c(1),
             c(2),
             c(2))

p <- gridExtra::grid.arrange(grobs = multiplot, layout_matrix = lay)

setwd(file.path(dir0, "tex"))
ggsave(filename = "GEO.pdf", plot = p, width = 160, height = 100, unit = "mm")

