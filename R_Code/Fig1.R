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
radians <- function (degree) {degree * (pi/180)}
degrees <- function (radian) {radian * (180/pi)}
# some plotting parameters
plot.size = 9; line.size = 0.15; point.size = 2; legend.size = 0.4; text.size = plot.size*5/14;
################################################################################

datetime <- lubridate::ymd_hm("2024-09-23 04:00")

# get nominal pixel coordinates
loc <- get_lon_lat(res = "4000M", directory = file.path(dir0, "Data"))
lat <- loc$lat
lon <- loc$lon
plot(raster(lat))
# identical(which(is.na(loc$lon)), which(is.na(loc$lat)))

# get AGRI scanning time from a HDF file
row_scan_time <- get_row_scan_time(datetime, directory = file.path(dir0, "Data"))
scan_time <- allocate_scan_time(row_scan_time$time_start, row_scan_time$time_end, lat = lat, lon = lon)
scan_time <- as.POSIXct(scan_time, origin = "1970-01-01", tz = "UTC")

# get auxiliary inputs for solar positioning, which include altitude, surface pressure, air temperature, and delta T
alt <- get_dem(lat, lon, directory = file.path(dir0, "Data"))
sp <- get_sp(lat, lon, directory = file.path(dir0, "Data"))
t2m <- get_t2m(lat, lon, directory = file.path(dir0, "Data"))
deltaT <- get_deltaT(datetime, directory = file.path(dir0, "Data"))

# perform the most rigorous positioning
sun_angles <- solpos(lat = lat, lon = lon, alt = alt, 
                     sp = sp/100, t2m = t2m, deltaT = deltaT, 
                     datetime = scan_time, tz = 0)
theta0 <- sun_angles$SolZen
phi0 <- sun_angles$SolAzi

# remove pixels with zenith over 90, which is unphysical
remove <- which(theta0 > 90, arr.ind = TRUE)
theta0[remove] <- 90

# performing differet versions of solar positioning over the full disk area using SPA-R
sun_angles1 <- solpos(lat = lat, lon = lon, datetime = datetime, tz = 0)
sun_angles2 <- solpos(lat = lat, lon = lon, datetime = scan_time, tz = 0)
sun_angles3 <- solpos(lat = lat, lon = lon, alt = alt, datetime = scan_time, tz = 0)
sun_angles4 <- solpos(lat = lat, lon = lon, alt = alt, sp = sp/100, datetime = scan_time, tz = 0)
sun_angles5 <- solpos(lat = lat, lon = lon, alt = alt, sp = sp/100, t2m = t2m, datetime = scan_time, tz = 0)

# prepare data
data.plot <- NULL
for(i in 1:5)
{
  # get the other versions of solar zenith and solar azimuth
  SPAx <- get(paste0("sun_angles", i))
  theta0_SPAx <- SPAx$SolZen
  phi0_SPAx <- SPAx$SolAzi
  
  # remove pixels with zenith over 90, which is unphysical
  remove <- which(theta0_SPAx > 90, arr.ind = TRUE)
  theta0_SPAx[remove] <- 90
  
  # compute the differences
  Delta_theta0 <- theta0 - theta0_SPAx
  Delta_phi0 <- phi0 - phi0_SPAx
  Delta_phi0 <- pmin(Delta_phi0, 360 - Delta_phi0)
  
  # store them for plotting
  ras_Delta_theta0 <- raster(Delta_theta0, xmn = 1, xmx = 2748, ymn = 1, ymx = 2748)
  ras_Delta_theta0 <- aggregate(ras_Delta_theta0, fact = 4)
  ras_Delta_theta0 <- rasterToPoints(ras_Delta_theta0)
  ras_Delta_theta0 <- as_tibble(ras_Delta_theta0) %>%
    mutate(y = rev(y)) 
  ras_Delta_phi0 <- raster(Delta_phi0, xmn = 1, xmx = 2748, ymn = 1, ymx = 2748)
  #ras_azi <- aggregate(ras_azi, fact = 4)
  ras_Delta_phi0 <- rasterToPoints(ras_Delta_phi0)
  ras_Delta_phi0 <- as_tibble(ras_Delta_phi0) %>%
    mutate(y = rev(y)) 
  
  # arrange in tibble for plotting
  data.plot <- data.plot %>%
    bind_rows(., tibble(x = ras_Delta_theta0$x, y = ras_Delta_theta0$y, 
                        diff = ras_Delta_theta0$layer, version = i, 
                        angle = "Solar zenith")) %>%
    bind_rows(., tibble(x = ras_Delta_phi0$x, y = ras_Delta_phi0$y, 
                        diff = ras_Delta_phi0$layer, version = i, 
                        angle = "Solar azimuth"))
}

################################################################################
# plot zenith angle
################################################################################
min_index <- which.min(theta0) 
position <- arrayInd(min_index, dim(theta0))
# make a point to plot the position of the sun
data.point <- tibble(x = position[,2], y = position[,1])

ras_theta0 <- raster(theta0)
extent(ras_theta0) <- c(1, 2748, 1, 2748)  
ras_theta0 <- aggregate(ras_theta0, fact=4)
ras_theta0 <- rasterToPoints(ras_theta0)
ras_theta0 <- as_tibble(ras_theta0) %>%
  mutate(y = rev(y)) 

p1 <- ggplot(ras_theta0) + 
  geom_raster(aes(x = x, y = y, fill = layer)) +
  geom_point(data = data.point, aes(x = x, y = y), shape = 1, size = point.size)+
  coord_fixed() +
  viridis::scale_fill_viridis(name = expression(paste(italic(theta)[0], " [", degree, "]")), direction = -1, option = "B", limits = c(0,90), breaks = c(0, 30,60,90)) +
  scale_x_continuous("Column index", breaks = c(1, 1000, 2000, 2748)) +
  scale_y_reverse("Row index", breaks = c(1, 1000, 2000, 2748)) + #
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0,0), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(family = "Times", size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size/2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1,0,0.1,0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0,0.1,0,0.1), "lines")), panel.spacing = unit(0 , "lines"), legend.position = "bottom", legend.key.width = unit(1.2, "lines"), legend.key.height = unit(0.4, "lines"), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7,0,0,0.0), "lines"), legend.direction = "horizontal")


p2 <- ggplot(data.plot %>% filter(version == 1, angle == "Solar zenith")) + 
  geom_raster(aes(x = x, y = y, fill = diff)) +
  geom_point(data = data.point, aes(x = x, y = y), shape = 1, size = point.size)+
  coord_fixed() +
  scale_fill_gradient2(name = expression(paste(Delta, theta[0]," [", degree, "]")), low = brewer.pal(11, "RdBu")[11:10], mid = "white", high = brewer.pal(11, "RdBu")[2:1], midpoint = 0) +
  scale_x_continuous("Column index", breaks = c(1, 1000, 2000, 2748)) +
  scale_y_reverse("Row index", breaks = c(1, 1000, 2000, 2748)) + #
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0,0), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(family = "Times", size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size/2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1,0,0.1,0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0,0.1,0,0.1), "lines")), panel.spacing = unit(0 , "lines"), legend.position = "bottom", legend.key.width = unit(1.2, "lines"), legend.key.height = unit(0.4, "lines"), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7,1.8,0,0.0), "lines"), legend.direction = "horizontal")

################################################################################
# quantify the effect on extraterrestrial irradiance
################################################################################

# compute extraterrestrial irradiance using both versions of zenith angles obtained from SPA and SPA1
doy <- insol::daydoy(datetime)
da <- (2 * pi/365) * (doy - 1)
re <- 1.00011 + 0.034221 * cos(da) + 0.00128 * sin(da) + 0.00719 * 
  cos(2 * da) + 7.7e-05 * sin(2 * da)
E0_SPA <- round(1361.1 * re * cos(radians(theta0)))
theta0_SPA1 <- sun_angles1$SolZen
remove <- which(theta0_SPA1 > 90, arr.ind = TRUE)
theta0_SPA1[remove] <- 90
E0_SPA1 <- round(1361.1 * re * cos(radians(theta0_SPA1)))
Delta_E0 <- E0_SPA - E0_SPA1 

ras_Delta_E0 <- raster(Delta_E0)
extent(ras_Delta_E0) <- c(1, 2748, 1, 2748)  
ras_Delta_E0 <- aggregate(ras_Delta_E0, fact=4)
ras_Delta_E0 <- rasterToPoints(ras_Delta_E0)
ras_Delta_E0 <- as_tibble(ras_Delta_E0) %>%
  mutate(y = rev(y)) 

p3 <- ggplot(ras_Delta_E0) + 
  geom_raster(aes(x = x, y = y, fill = layer)) +
  geom_point(data = data.point, aes(x = x, y = y), shape = 1, size = point.size)+
  coord_fixed() +
  scale_fill_gradient2(name = expression(paste(Delta, italic(E)[italic(0)], " [W/", m^2, "]")), low = brewer.pal(11, "PiYG")[11:10], mid = "white", high = brewer.pal(11, "PiYG")[2:1], midpoint = 0, limits = c(-50,50), breaks = c(-50, -25, 0, 25, 50)) +
  scale_x_continuous("Column index", breaks = c(1, 1000, 2000, 2748)) +
  scale_y_reverse("Row index", breaks = c(1, 1000, 2000, 2748)) + #
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0,0), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(family = "Times", size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size/2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1,0,0.1,0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0,0.1,0,0.1), "lines")), panel.spacing = unit(0 , "lines"), legend.position = "bottom", legend.key.width = unit(1.2, "lines"), legend.key.height = unit(0.4, "lines"), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7,2,0,0.0), "lines"), legend.direction = "horizontal", legend.title = element_text(margin = margin(r = 10)))

pA <- ggpubr::ggarrange(p1, p2, p3, nrow = 1, labels = c("(a)", "(b)", "(c)"), widths = c(1, 1, 1), font.label = list(size = plot.size, color = "black", face = "plain", family = "Times"))

################################################################################
# plot azimuth angle
################################################################################

ras_phi0 <- raster(phi0)
extent(ras_phi0) <- c(1, 2748, 1, 2748)
ras_phi0 <- aggregate(ras_phi0, fact = 4)
ras_phi0 <- rasterToPoints(ras_phi0)
ras_phi0 <- as_tibble(ras_phi0) %>%
  mutate(y = rev(y)) 

p4 <- ggplot(ras_phi0) + 
  geom_raster(aes(x = x, y = y, fill = layer)) +
  geom_point(data = data.point, aes(x = x, y = y), shape = 1, size = point.size)+
  coord_fixed() +
  viridis::scale_fill_viridis(name = expression(paste(phi[0], " [", degree, "]")), direction = -1, option = "D", limits = c(0,360), breaks = c(0, 100, 200, 300)) +
  scale_x_continuous("Column index", breaks = c(1, 1000, 2000, 2748)) +
  scale_y_reverse("Row index", breaks = c(1, 1000, 2000, 2748)) + #
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0,0), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(family = "Times", size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size/2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1,0,0.1,0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0,0.1,0,0.1), "lines")), panel.spacing = unit(0 , "lines"), legend.position = "bottom", legend.key.width = unit(1.2, "lines"), legend.key.height = unit(0.4, "lines"), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7,0,0,0.0), "lines"), legend.direction = "horizontal")

p5 <- ggplot(data.plot %>% filter(version == 1, angle == "Solar azimuth")) + 
  geom_raster(aes(x = x, y = y, fill = diff)) +
  geom_point(data = data.point, aes(x = x, y = y), shape = 1, size = point.size)+
  coord_fixed() +
  scale_fill_gradient2(name = expression(paste(Delta, phi[0], " [", degree, "]")), low = brewer.pal(11, "RdBu")[11:10], mid = "white", high = brewer.pal(11, "RdBu")[2:1], midpoint = 0, transform = "pseudo_log", limits = c(-180,180), breaks = c(-100, -10, 0, 10, 100)) +
  scale_x_continuous("Column index", breaks = c(1, 1000, 2000, 2748)) +
  scale_y_reverse("Row index", breaks = c(1, 1000, 2000, 2748)) + #
  theme_bw() +
  theme(plot.margin = unit(c(0.2,0.2,0,0), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(family = "Times", size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size/2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1,0,0.1,0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0,0.1,0,0.1), "lines")), panel.spacing = unit(0 , "lines"), legend.position = "bottom", legend.key.width = unit(1.2, "lines"), legend.key.height = unit(0.4, "lines"), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7,1,0,0.0), "lines"), legend.direction = "horizontal", legend.title = element_text(margin = margin(r = 10)))

pB <- ggpubr::ggarrange(p4, p5, ncol = 1, labels = c("(d)", "(e)"), heights = c(1, 1), font.label = list(size = plot.size, color = "black", face = "plain", family = "Times"))

################################################################################
# plot zenith angle, using other inputs, i.e., SPA2~SPA5
################################################################################

tmp <- data.plot %>% 
  filter(version == 2, angle == "Solar zenith") 
# Calculate quantile breaks
no_classes <- 10
quantiles <- quantile(tmp$diff, probs = seq(0.003, 0.997, length.out = no_classes + 1), na.rm = TRUE)

p6 <- ggplot(data.plot %>% 
               filter(version %in% c(2:5), angle == "Solar zenith") %>% 
               mutate(diff > quantiles[length(quantiles)], quantiles[length(quantiles)], diff) %>%
               mutate(diff < quantiles[1], quantiles[1], diff)) + 
  geom_raster(aes(x = x, y = y, fill = diff)) +
  geom_point(data = data.point, aes(x = x, y = y), shape = 1, size = point.size)+
  facet_wrap(~version, labeller = labeller(version = c("2" = "SPA 2", "3" = "SPA 3", "4" = "SPA 4", "5" = "SPA 5"))) +
  coord_fixed() +
  scale_fill_gradient2(name = expression(paste(Delta, theta[0], " [", degree, "]")), low = brewer.pal(11, "RdBu")[11:8], mid = "white", high = brewer.pal(11, "RdBu")[4:1], midpoint = 0, limits = c(-0.02, 0.02)) +
  scale_x_continuous("Column index", breaks = c(1, 1000, 2000, 2748)) +
  scale_y_reverse("Row index", breaks = c(1, 1000, 2000, 2748)) + #
  theme_bw() +
  theme(plot.margin = unit(c(1, 1, 1, 1), "lines"), text = element_text(family = "Times", size = plot.size), axis.text = element_text(family = "Times", size = plot.size), panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "gray80", linewidth = line.size/2), panel.background = element_rect(fill = "transparent"), plot.background = element_rect(fill = "transparent", color = NA), strip.text.x = element_text(size = plot.size, margin = unit(c(0.1,0,0.1,0), "lines")), strip.text.y = element_text(size = plot.size, margin = unit(c(0,0.1,0,0.1), "lines")), panel.spacing = unit(0.5, "lines"), legend.position = "bottom", legend.key.width = unit(1.1, "lines"), legend.key.height = unit(0.4, "lines"), legend.text = element_text(size = plot.size), legend.background = element_rect(fill = "transparent"), legend.key = element_rect(fill = "transparent"), legend.box.margin = unit(c(-0.7,0,0,0), "lines"), legend.direction = "horizontal", legend.title = element_text(margin = margin(r = 15)))

pC <- ggpubr::ggarrange(p6, ncol = 1, labels = c("(f)"), font.label = list(size = plot.size, color = "black", face = "plain", family = "Times"))

multiplot <- list(pA, pB, pC)

lay <- rbind(c(1, 1, 1),
             c(2, 3, 3),
             c(2, 3, 3))

p <- gridExtra::grid.arrange(grobs = multiplot, layout_matrix = lay)

setwd(file.path(dir0, "tex"))
ggsave(filename = "zenith.pdf", plot = p, width = 160, height = 160, unit = "mm")

