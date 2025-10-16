################################################################################
# This code is written by Dazhi Yang
# School of Electrical Engineering and Automation
# Harbin Institute of Technology
# emails: yangdazhi.nus@gmail.com
################################################################################

################################################################################
# The following functions are related to reading and writing of (processed) Fengyun L1 data
################################################################################

# get nominal longitudes and latitudes values from "FY-4B数据行列号和经纬度查找表(105°E)"
get_lon_lat <- function(res = "4000M", directory = "/Users/seryangd/Library/CloudStorage/Dropbox/Working papers/Angles/Data")
{
  ######################################################################
  # INPUT
  # res (character string): resolution of the L1b data file, "0500M", "1000M", "2000M", "4000M"
  # directory (character string): the directory that holds the .raw files from http://www.nsmc.org.cn/nsmc/cn/satellite/FY4B.html
  
  # OUTPUT
  # lon (numeric): longitudes of each pixel, in matrix format
  # lat (numierc): latitudes of each pixel, in matrix format
  ######################################################################
  setwd(directory)
  file <- switch(res, 
                 "0500M" = "FY4B-_DISK_1050E_GEO_NOM_LUT_20240227000000_0500M_V0001.raw",
                 "1000M" = "FY4B-_DISK_1050E_GEO_NOM_LUT_20240227000000_1000M_V0001.raw",
                 "2000M" = "FY4B-_DISK_1050E_GEO_NOM_LUT_20240227000000_2000M_V0001.raw",
                 "4000M" = "FY4B-_DISK_1050E_GEO_NOM_LUT_20240227000000_4000M_V0001.raw")
  gridline <- switch(res, "0500M" = 21984, "1000M" = 10992, "2000M" = 5496, "4000M" = 2748)
  
  # read the binary file
  fid <- file(file, "rb")
  # read the data as doubles with little-endian byte order
  # note: R reads by column by default
  data <- readBin(fid, what = "double", n = 2 * gridline * gridline, size = 8, endian = "little")
  data <- matrix(data, nrow = 2*gridline, ncol = gridline, byrow = FALSE)
  # close the file
  close(fid)
  
  # extract latitude and longitude
  # MATLAB's data(1:2:end,:) equivalent in R
  lat <- data[seq(1, nrow(data), by = 2), ]  # Latitude (纬度)
  lon <- data[seq(2, nrow(data), by = 2), ]  # Longitude (经度)
  lat <- ifelse(lat > 1000, NA, lat) # remove missing values (圆盘外)
  lon <- ifelse(lon > 1000, NA, lon) # remove missing values (圆盘外)
  
  # rotate the matrix, such that the latitude and longitude follows the Earth
  lat <- t(lat)
  lon <- t(lon)
  
  return(list(lon = lon, lat = lat))
}

# read scanning time of AGRI
get_row_scan_time <- function(datetime, directory = "/Volumes/Macintosh Research/Data/FY-4B/L1")
{
  ######################################################################
  # INPUT
  # datetime (dttm): the STARTING time stamp of the desired FDI file
  # directory (character string): the directory that holds the .HDF files
  
  # OUTPUT
  # time_start (dttm): the starting time for each row scan
  # time_end (dttm): the ending time for each row scan
  ######################################################################
  
  if(missing(datetime))
  {
    stop("Argument 'datetime' must be provided!")
  }
  
  # construct file name
  time_start <- format(datetime, "%Y%m%d%H%M%S")
  time_end <- format(datetime + 14*60+59, "%Y%m%d%H%M%S")
  file <- paste0("FY4B-_AGRI--_N_DISK_1050E_L1-_FDI-_MULT_NOM_", time_start, "_", time_end, "_4000M_V0001.HDF")
  
  require("raster")
  
  # read file and extract data
  setwd(directory)
  
  # MONObsTime is store in another group in HDF5, so read it separately
  NOMObsTime <- raster::raster(paste0('HDF5:"', file, '"://NOMObs/NOMObsTime'))
  NOMObsTime <- flip(NOMObsTime, "y")
  
  # now convert it to two dttm objects, representing the start and end time for each row scan
  NOMObsTime <- as.matrix(NOMObsTime)
  time_start <- as.character(NOMObsTime[,1])
  second <- paste0(substr(time_start, 13, 14), ".", substr(time_start, 15, 17))
  time_start <- paste0(substr(time_start, 1, 12), second)
  time_start <- strptime(time_start, "%Y%m%d%H%M%OS", tz = "UTC") 
  time_end <- as.character(NOMObsTime[,2])
  second <- paste0(substr(time_end, 13, 14), ".", substr(time_end, 15, 17))
  time_end <- paste0(substr(time_end, 1, 12), second)
  time_end <- strptime(time_end, "%Y%m%d%H%M%OS", tz = "UTC") 
  
  return(list(time_start = time_start, time_end = time_end))
}

# allocate the exact scanning time to lon/lat positions, for accurate solar positioning
allocate_scan_time <- function(time_start, time_end, lat, lon) 
{
  ######################################################################
  # INPUT
  # time_start (dttm): the starting time for each row scan
  # time_end (dttm): the ending time for each row scan
  # lat (numeric): a matrix holding the latitudes of query points, as obtained from "get_lon_lat"
  # lon (numeric): a matrix holding the longitudes of query points, as obtained from "get_lon_lat"
  
  # OUTPUT
  # scan_time (numeric): scanning time of each pixel, in matrix format (use as.POSIXct(scan_time, origin = "1970-01-01") to convert back to dttm)
  ######################################################################
  
  # validate inputs
  if (length(time_start) != nrow(lat)) {
    stop("time_start length must match number of rows in latitude matrix")
  }
  if (length(time_end) != nrow(lat)) {
    stop("time_end length must match number of rows in latitude matrix")
  }
  
  # Initialize result matrix
  scan_time <- matrix(NA, nrow = nrow(lat), ncol = ncol(lat))
  
  for (i in 1:nrow(lat)) {
    # Skip rows with missing time information
    if (is.na(time_start[i]) || is.na(time_end[i])) {
      next
    }
    
    # Find valid columns (non-NA in both lat and lon)
    valid_cols <- which(!is.na(lat[i, ]) & !is.na(lon[i, ]))
    
    if (length(valid_cols) == 0) {
      next
    }
    
    # Convert to numeric for calculations
    start_num <- as.numeric(time_start[i])
    end_num <- as.numeric(time_end[i])
    row_duration <- end_num - start_num
    
    # Linear interpolation between first and last valid column
    if (length(valid_cols) > 1) {
      time_increment <- row_duration / (length(valid_cols) - 1)
      for (j in 1:length(valid_cols)) {
        time_offset <- (j - 1) * time_increment
        scan_time[i, valid_cols[j]] <- start_num + time_offset
      }
    } else {
      # Single valid column
      scan_time[i, valid_cols] <- start_num
    }
  }
  
  return(scan_time)
}

# main function for reading HDF
hdf_read <- function(datetime, name, directory = "/Users/seryangd/Library/CloudStorage/Dropbox/Working papers/Angles/Data")
{
  ######################################################################
  # INPUT
  # datetime (dttm): the STARTING time stamp of the desired file
  # name (character): the data name, which includes GEO and FDI
  # directory (character string): the directory that holds the .HDF files
  
  # OUTPUT
  # lon (numeric): longitudes of each pixel, in matrix format
  # lat (numierc): latitudes of each pixel, in matrix format
  ######################################################################
  
  if(missing(datetime) | missing(name))
  {
    stop("Argument 'datetime' and 'name' must be provided!")
  }

  # construct file name
  time_start <- format(datetime, "%Y%m%d%H%M%S")
  time_end <- format(datetime + 14*60+59, "%Y%m%d%H%M%S")
  file <- paste0("FY4B-_AGRI--_N_DISK_1050E_L1-_", name, "-_MULT_NOM_", time_start, "_", time_end, "_4000M_V0001.HDF")
  
  require("raster")
  
  # read file and extract data
  setwd(directory)
  
  # the file can be either GEO or FDI
  if(name == "GEO")
  {
    # 7 variables, with 5 being positioning angles
    vars <- c("ColumnNumber", "LineNumber", "NOMSatelliteAzimuth", "NOMSatelliteZenith", "NOMSunAzimuth", "NOMSunGlintAngle", "NOMSunZenith")
  }else if(name == "FDI")
  {
    # 15 channels and nominal observation time
    vars <- c(paste0("NOMChannel", sprintf("%02d", 1:15)))
  }
    
  # read the L1 data, and flip about equator (this file is as such, but need to check others)
  data <- raster::stack(file)
  data <- flip(data, direction = "y")
  
  for(i in 1:length(vars))
  {
    # get the variables, each being a raster
    assign(vars[i], raster::raster(data, layer = i)) 
    
    # store the SpatRaster into a temporary variable
    tmp <- get(vars[i]) 
    
    # convert tmp into a matrix, and remove the fill values
    if(vars[i] %in% c("ColumnNumber", "LineNumber"))
    {
      tmp[tmp == -1] <- NA
    }else{
      tmp[tmp == 65535] <- NA
    }
    
    # assign back to the corresponding variable name
    assign(vars[i], tmp)
  }
  
  # prepare output list
  out <- lapply(1:length(vars), function(x) get(vars[x]))
  names(out) <- vars

  return(out)
}

################################################################################
# The following functions are related to computing the angle parameters
################################################################################

# get digital elevation model at 0.05 deg resolution
get_dem <- function(lat, lon, directory = "/Users/seryangd/Library/CloudStorage/Dropbox/Working papers/Angles/Data")
{
  ######################################################################
  # INPUT
  # lat (numeric): a matrix holding the latitudes of query points, as obtained from "get_lon_lat"
  # lon (numeric): a matrix holding the longitudes of query points, as obtained from "get_lon_lat"
  # directory (character string): the directory that holds the DEM data, which is called world_ll_elev_0.05deg.nc4
  
  # OUTPUT
  # alt (numeric): altitude of each pixel, in matrix format
  ######################################################################
  
  require(ncdf4)
  require(RANN) # for fast match of coordinates
  
  setwd(directory)
  
  # read the DEM data
  ncin <- ncdf4::nc_open("world_ll_elev_0.05deg.nc4")
  alt <- ncdf4::ncvar_get(ncin)
  lon_alt <- ncin$var$elevation$dim[[1]]$vals
  lat_alt <- ncin$var$elevation$dim[[2]]$vals
  ncdf4::nc_close(ncin)
  
  # Create a data frame of all elevation grid points
  alt_grid <- expand.grid(lon = lon_alt, lat = lat_alt)
  alt_grid$alt <- as.vector(alt)
  
  # Flatten matrices holding the query points
  lat_vector <- as.vector(lat)
  lon_vector <- as.vector(lon)
  
  # Identify non-NA positions
  non_na_mask <- !is.na(lat_vector) & !is.na(lon_vector)
  
  # Convert non-NA target points to data frame for nn2 search
  query_points <- data.frame(
    lon = lon_vector[non_na_mask],
    lat = lat_vector[non_na_mask]
  )
  
  # Initialize result matrix with NAs
  alt_matrix <- matrix(NA, nrow = nrow(lat), ncol = ncol(lat))
  
  # fill the non_na locations with 0 first, for oceans
  alt_matrix[non_na_mask] <- 0
    
  # Only perform search if there are non-NA points
  if (nrow(query_points) > 0) {
    # Perform nearest neighbor search
    nn_result <- RANN::nn2(
      data = alt_grid[, c("lon", "lat")],
      query = query_points,
      k = 1
    )
    
    # Extract elevation values using the indices
    alt_values <- alt_grid$alt[nn_result$nn.idx]
    
    # Assign values to non-NA positions in the result matrix
    alt_matrix[non_na_mask] <- alt_values
  }
  
  # output a matrix same size as lat (lon)
  return(alt_matrix)
}

# get air pressure from monthly ERA5 reanalysis
get_sp <- function(lat, lon, directory = "/Users/seryangd/Library/CloudStorage/Dropbox/Working papers/Angles/Data")
{
  ######################################################################
  # INPUT
  # lat (numeric): a matrix holding the latitudes of query points, as obtained from "get_lon_lat"
  # lon (numeric): a matrix holding the longitudes of query points, as obtained from "get_lon_lat"
  # directory (character string): the directory that holds the monthly ERA5 data. 
  
  # OUTPUT
  # sp (numeric): surface pressure of each pixel, in matrix format (unit: Pa)
  ######################################################################
  # NOTE: ERA5 data can be downloaded from https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-monthly-means?tab=overview, 
  # NOTE: The format should be in NetCDF, and should contain the monthly surface pressure (sp) and air temperature (t2m) over a full year
  # NOTE: In other words, all sp and t2m data should be held in a single .nc file
  
  require(ncdf4)
  require(RANN) # for fast match of coordinates
  
  setwd(directory)
  
  # read the ERA5 monthly data
  ncin <- ncdf4::nc_open("2024monthly.nc")
  sp <- ncdf4::ncvar_get(ncin, "sp")
  lon_sp <- ncdf4::ncvar_get(ncin, "longitude")
  lat_sp <- ncdf4::ncvar_get(ncin, "latitude")
  ncdf4::nc_close(ncin)
  
  # average 12 months data for each lat/lon point
  sp <- apply(sp, c(1, 2), mean, na.rm = TRUE)
  
  # Create a data frame of all surface pressure grid points
  sp_grid <- expand.grid(lon = lon_sp, lat = lat_sp)
  sp_grid$sp <- as.vector(sp)
  
  # Flatten matrices holding the query points
  lat_vector <- as.vector(lat)
  lon_vector <- as.vector(lon)
  
  # Identify non-NA positions
  non_na_mask <- !is.na(lat_vector) & !is.na(lon_vector)
  
  # Convert non-NA target points to data frame for nn2 search
  query_points <- data.frame(
    lon = lon_vector[non_na_mask],
    lat = lat_vector[non_na_mask]
  )
  
  # Initialize result matrix with NAs
  sp_matrix <- matrix(NA, nrow = nrow(lat), ncol = ncol(lat))
  
  # fill the non_na locations with 0 first, for oceans
  sp_matrix[non_na_mask] <- 0
  
  # Only perform search if there are non-NA points
  if (nrow(query_points) > 0) {
    # Perform nearest neighbor search
    nn_result <- RANN::nn2(
      data = sp_grid[, c("lon", "lat")],
      query = query_points,
      k = 1
    )
    
    # Extract elevation values using the indices
    sp_values <- sp_grid$sp[nn_result$nn.idx]
    
    # Assign values to non-NA positions in the result matrix
    sp_matrix[non_na_mask] <- sp_values
  }
  
  # output a matrix same size as lat (lon)
  return(sp_matrix)
}

# get air temperature from monthly ERA5 reanalysis
get_t2m <- function(lat, lon, directory = "/Users/seryangd/Library/CloudStorage/Dropbox/Working papers/Angles/Data")
{
  ######################################################################
  # INPUT
  # lat (numeric): a matrix holding the latitudes of query points, as obtained from "get_lon_lat"
  # lon (numeric): a matrix holding the longitudes of query points, as obtained from "get_lon_lat"
  # directory (character string): the directory that holds the monthly ERA5 data. 
  
  # OUTPUT
  # t2m (numeric): 2-m air temperature of each pixel, in matrix format (unit: deg C)
  ######################################################################
  # NOTE: ERA5 data can be downloaded from https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-monthly-means?tab=overview, 
  # NOTE: The format should be in NetCDF, and should contain the monthly surface pressure (sp) and air temperature (t2m) over a full year
  # NOTE: In other words, all sp/t2m data should be held in a single .nc file
  
  require(ncdf4)
  require(RANN) # for fast match of coordinates
  
  setwd(directory)
  
  # read the ERA5 monthly data
  ncin <- ncdf4::nc_open("2024monthly.nc")
  t2m <- ncdf4::ncvar_get(ncin, "t2m")
  lon_t2m <- ncdf4::ncvar_get(ncin, "longitude")
  lat_t2m <- ncdf4::ncvar_get(ncin, "latitude")
  ncdf4::nc_close(ncin)
  
  # average 12 months data for each lat/lon point
  t2m <- apply(t2m, c(1, 2), mean, na.rm = TRUE) 
  
  # convert from Kelvin to Celsius degree
  t2m <- t2m - 273.15
  
  # Create a data frame of all surface pressure grid points
  t2m_grid <- expand.grid(lon = lon_t2m, lat = lat_t2m)
  t2m_grid$t2m <- as.vector(t2m)
  
  # Flatten matrices holding the query points
  lat_vector <- as.vector(lat)
  lon_vector <- as.vector(lon)
  
  # Identify non-NA positions
  non_na_mask <- !is.na(lat_vector) & !is.na(lon_vector)
  
  # Convert non-NA target points to data frame for nn2 search
  query_points <- data.frame(
    lon = lon_vector[non_na_mask],
    lat = lat_vector[non_na_mask]
  )
  
  # Initialize result matrix with NAs
  t2m_matrix <- matrix(NA, nrow = nrow(lat), ncol = ncol(lat))
  
  # fill the non_na locations with 0 first, for oceans
  t2m_matrix[non_na_mask] <- 0
  
  # Only perform search if there are non-NA points
  if (nrow(query_points) > 0) {
    # Perform nearest neighbor search
    nn_result <- RANN::nn2(
      data = t2m_grid[, c("lon", "lat")],
      query = query_points,
      k = 1
    )
    
    # Extract elevation values using the indices
    t2m_values <- t2m_grid$t2m[nn_result$nn.idx]
    
    # Assign values to non-NA positions in the result matrix
    t2m_matrix[non_na_mask] <- t2m_values
  }
  
  # output a matrix same size as lat (lon)
  return(t2m_matrix)
}

# get deltaT from USNO delta T: https://maia.usno.navy.mil/products/deltaT
get_deltaT <- function(datetime, directory = "/Users/seryangd/Library/CloudStorage/Dropbox/Working papers/Angles/Data")
{
  ######################################################################
  # INPUT
  # datetime (dttm): Datetime format, which can be arrived at using lubridate::ymd_hm or similar
  # directory (character string): the directory that holds the USNO delta T .txt file. 
  
  # OUTPUT
  # deltaT (numeric): difference between terrestrial time and UT1.
  ######################################################################
  
  require("dplyr")
  
  setwd(directory)
  
  # read the list of deltaT values from file
  deltaT_all <- tibble(read.table("deltat_since1973.txt", header = FALSE))
  # manipulate the tibble for easy access
  deltaT_all <- deltaT_all %>%
    rename(year = V1, month = V2, day = V3, deltaT = V4) %>%
    mutate(date = paste(year, month, day, sep = "-")) %>%
    mutate(date = ymd(date)) %>%
    dplyr::select(one_of("date", "deltaT"))
  
  # change the datetime to date formate
  matched_index <- which(deltaT_all$date == floor_date(datetime, unit = "1 mon"))
  if(length(matched_index) == 0)
  {
    deltaT <- 0.69
  }else{
    deltaT <- deltaT_all$deltaT[matched_index]
  }
  
  return(deltaT)
}

# calculate solar zenith and azimuth angles in R solarPos package
# angle calculation follows SPA by Reda and Andreas, 2004, Solar Energy, 76(5), 577-589.
solpos <- function(lat, lon, alt = 0, sp = 1013.25, t2m = 16, deltaT = 69, datetime, tz = 0) 
{
  ######################################################################
  # INPUT
  # lat (numeric): a matrix holding the latitudes
  # lon (numeric): a matrix holding the longitudes
  # alt (numierc): a matrix holding the altitudes (m)
  # sp (numeric): a matrix holding the surface pressure (hPa)
  # t2m (numierc): a matrix holding the air temperature (degC)
  # deltaT (numeric): a scalar holding the deltaT value, see pvlib documentation
  # datetime (dttm): Datetime format, which can be arrived at using lubridate::ymd_hm or similar
  # tz (numeric): time zone, but since Fengyun data is in UTC, default it to 0 is fine
  
  # OUTPUT
  # SolZen (numeric): solar zenith angle in degrees
  # SolAzi (numeric): solar azimuth angle in degrees, zero north east positive, range 0 -- 360
  ######################################################################
  
  # require the solarPos package for solar positioning and the lubridate package to handle datetime data
  require("solarPos")
  require("lubridate")
  
  jd <- solarPos::julianDay(lubridate::year(datetime), lubridate::month(datetime), 
                            lubridate::day(datetime), lubridate::hour(datetime), 
                            lubridate::minute(datetime), lubridate::second(datetime), tz = tz)

  if(length(jd) == length(lat))
  {
    # convert lat, lon, ..., sp to vector
    jd_mat <- matrix(jd, nrow = nrow(lat), ncol = ncol(lat))
    
    # calculate mean of non-NA values
    jd_mean <- mean(jd_mat, na.rm = TRUE)
    
    # replace NAs with the mean (otherwise the vector computation in solarPosition will fail)
    jd_mat[is.na(jd_mat)] <- jd_mean
    
    # change the name back to jd
    jd <- jd_mat
  }
  
  # this code auto convert matrix input to vector, and output vector
  result <- solarPos::solarPosition(jd = jd, lon = lon, lat = lat, elev = alt, 
                                    delta_t = deltaT, temp = t2m, pres = sp)
  
  # convert azi and zen back to matrices
  theta0 <- matrix(result[,'zenith'], nrow = nrow(lat), ncol = ncol(lat))
  phi0 <- matrix(result[,'azimuth'], nrow = nrow(lat), ncol = ncol(lat))
  
  return(list(SolZen = theta0, SolAzi = phi0))
}

################################################################################
# The following functions are related to the calculation of satellite angles
################################################################################

satpos <- function (lat, lon, alt = 0, ar_cor = FALSE, sp = 1013.25, t2m = 16) 
{
  ######################################################################
  # INPUT
  # lat (numeric): a matrix holding the latitudes
  # lon (numeric): a matrix holding the longitudes
  # alt (numierc): a matrix holding the altitudes (m)
  # sp (numeric): a matrix holding the surface pressure (hPa)
  # t2m (numierc): a matrix holding the air temperature (degC)
  
  # OUTPUT
  # SatZen (numeric): view zenith angle, in degrees 卫星天顶角
  # SatAzi (numeric): view azimuth angle, in degrees 卫星方位角
  ######################################################################
  
  radians <- function (degree) {degree * (pi/180)}
  degrees <- function (radian) {radian * (180/pi)}
  
  # lon lat in radians
  lambda <- radians(lon)
  varphi <- radians(lat)
  
  # computation according to https://geodesy.noaa.gov/CORS/Articles/SolerEisemannJSE.pdf
  h <- alt # 高度
  a <- 6378137  # 地球在赤道的半径 (m)
  f <- 1/298.257222101 # flattening
  hs <- 35786000 # 传感器到地表的距离 (m)
  e2 <- 2*f-f^2 # 椭球偏心率的平方
  N <- a/sqrt(1-e2*(sin(varphi))^2) # 主垂直曲率半径
  lambdas <- radians(105) # 风云点位经度, 弧度
  varphis <- radians(0) # 风云点位纬度, 弧度
  
  # satellite’s position in the geocentric Cartesian coordinate system
  xs <- (a + hs)*cos(lambdas)
  ys <- (a + hs)*sin(lambdas)
  zs <- 0
  
  #  ground point’s coordinates in the geocentric Cartesian system
  xp <- (N + h)*cos(varphi)*cos(lambda)
  yp <- (N + h)*cos(varphi)*sin(lambda)
  zp <- ((1-e2)*N + h)*sin(varphi)
  
  # direction vector to satellite
  x <- xs - xp
  y <- ys - yp
  z <- zs - zp
  
  # compute e, n, and u
  e_vec <- -sin(lambda)*x + cos(lambda)*y
  n_vec <- -sin(varphi)*cos(lambda)*x - sin(varphi)*sin(lambda)*y + cos(varphi)*z
  u_vec <- cos(varphi)*cos(lambda)*x + cos(varphi)*sin(lambda)*y+ sin(varphi)*z
    
  # compute viewing geometry angles
  theta <- 90 - degrees(atan2(u_vec, sqrt(e_vec^2 + n_vec^2)))
  phi <- degrees(atan2(e_vec, n_vec))
  phi <- ifelse(phi < 0, 360 + phi, phi) # convert to 0-360
  
  if(ar_cor)
  {
    gamma <- 90 - theta
    delta_gamma <- sp/1010 * 283/(273+t2m) * 1.02/(60*tan(radians(gamma + 10.3/(gamma+5.1))))
    delta_gamma <-ifelse(gamma > - 0.83337, delta_gamma, 0)
    theta <- 90 - gamma - delta_gamma
  }
  
  return(list(SatZen = theta, SatAzi = phi))
}


################################################################################
# The following function calculates the sun glint angle
################################################################################

sunglint <- function(SolZen, SatZen, SolAzi, SatAzi) 
{
  ######################################################################
  # INPUT
  # SolZen (numeric): a matrix holding the solar zenith angles
  # SatZen (numeric): a matrix holding the viewing zenith angles
  # SolAzi (numierc): a matrix holding the solar azimuth angles
  # SatAzi (numeric): a matrix holding the viewing azimuth angles
  
  # OUTPUT
  # SunGlint (numeric): sun glint angle, in degrees 天阳耀斑角
  ######################################################################
  
  radians <- function (degree) {degree * (pi/180)}
  degrees <- function (radian) {radian * (180/pi)}
  
  # Calculating sunglint angle
  Temp1 <- cos(radians(SolZen)) * cos(radians(SatZen))
  Temp2 <- sin(radians(SolZen)) * sin(radians(SatZen))
  Temp3 <- cos(radians(180.0 - (SolAzi - SatAzi)))
  
  SunGlintA <- degrees(acos(Temp1 + Temp2 * Temp3))
  return(SunGlintA)
}




