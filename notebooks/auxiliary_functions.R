computeTrend = function(ts) {  
  # Function to compute temporal trends
  df = data.frame(x = 1:length(ts), y = ts)
  if (sum(is.na(ts)) < round(0.75*length(ts))) {  # ask for a minimum of 75% of non-missing data to compute the trend
    reg = lm(y ~ x, df)
    return(reg$coefficients[2])
  } else {
    return(NA)
  }
}

computePvalTrend = function(ts) {  
  # Function to compute the p-value of temporal trends
  df = data.frame(x = 1:length(ts), y = ts)
  if (sum(is.na(ts)) < round(0.75*length(ts))) {  # ask for a minimum of 75% of non-missing data to compute the trend
    reg = lm(y ~ x, df)
    return(summary.lm(reg)$coefficients[2,4])
  } else {
    return(NA)
  }
}

computePSS <- function(data1, data2, bins = 100) {
  # Function to compute the Perkins Skill Score (PSS)
  
  ## removing NA values (for safety)
  data1 <- na.omit(as.numeric(data1))
  data2 <- na.omit(as.numeric(data2))
  
  ## define global common limits so bins align perfectly
  min_limit <- min(c(data1, data2))
  max_limit <- max(c(data1, data2))
  
  ## create common break points for the histogram
  breaks <- seq(min_limit, max_limit, length.out = bins + 1)
  
  ## compute the histogram for each dataset using the same breaks
  ## setting plot = FALSE prevents R from drawing the charts during calculation
  hist1 <- hist(data1, breaks = breaks, plot = FALSE)
  hist2 <- hist(data2, breaks = breaks, plot = FALSE)
  
  ## convert absolute frequencies into relative frequencies (empirical PDF)
  ## the sum of relative frequencies for each dataset equals 1
  rel_freq1 <- hist1$counts / sum(hist1$counts)
  rel_freq2 <- hist2$counts / sum(hist2$counts)
  
  ## calculate the Perkins Score by summing the minimum frequency per bin
  perkins_score <- sum(pmin(rel_freq1, rel_freq2))
  
  ## return the final numeric score
  return(perkins_score)
}

maskGrid <- function(grid, shape, crs = 4326) {
  # Function to mask data to a specific region, 
  # delimited by a shapefile
  
  ## reading shapefile
  if (is.character(shape)) {
    shape <- st_read(shape, quiet = TRUE)
  }
  
  ## transform to the grid's CRS
  shape <- st_transform(shape, crs)
  
  lon <- grid$xyCoords$x
  lat <- grid$xyCoords$y
  
  ## creating spatial grid
  pts <- expand.grid(
    lon = lon,
    lat = lat
  )
  
  pts.sf <- st_as_sf(
    pts,
    coords = c("lon", "lat"),
    crs = crs
  )
  
  ## selecting points within the polygon
  inside <- lengths(st_within(pts.sf, shape)) > 0
  
  ## building mask (lat × lon)
  mask <- matrix(
    inside,
    nrow = length(lon),
    ncol = length(lat),
    byrow = FALSE
  )
  
  mask <- t(mask)
  
  ## applying mask
  out <- grid
  
  nt <- dim(grid$Data)[1]
  
  for (i in seq_len(nt)) {
    tmp <- out$Data[i, , ]
    tmp[!mask] <- NA
    out$Data[i, , ] <- tmp
  }
  
  return(out)
}

harmonizeDatesCMIP5 = function(grid, dates.ref) {
  # Function to eliminate repeated dates in a grid 
  # (data corresponding to repeated dates are averaged)
  
  grid.out = grid  # template
  grid.out$Dates$start = dates.ref
  grid.out$Dates$end = dates.ref
  grid.out$Data = array(NA, c(length(dates.ref), 
                              getShape(grid)["lat"],
                              getShape(grid)["lon"]))
  attributes(grid.out$Data)$dim = c(length(dates.ref),
                                    as.numeric(getShape(grid)["lat"]),
                                    as.numeric(getShape(grid)["lon"]))
  attributes(grid.out$Data)$dimensions = c("time", "lat", "lon")
  
  for (d in 1:length(dates.ref)) {
    ind = which(is.element(grid$Dates$start, dates.ref[d]))
    grid.out$Data[d, , ] = suppressMessages(climatology(subsetDimension(grid, 
                                                                        dimension = "time", 
                                                                        indices = ind))$Data)
  }
  return(grid.out)
}

plot_pdfs <- function(grid1, grid2, text1, text2, xlab = NULL, ylab = NULL, xlim = NULL) {
  # Function to plot two PDFs in the same figure
  pdf1 <- density(grid1$Data, na.rm = TRUE)
  pdf2 <- density(grid2$Data, na.rm = TRUE)

  plot(pdf1, col = "blue", lwd = 2,
       main = "PDFs",
       xlab = xlab, ylab = ylab,
       xlim = xlim)
  lines(pdf2, col = "red", lwd = 2)
  legend("topleft",
         legend = c(text1, text2),
         col = c("blue", "red"),
         lwd = 2)
  grid()
  invisible(NULL)
}

binData <- function(x, threshold, direction) {
  ## Function to convert a vector to 0/1 when a certain threshold is (is not) surpassed
  ## x: vector to binarize
  ## threshold: critical value for establishing the binarization
  ## direction: must be "GE" (greater or equal) or "LE" (lower or equal)
  
  x.bin = x
  
  if (direction == "GE") {
    x.bin[data >= threshold & !is.na(x)] = 1
    x.bin[data < threshold & !is.na(x)] = 0
  } else if (direction == "LE") {
    x.bin[data <= threshold & !is.na(x)] = 1
    x.bin[data > threshold & !is.na(x)] = 0
  }
  
  return(x.bin)
}

binGrid <- function(grid, threshold, direction) {
  ## Function to convert a vector to 0/1 when a certain threshold is (is not) surpassed
  ## data: vector to binarize
  ## threshold: critical value for establishing the binarization
  ## direction: must be "GE" (greater or equal) or "LE" (lower or equal)
  
  grid.bin = grid
  
  if (direction == "GE") {
    grid.bin$Data[grid$Data >= threshold & !is.na(grid$Data)] = 1
    grid.bin$Data[grid$Data < threshold & !is.na(grid$Data)] = 0
    
  } else if (direction == "LE") {
    grid.bin$Data[grid$Data <= threshold & !is.na(grid$Data)] = 1
    grid.bin$Data[grid$Data > threshold & !is.na(grid$Data)] = 0
  }
  
  attributes(grid.bin$Data) = attributes(grid$Data)
  return(grid.bin)
}
