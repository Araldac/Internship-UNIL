################################################################################
## Project: Intership- Relative influence invasion factor
## Author: Aitana Ralda Corfas
## Date: 05/11/2025
## Description: Functions of the project
################################################################################

################################################################################
#A-cleaning and thinning####
################################################################################

####	####	Cleaning ####	####
pkg <- c("magrittr", "CoordinateCleaner", "countrycode", "rnaturalearth", "sf", "sfheaders")
lapply(pkg, function(x) require(x, quietly = TRUE, character.only = TRUE))

load("Internship UNIL/Data/global_rasters/sf/sea.RData")	# oject w_gadm (in fact w_gadm_41): used to find points in the sea
load("Internship UNIL/Data/global_rasters/sf/lak.RData")	# object l_esri (in fact l_esri_100): used to find points in the lakes
load("Internship UNIL/Data/global_rasters/sf/land.RData")	# object land (in fact land f naturalearth): used to draw sea and land when mapping results
load("Internship UNIL/Data/global_rasters/sf/unland.RData")	# object unland (in fact l_esri_100): used to draw lakes when mapping results

"occ_cleaning" <- function(occ, lon = "decimallongitude", lat = "decimallatitude",
                           GBIF = TRUE, occ_precision = 0.01, occ_uncertainity1 = 50*1000,
                           occ_uncertainity2 = c(301,3036,999,9999), occ_count = 0,
                           occ_status = c("PRESENT"), occ_record = c("FOSSIL_SPECIMEN"),
                           occ_year = 1700, sea = NULL, lak = NULL, map = TRUE, details = TRUE){
  ## pipeline to clean occurrences of one species - the file could come from GBIF or not
  # occ: data.frame with two columns called "decimallongitude" and "decimallatitude" for coordinates in geographical units ; if no species field, a column will be added with an arbitrary species name
  # lon: name of the longitude column
  # lat: name of latitude column
  # GBIF: FALSE when occurrences do not come from GBIF and all the arguments are facultative
  # occ_precision: threshold linked to "coordinatePrecision" field ; if the field is absent, specify NULL
  # occ_uncertainity1: threshold linked to "coordinateUncertaintyInMeters" field ; if the field is absent, specify NULL
  # occ_uncertainity2: false values in "coordinateUncertaintyInMeters" field ; if the field is absent, specify NULL
  # occ_count: threshold linked to "individualCount" field : eliminate occurences if count is equal or under threshold ; in the field is absent, specify NULL
  # occ_status: attributes kept in the "occurrenceStatus" field ; in the field is absent, specify NULL
  # occ_record: attributes kept in the "basisOfRecord" field ; in the field is absent, specify NULL
  # occ_year: threshold linked to "year" field : eliminate occurences if year is under threshold ; in the field is absent, specify NULL
  # sea: if NULL, do not check if points fall in the sea or not; if TRUE, keep points in the sea using w_gam (in fact w_gadm_41); if FALSE, remove points in the sea using w_gadm (in fact w_gadm_41)
  # lak: if NULL, do not check if points fall in lakes or not; if TRUE, keep points in lakes using l_esri (in fact l_esri_100); if FALSE, remove points in lakes using l_esri (in fact l_esri_100)
  # map: if TRUE, map deleted occurences by filter
  # details: if TRUE, return rem, kip and ind ; else return directly index (ind) of points to keep
  
  #### Prepare and check dataset
  ### occ must be a data.frame
  if(!(is.data.frame(occ) | is.matrix(occ)))
    return(NULL)
  
  if(!is.matrix(occ))
    occ <- as.data.frame(occ)
  nr <- nrow(occ)
  nam <- names(occ)	
  
  ### occ must have two columns named "decimallongitude", "decimallatitude"
  i1 <- which(names(occ) == lon)
  i2 <- which(names(occ) == lat)
  names(occ)[c(i1, i2)] <- c("decimallongitude", "decimallatitude")
  
  ### occ must have a species column with only one species name
  test <- c("species") %in% nam
  if(test){
    l <- unique(occ$species)
    if(length(l) > 1){
      print("occ must concern only one species: it is probably not the case for your data, you should check")
      species <- rep("no species name", nr)
      occ$species <- species
    }
  }	else{
    species <- rep("no species name", nr)
    occ <- cbind.data.frame(occ, species)
    nam <- c(nam, "species")
  }
  
  ### Assign facultative argument if GBIF == FALSE
  if(!GBIF)
    occ_precision <- occ_uncertainity1 <- occ_uncertainity2 <- occ_count <- occ_status <- occ_record <- occ_year <- NULL
  
  
  ### Check other column names according to argument
  ## All facultative arguments
  if(is.null(occ_uncertainity1) & is.null(occ_uncertainity2)){
    local1 <- NULL
  }	else{
    local1 <- "UNNULL"
  }
  arg <- list(occ_precision, local1, occ_count, occ_status, occ_record, occ_year)
  colnam <- c("coordinatePrecision", "coordinateUncertaintyInMeters", "individualCount", "occurrenceStatus", "basisOfRecord", "year")
  null <- unlist(lapply(arg, is.null))
  test <- colnam[!null] %in% nam
  if(any(!test)){
    txt <- paste(colnam[!null][!test], collapse = ", ")
    txt <- paste0("According to the arguments, you need to have column(s) named ", txt, "; be careful to letter case")
    stop(txt)
  }
  
  #### Unvalid and duplicated
  ### Index of points
  ind <- 1:nr
  
  ### Unvalid
  val <- cc_val(occ, lon = lon, lat = lat, value = "flagged", verbose = FALSE)
  ind <- ind[val]
  occ <- occ[val,]
  
  ### Duplicated
  dup <- cc_dupl(occ, lon = lon, lat = lat, value = "flagged", verbose = FALSE)
  ind <- ind[dup]
  occ <- occ[dup,]
  
  nr_new <- nrow(occ)
  
  #--> to be tested
  if(nr_new == 0){
    res <- NULL
    print("No more occurences after cleaning unvalid and duplicated points")
    return(res)
  }
  
  
  #### Intrinsic
  ### Zero coordinates
  nul <- suppressWarnings(cc_zero(occ, lon = lon, lat = lat, buffer = 0.001, value = "flagged", verbose = FALSE))
  
  ### Equal coordinates
  equ <- cc_equ(occ, lon = lon, lat = lat, value = "flagged", verbose = FALSE)
  
  ### Precision
  if(!is.null(occ_precision))
    pre <- (occ$coordinatePrecision < occ_precision) | (is.na(occ$coordinatePrecision))
  
  ### Uncertainity
  if(!is.null(occ_uncertainity1)){
    unc1 <- (occ$coordinateUncertaintyInMeters < occ_uncertainity1) | (is.na(occ$coordinateUncertaintyInMeters))
  }	else{
    unc1 <- rep(TRUE, nr_new)
  }
  
  if(!is.null(occ_uncertainity1)){
    unc2 <- !(occ$coordinateUncertaintyInMeters %in% occ_uncertainity2)
  }	else{
    unc2 <- rep(TRUE, nr_new)
  }
  
  if(!is.null(occ_uncertainity1) | !is.null(occ_uncertainity2))
    unc <- unc1 & unc2
  
  
  #### Extrinsic
  ### Occurrences that fall in the sea
  if(is.null(sea)){
    wat1 <- rep(TRUE, nr_new)
    names(wat1) <- row.names(occ)
  }	else{
    lonlat <- sf_point(occ[,c("decimallongitude", "decimallatitude")])
    lonlat <- lonlat$geometry
    st_crs(lonlat) <- st_crs(w_gadm)
    
    local <- suppressMessages(st_intersects(lonlat, w_gadm))
    local <- unclass(local)				
    wat1 <- (unlist(lapply(local, length)) == 0)
    names(wat1) <- row.names(occ)
    
    if(!sea)	wat1 <- !wat1		
  }
  
  ### Occurrences that fall in lakes
  if(is.null(lak)){
    wat2 <- rep(TRUE, nr_new)
    names(wat2) <- row.names(occ)
  }	else{
    lonlat <- sf_point(occ[,c("decimallongitude", "decimallatitude")])
    lonlat <- lonlat$geometry
    st_crs(lonlat) <- st_crs(l_esri)
    
    local <- suppressMessages(st_intersects(lonlat, l_esri))
    local <- unclass(local)				
    wat2 <- !(unlist(lapply(local, length)) == 0)
    names(wat2) <- row.names(occ)
    
    if(!lak)	wat2 <- !wat2		
  }
  
  ### Vincinity of country capitals: error in the code (limits <- raster::extent(dat) + buffer #--> buffer is is meters and not in lonlat: what is the meaning of the crop step? buffer has been changed in 250 (it was 100 before the 22/03/2023)
  cap <- suppressWarnings(cc_cap(occ, lon = lon, lat = lat, buffer = 250, value = "flagged", verbose = FALSE))
  
  ### Vincinity of country and province centroids
  cen <- suppressWarnings(cc_cen(occ, lon = lon, lat = lat, buffer = 250, value = "flagged", verbose = FALSE))
  
  ### Vincinity of biodiversity institutions
  ins <- suppressWarnings(cc_inst(occ, lon = lon, lat = lat, buffer = 250, value = "flagged", verbose = FALSE))
  
  ### GBIF headquarters
  suppressWarnings(gbi <- cc_gbif(occ, lon = lon, lat = lat, buffer = 250, value = "flagged", verbose = FALSE))
  
  
  #### Attributes
  ### Status
  if(!is.null(occ_count)){
    sta1 <- (occ$individualCount > occ_count) | is.na(occ$individualCount)
  }	else{
    sta1 <- rep(TRUE, nr_new)
  }
  
  if(!is.null(occ_status)){
    sta2 <- (occ$occurrenceStatus %in% occ_status) | is.na(occ$occurrenceStatus)
  }	else{
    sta2 <- rep(TRUE, nr_new)
  }
  
  if(!is.null(occ_count) | !is.null(occ_status)){
    sta <- sta1 & sta2
  }
  
  ### Record
  if(!is.null(occ_record))
    rec <- !(occ$basisOfRecord %in% occ_record) | is.na(occ$basisOfRecord)
  
  
  ### Year
  if(!is.null(occ_year))
    dat <- (occ$year >= occ_year) | is.na(occ$year)
  
  
  #### Prepare results
  kip <- cbind.data.frame(nul, equ, cap, cen, gbi, wat1, wat2, ins)
  
  if(is.null(occ_uncertainity1) & is.null(occ_uncertainity2)){
    local1 <- NULL
  }	else{
    local1 <- "UNNULL"
  }
  
  if(is.null(occ_count) & is.null(occ_status)){
    local2 <- NULL
  }	else{
    local2 <- "UNNULL"
  }
  
  arg <- list(occ_precision, local1, local2, occ_record, occ_year)
  null <- unlist(lapply(arg, is.null))
  if(any(!null)){
    colnam <- c("pre", "unc", "sta", "rec", "dat")[!null]
    txt <- paste(colnam, collapse = ", ")
    exp <- paste0("cbind.data.frame(", txt, ")")
    df <- eval(parse(text = exp))
    kip <- cbind.data.frame(kip, df)
  }
  
  rem <- apply(!kip, 2, sum)
  rem <- c(sum(!val), sum(!dup), rem)
  names(rem)[1:2] <- c("val", "dup")
  
  
  #### Map results
  if(map){
    index <- apply(!kip, 2, sum)
    index <- which(index > 0)
    nc <- length(index) + 1
    namcol <- names(index)
    par(mfrow = n2mfrow(nc))
    par(mar = c(2.5, 0.5, 2.5, 0.5))
    
    ### Map kept points
    # plot(w_gadm, col = "grey", border = "grey")
    plot(land, col = "grey", border = "grey", bg = "cyan3")
    # plot(l_esri, col = "cyan1", border = "cyan1", add = T)
    plot(unland, col = "cyan1", border = "cyan1", add = T)
    poi <- which(apply(kip, 1, sum) == ncol(kip))
    points(occ[poi,c("decimallongitude", "decimallatitude")], col = "black", pch = "+", cex = 0.1)
    
    ### Map deleted points by filter
    for(i in namcol){
      # plot(w_gadm, col = "grey", border = "grey")
      plot(land, col = "grey", border = "grey", bg = "cyan3")
      # plot(l_esri, col = "cyan1", border = "cyan1", add = T)
      plot(unland, col = "cyan1", border = "cyan1", add = T)
      points(occ[!kip[,i],c("decimallongitude", "decimallatitude")], col = "red", pch = "+", cex = 0.5)
      txt <- paste(i, rem[i], sep = " - ")
      title(txt)
    }
  }
  
  
  #### Return
  if(details){
    res <- list(rem = rem, kip = kip, ind = ind)
  }	else{
    test <- apply(!kip, 1, sum)
    index <- which(test == 0)
    res <- ind[index]
  }
  return(res)
}

####	#### Thinning	####	#### 
#"splitcloud: see thinning.docx 
#thin split : see thinning.docx
#Thin breaks: apply thinning every X years to avoid thinning different temporal data

pkg <- c("spacesXYZ", "geosphere", "parallel", "raster", "readr", "colorspace", "rgl", "ade4", "farver")
lapply(pkg, function(x) require(x, quietly = TRUE, character.only = TRUE))

"split_cloud" <- function(cloud, threshold, option = 1){
  #--> split a cloud of points according to a threshold around median
  
  #--> cloud : coordinates of points cloud (not restricted to 2D)
  #--> threshold : threshold for splitting (maximum number of points by subclouds)
  #--> method : 1 for splitting recursively around medians, axe by axe ; 2 for splitting recursively in 2**D parts around the median centroid
  
  nc <- ncol(cloud)
  nr <- nrow(cloud)
  test <- (nr >= threshold)
  w <- 0
  subcloud <- list(cloud)
  
  if(option == 1){
    while(any(test)){
      tosplit <- subcloud[test]
      subcloud <- subcloud[!test]
      local <- list()
      d <- 1
      while((d <= nc) & (any(test))){
        for(i in 1:length(tosplit)){
          me <- median(tosplit[[i]][,d])
          t <- (tosplit[[i]][,d] > me)
          local[[i*2 - 1]] <- tosplit[[i]][t,]
          local[[i*2]] <- tosplit[[i]][!t,]
        }
        tosplit <- local	
        d <- d + 1	
        subcloud <- c(subcloud, tosplit)
        nr <- unlist(lapply(subcloud, nrow))
        test <- (nr >= threshold)
        tosplit <- subcloud[test]
        subcloud <- subcloud[!test]
      }
      subcloud <- c(subcloud, tosplit)
      nr <- unlist(lapply(subcloud, nrow))
      test <- (nr >= threshold)
      w <- w + 1
      # print(w)
    }
  }
  
  if(option == 2){
    while(any(test)){
      tosplit <- subcloud[test]
      subcloud <- subcloud[!test]
      
      for(i in 1:length(tosplit)){
        me <- apply(tosplit[[i]], 2, median)
        #--> could use means or half of ranges
        
        for(j in 1:nc){
          assign(paste0("t", j), (tosplit[[i]][,j] < me[j]))
          assign(paste0("t_", j), !(tosplit[[i]][,j] < me[j]))
          assign(paste0("fac", j), as.factor(c(paste0("t", j), paste0("t_", j))))
        }
        
        exp <- paste(paste0("fac", 1:nc), collapse = ":")
        t <- levels(eval(parse(text = exp)))
        t <- strsplit(t, ":")
        
        for(k in 1:length(t)){
          exp <- paste0("get(t[[", k, "]][", 1:nc, "])")
          exp <- paste(exp, collapse = "&")
          t[[k]] <- eval(parse(text = exp))
          t[[k]] <- tosplit[[i]][t[[k]],]
        }
        
        subcloud <- c(subcloud, t)
      }
      
      nr <- unlist(lapply(subcloud, nrow))
      test <- (nr >= threshold)
      w <- w + 1
      # print(w)
    }
  }
  
  subcloud <- subcloud[nr > 0]
  return(subcloud)
} 

"thin_split" <- function(cloud, threshold_s = Inf, fun = "dist", method = "euclidean", threshold_d, ncores = 1){
  #--> thinning (rep = 1) big clouds of points by splitting it in many subclouds according to threshold_s
  #--> allows to use other distances than the euclidean ones (actually, spacesXYZ::DeltaE for color, geosphere::distGeo for geographic distance
  
  #--> cloud : points cloud to split and thin (not restricted to 2D)
  #--> threshold_s : threshold for splitting (by default, Inf, no splitting)
  #--> fun : distance function, by default stats::dist
  #--> method : distance method associated to fun function (by default, euclidean function associated to dist but could be metric associated to spacesXYZ::DeltaE() or farver:compare_colour() [faster than DeltaE, do not need parallelization])
  #--> threshold_d : threshold distance for thinning
  #--> ncores : if you want to parrallelize your code on more cores
  
  ## Dimensions
  nc <- ncol(cloud)
  nr <- nrow(cloud)
  
  ## Memory limit size
  #--> help(Memory)
  mem_max <- Inf	#--> mem_max <- memory.size(max = NA) : deprecated
  if(!is.infinite(threshold_s)){
    mat_size <- unclass(object.size(matrix(1, threshold_s, threshold_s))) / 10**6
    if(mat_size > (mem_max / ncores)){
      print("Your threshold is too big regarding your memory, you should choose a lower one")
      break
    }	
  }	else{
    possibleError <- tryCatch(
      mat_size <- format(object.size(matrix(NA, nr, nr)), "Mb"
      ), error = function(e) e)
    if(inherits(possibleError, "error")){
      print("You should define a threshold, your cloud is too big to store a distance matrix")
    }			
  }
  
  ## Split the cloud in subclouds according to threshold
  subclouds <- split_cloud(cloud, threshold_s)
  
  ## Spatial thining on subclouds
  k <- 1
  nsub <- length(subclouds)
  if(nsub == 1){
    print("Your cloud has not been splitted")
  }	else{
    print(paste0("Your cloud has been splitted in ", nsub, " subclouds"))
  }
  res <- NULL
  
  repeat{
    # Calculate matrix of distances until memory size of submats is reached
    if(fun == "dist"){
      submats <- vector("list", nsub)
      k <- s <- nsub + 1
    }	else{
      gc()
      submats <- as.list(NULL)
      s <- 1
      mem <- unclass(object.size(submats)) / 10**6
      
      while((mem < (mem_max * 4 / 5 / ncores)) & (k <= nsub)){
        subcloud <- subclouds[[k]]
        cl <- makeCluster(ncores)
        
        if(fun == "DeltaE"){
          clusterExport(cl, varlist = c("DeltaE", "method"), envir = environment())
          submats[[s]] <- parApply(cl, subcloud, 1, function(x) DeltaE(x, subcloud, metric = method))
          
        }
        
        if(fun == "compare_colour"){
          submats[[s]] <- compare_colour(subcloud, subcloud, "luv", "luv", method = method)
          
        }
        
        if(fun == "distGeo"){
          clusterExport(cl, varlist = c("distGeo"))
          submats[[s]] <- parApply(cl, subcloud, 1, function(x) distGeo(x, subcloud))
        }
        
        stopCluster(cl)
        mem <- unclass(object.size(submats)) / 10**6
        # print(format(mem, "Mb"))
        k <- k + 1
        print(k)
        s <- s + 1
      }
    }
    
    # Execute spatial thinning
    ij <- mapply(function(x,y) c(x,y), as.list((k-(s-1)):(k-1)), as.list(1:(s-1)), SIMPLIFY = F)
    cl <- makeCluster(ncores)
    clusterExport(cl, varlist = c("thin_algo", "subclouds", "submats", "method", "threshold_d"), envir = environment())
    if(fun == "dist"){
      subres <- parLapply(cl, ij, function(x) thin_algo(subclouds[[x[1]]], submats[[x[2]]], method = method, threshold = threshold_d)[[1]])
    }	else{
      subres <- parLapply(cl, ij, function(x) thin_algo(subclouds[[x[1]]], submats[[x[2]]], threshold = threshold_d)[[1]])
    }
    stopCluster(cl)
    
    res <- rbind(res, do.call(rbind, subres))
    # save(k, res, file = "test.RData")
    
    if(k >= nsub){
      break
    }
  }
  return(res)
}

"thin_breaks"<- function(data, timestep=20 ){
  set.seed(1)
  list_index<-c()
  breaks<-seq(min(data$year), 2030, by=timestep)
  for (i in 1:length(breaks)){
    if (i < length(breaks)){
      d<- data[between (data$year, breaks[i], breaks[i+1]), ]
      system.time(
        occ_t <- thin_split(d[,c("decimallongitude", "decimallatitude")], 500, fun = "distGeo", threshold_d = 10*1000, ncores = 12)
      )
      list_index<-c(list_index, rownames(occ_t))
      
    }
  }
  return(data[list_index,])}



################################################################################
#D-Trade ####
################################################################################

"assign_trade"<-function(data, sum_import){
  country_data<-unique(data$ISO)

  for (c in country_data){
    if (c %in% sum_import$ISO3){
    data$trade[data$ISO == c & data$S =="S1"]<-sum_import$Sum[sum_import$year == 1970 & sum_import$ISO3 ==c]
    data$trade[data$ISO == c & data$S =="S2"]<-sum_import$Sum[sum_import$year == 2014 & sum_import$ISO3 ==c]
    }
    else {
      data$trade[data$ISO == c & data$S =="S1"]<-NA
      data$trade[data$ISO == c & data$S =="S2"]<-NA 
      
    }
  }
  
  return(data)
}

################################################################################
#E-Human modifications ####
################################################################################

load_HM<- function(data){
  # Convert to SpatRaster
  d_s1<-data %>% filter(S=="S1") %>% dplyr::select(decimallongitude, decimallatitude)
  d_s2<-data %>% filter(S=="S2") %>% dplyr::select(decimallongitude, decimallatitude)
  
  
  #climate_data_terra <- terra::rast(hm_1970)
  # Extract values all at once
  values_S1 <- terra::extract(hm_1970, d_s1, cells = TRUE)
  values_S2 <- terra::extract(hm_2025, d_s2, cells = TRUE)
  
  # Join results
  data_hm<-data
  data_hm$HM[data_hm$S =="S1"]<-values_S1
  data_hm$HM[data_hm$S =="S2"]<-values_S2
  
  # Remove the NAs
  #data_pca <- data_pca %>% drop_na(PCA_1, PCA_2) # from 408158 to 407807 rows
  
  return(data_hm)
}
