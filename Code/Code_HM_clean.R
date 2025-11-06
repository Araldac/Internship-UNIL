##hm code cleaned with the same way of measuring anthromes and pop density as GSC human modif

#download/install packages
pkg <- c( "countrycode", "zoo", "dplyr", "tidyr", "purrr", "pbapply", "terra","raster", "popgrids", "sp", "scales", "RColorBrewer")
lapply(pkg, function(x) require(x, quietly = TRUE, character.only = TRUE))



#load data
d_pop1970<-raster("Internship UNIL/PREV 05-11/D-Human_modif/popd_1970AD.asc")
d_pop2025<-raster("Internship UNIL/PREV 05-11/D-Human_modif/popd_2025AD.asc")

d_ant1970<-raster("Internship UNIL/PREV 05-11/D-Human_modif/anthromes1970AD.asc")
d_ant2025<-raster("Internship UNIL/PREV 05-11/D-Human_modif/anthromes2025AD.asc")

g_1980<-subset(gdp_grid, match("g_1980", names(gdp_grid)))
g3_2020<-subset(gdp_grid, match("g3_2020", names(gdp_grid)))


##change anthromes to specific value
anth_values<-c("11", "12", "21", "22", "23", "24", "31", "32", "33", "34", "41", "42", "43", "51", "52", "53", "54", "61", "62", "63", "70")
anth_values<-as.numeric(anth_values)
lof_values<-c(rep(10,2), 7, 10, 7, 7, 7, 7, 4, 4, 4, 4, 0, 4, 4, 0, 2, 0, 0, 0, 0)
anth_to_lw<-data.frame(anth_values, lof_values)

#substitue values 
landu_1970<-subs(d_ant1970, anth_to_lw)
landu_2025<-subs(d_ant2025, anth_to_lw)


##rearragne pop density
"rearrange_popdensity"<-function(raster_data){
  values<-raster::getValues(raster_data)
  # Identify the 99.9th percentile threshold
  cutoff <- quantile(values, 0.999, na.rm = TRUE)
  
  # Cap values above the 99.9th percentile
  gpw_values_capped <- ifelse(values > cutoff, cutoff, values)

  #Log-transform the capped values using log(X + 1)
  gpw_log_transformed <- log1p(gpw_values_capped)

  #Max-normalize the log-transformed values
  gpw_normalized <- gpw_log_transformed / max(gpw_log_transformed, na.rm = TRUE)

  #Create a new raster with the normalized values
  gpw_raster_normalized <- terra::setValues(raster_data, gpw_normalized)
  
  return(gpw_raster_normalized)

}

upd_pop1970<-rearrange_popdensity(d_pop1970)
upd_pop2025<-rearrange_popdensity(d_pop2025)

##change resolution of GDP raster to match the other two
fin_g_2020<-resample(g3_2020, upd_pop2025)
fin_g_1980<-resample(g_1980,upd_pop2025 )

#transformt them to values between 0-10
fin_g_2020 = calc(fin_g_2020, fun=function(x, max=g3_2020@data@max) { (x/max)*10 })
fin_g_1980 = calc(fin_g_1980, fun=function(x, max=g_1980@data@max) { (x/max)*10 })

##Sum all the three values
hm_1970<- overlay(fin_g_1980, upd_pop1970, landu_1970, fun=function(x,y,z)x+y+z)
hm_2025<- overlay(fin_g_2020, upd_pop2025, landu_2025, fun=function(x,y,z)x+y+z)

##plot the two maps
my.palette <- rev(brewer.pal(n = 30, name = "Spectral"))
par(mfrow=c(1,2))
plot(hm_1970, col=my.palette, main="Human footprint 1970")
plot(hm_2025, col=my.palette, main="Human footprint 2025")

save(hm_1970, file="Internship UNIL/Data/global_rasters/hm_1970.RData")
save(hm_2025, file="Internship UNIL/Data/global_rasters/hm_2025.RData")




