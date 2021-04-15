# ==========================================================================================================================================
#
# HIA, DC, lc, DOH data
# *** only need to change the date in 4 places and rerun the loops *** 
# Started: 08032020
# 
# ==========================================================================================================================================

install.packages("rJava")
install.packages("OpenStreetMap")
install.packages('extrafont')
install.packages("plyr")
remove.packages("dplyr")
install.packages("dplyr")

library(raster)
library(rgdal)
library(ggplot2)
library(ggspatial)     # no package called 'sf'
library(scales)
library(extrafont)
library(maptools)
library(extrafont)
loadfonts()
library(magrittr)
library(plyr)
library(dplyr)
library(sf)     # no package called 'sf'
library(rJava)
library(OpenStreetMap)
library(gridExtra)
library(rgeos)
library(ggpubr)
library(spatialEco)  # no package called 'spatialEco'
library(grid)

# ===============================================================================================================
# Set working directory and load files

# setwd("C:/Users/mdcastillo/Desktop/dc_hia/HIA_calculations_mc/run/") # at Milken
setwd("/Users/mdcc/Box/dc_hia/HIA_calculations_mc/run/")  # box wd

# setwd('/GWSPH/home/mdcastillo/run/')    
getwd()

# Specify where the files are:
rates <- '/Users/mdcc/Box/dc_hia/HIA_calculations_mc/run/rates_dc/'
concs <- '/Users/mdcc/Box/dc_hia/HIA_calculations_mc/run/conc_dc/vdk19/'    #  updated with concentrations from 2019 dataset
pops <- '/Users/mdcc/Box/dc_hia/HIA_calculations_mc/run/pop_dc/'
shps <- '/Users/mdcc/Box/dc_hia/HIA_calculations_mc/run/clip_dc/'

#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

theme_map <- function(...) {
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(hjust = 0, size=13), #,family="DejaVu Sans Light"),
    plot.subtitle = element_text(hjust=0, size=12), #,family="DejaVu Sans Light"),
    plot.caption = element_text(hjust=0, size=11), #,family="DejaVu Sans Light"),
    legend.title = element_text(size=12), #, family="DejaVu Sans Light"),
    legend.text = element_text(size=12), #, family="DejaVu Sans Light"),
    axis.title = element_blank(),
    legend.position = 'bottom',
    legend.justification = 'center',
    legend.spacing = unit(c(-.1,0.2,.2,0.2), "cm"),
    panel.border = element_rect(colour = "black", fill=NA, size=0.5),
    rect = element_blank())
}

myZonal <- function (x, z, stat, digits = 0, na.rm = TRUE, 
                     ...) {
  library(data.table)
  fun    <- match.fun(stat) 
  vals   <- getValues(x) 
  zones  <- round(getValues(z), digits = digits) 
  rDT    <- data.table(vals, z=zones) 
  setkey(rDT, z) 
  rDT[, lapply(.SD, fun, na.rm = TRUE), by=z] 
} 

ZonalPipe<- function (zone.in, raster.in, shp.out=NULL, stat){
  require(raster)
  require(rgdal)
  require(plyr)
  
  # Load raster
  r <- raster.in
  # Load zone shapefile
  shp <- zone.in
  # Project 'zone' shapefile into the same coordinate system than the input raster
  shp <- spTransform(shp, crs(r))
  # Add ID field to Shapefile
  shp@data$ID <- c(1:length(shp@data[,1]))
  # Crop raster to 'zone' shapefile extent
  r <- crop(r, extent(shp))	
  # Rasterize shapefile
  zone <- rasterize(shp, r, field="ID", dataType = "INT1U") # Change dataType if nrow(shp) > 255 to INT2U or INT4U
  # Zonal stats
  Zstat <- data.frame(myZonal(r, zone, stat))
  colnames(Zstat) <- c("ID", paste(names(r), "_", c(1:(length(Zstat)-1)), "_",stat))
  # Merge data in the shapefile and write it
  shp@data <- plyr::join(shp@data, Zstat, by="ID")
  
  if (is.null(shp.out)){
    return(shp)
  }else{
    writeOGR(shp, shp.out, layer= sub("^([^.]*).*", "\\1", basename(zone.in)), driver="ESRI Shapefile")
  }
}


# ///////////////////////////////////////////////////////////////////////////////////////////////////
# bring in all data inputs for the analyses

# BETAs: Lower, beta, upper
# beta.groups <- c(0.001421475836, 0.005602959022, 0.009665521211)    # Mar et al 2010 (0-99) for Asthma
# beta.groups <- c(0.0009950330853, 0.00953101798, 0.01739533071)     # Turner et al 2016 (30-99) for COPD
# beta.groups <- c(0.00198026273, 0.01310282624, 0.01988508587)       # Turner et al 2016 (30-99) for IHD
# beta.groups <- c(0.002955880224, 0.008617769624, 0.01484200051)     # Turner et al 2016 (30-99) for LC
# beta.groups <- c(0.004879016417, 0.01043600153, 0.01570037488)      # Turner et al 2016 (30-99) for Stroke
# beta.groups <- c(0.003922071315, 0.005826890812, 0.007696104114)    # Turner et al 2016 (30-99) for All-cause

# comment out pop.groups for the HIA calculation (delta mortality) when using counts instead of rates
pop.groups <- c('pop_all_div')
  names(pop.groups) <- c('GPW_All')
city.groups <- ('dc_boundary')
  names(city.groups) <- ('DC')
clip.groups <- c('png_dc') #,'ward_dc','tract_dc')
  names(clip.groups) <- c('PNGs') #,'Wards','Tracts')
beta.groups <- c(0.005826890812) 
  names(beta.groups) <- c('PE_all') 
rate.groups <- c('all_00_doh','all_01_doh','all_02_doh','all_03_doh', 'all_04_doh','all_05_doh',
                  'all_06_doh','all_07_doh','all_08_doh','all_09_doh', 'all_10_doh',
                  'all_11_doh','all_12_doh','all_13_doh','all_14_doh','all_15_doh')             
  names(rate.groups) <- c('ALL_00','ALL_01','ALL_02','ALL_03','ALL_04','ALL_05',
                         'ALL_06','ALL_07','ALL_08','ALL_09','ALL_10',
                         'ALL_11','ALL_12','ALL_13','ALL_14','ALL_15')    
conc.groups <- c('vdk_2018_dc') 
                # ('vdk_2000_dc','vdk_2001_dc','vdk_2002_dc','vdk_2003_dc','vdk_2004_dc','vdk_2005_dc',
                #  'vdk_2006_dc','vdk_2007_dc','vdk_2008_dc','vdk_2009_dc','vdk_2010_dc',
                #  'vdk_2011_dc','vdk_2012_dc','vdk_2013_dc','vdk_2014_dc','vdk_2015_dc')             
  names(conc.groups) <- c('DC_18')
                        # ('DC_00','DC_01','DC_02','DC_03','DC_04','DC_05',
                        # 'DC_06','DC_07','DC_08','DC_09','DC_10',
                        # 'DC_11','DC_12','DC_13','DC_14','DC_15')          

pdf(NULL)

getwd()
setwd("/Users/mdcc/Box/dc_hia/HIA_calculations_mc/results/all_doh/") # all_pm18/")

# ///////////////////////////////////////////////////////////////////////////////////////////////////

for (i in 1:length(beta.groups)){
  print(beta.groups[i])
  
  for (j in 1:length(conc.groups)){
  print(conc.groups[j])
    
    for (k in 1:length(rate.groups)){
      print(rate.groups[k]) 
      
      if(j == 1){
        b = raster(paste(rates, rate.groups[1], '.tif', sep=''))
        rate.names <- names(rate.groups[1])
      } 
      if (j == 2) {
        b = raster(paste(rates, rate.groups[2], '.tif', sep=''))
        rate.names <- names(rate.groups[2])
      } 
      if (j == 3) {
        b = raster(paste(rates, rate.groups[3], '.tif', sep=''))
        rate.names <- names(rate.groups[3])
      } 
      if (j == 4) {
        b = raster(paste(rates, rate.groups[4], '.tif', sep=''))
        rate.names <- names(rate.groups[4])
      } 
      if (j == 5) {
        b = raster(paste(rates, rate.groups[5], '.tif', sep=''))
        rate.names <- names(rate.groups[5])
      } 
      if (j == 6) {
        b = raster(paste(rates, rate.groups[6], '.tif', sep=''))
        rate.names <- names(rate.groups[6])
      } 
      if (j == 7) {
        b = raster(paste(rates, rate.groups[7], '.tif', sep=''))
        rate.names <- names(rate.groups[7])
      } 
      if (j == 8) {
        b = raster(paste(rates, rate.groups[8], '.tif', sep=''))
        rate.names <- names(rate.groups[8])
      } 
      if (j == 9) {
        b = raster(paste(rates, rate.groups[9], '.tif', sep=''))
        rate.names <- names(rate.groups[9])
      } 
      if (j == 10) {
        b = raster(paste(rates, rate.groups[10], '.tif', sep=''))
        rate.names <- names(rate.groups[10])
      } 
      if (j == 11) {
        b = raster(paste(rates, rate.groups[11], '.tif', sep=''))
        rate.names <- names(rate.groups[11])
      } 
      if (j == 12) {
        b = raster(paste(rates, rate.groups[12], '.tif', sep=''))
        rate.names <- names(rate.groups[12])
      } 
      if (j == 13) {
        b = raster(paste(rates, rate.groups[13], '.tif', sep=''))
        rate.names <- names(rate.groups[13])
      } 
      if (j == 14) {
        b = raster(paste(rates, rate.groups[14], '.tif', sep=''))
        rate.names <- names(rate.groups[14])
      } 
      if (j == 15) {
        b = raster(paste(rates, rate.groups[15], '.tif', sep=''))
        rate.names <- names(rate.groups[15])
      }

      for (h in 1:length(pop.groups)){
        print(pop.groups[h])
        
        for (m in 1:length(city.groups)){
          print(city.groups[m]) 
          
          a = raster(paste(concs,conc.groups[j],'.tif',sep=''))
          
          print(concs)
          print(conc.groups[j])
          
          a[a==0]<-NA
          
          af <- 1-exp(-beta.groups[i]*a)
          
          af2 <- af
          af <- af*100
          
          print(rates)
          print(rate.groups[k])
          b = raster(paste(rates,rate.groups[k],'.tif',sep=''))	
          
          # file.exists("/Users/mdcc/Box/MC_R/HIA/run/rates_dc/all_25_2011_dc.tif") # check whether file exists in folder
          
          mr <- af2*b     # betas*conc * rates
          mr[mr==0]<-NA
          
          c = raster(paste(pops,pop.groups[h],'.tif',sep=''))
          c[c==0]<- NA
          
          hia = overlay(c, b, a, fun=function(r1, r2, r3){return(r1*r2*(10^-5)*(1-exp(-beta.groups[i]*r3)))})
          hia[hia==0]<-NA
          
          require(rgdal)
          shp <- readOGR(dsn=path.expand("/Users/mdcc/Box/dc_hia/HIA_calculations_mc/run/clip_dc"), layer=paste(city.groups[m]))
          list.files("/Users/mdcc/Box/dc_hia/HIA_calculations_mc/run/clip_dc")
          
          # the data source name (dsn= argument) is the folder (directory) where the shapefile is, 
          # and the layer is the name of the shapefile (without the .shp extension)
          list.files("/Users/mdcc/Box/dc_hia/HIA_calculations_mc/run/clip_dc/", pattern='\\.shp$')
          file.exists("/Users/mdcc/Box/dc_hia/HIA_calculations_mc/run/clip_dc/png.dc.shp")
          crs(shp) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
          
          require(ggplot2)
          require(dplyr)
          shp.f <- fortify(shp) %>% 
            mutate(id = as.numeric(id))
          
          af <- crop(af, shp)
          af <- mask(af, shp)
          
          #==============================================================================================================================================================    
          # f1 = paste('/GWSPH/home/mdcastillo/results/dc/pm25/all_doh/adults/af/',names(city.groups[m]),names(conc.groups[j]),', ',names(rate.groups[k]),', ',names(pop.groups[h]),', ',names(beta.groups[i]),', ',names(city.groups[m]),'.tif',sep='')
          f1 = paste("/Users/mdcc/Box/dc_hia/HIA_calculations_mc/results/all_doh/af/",names(conc.groups[j]),', ',names(rate.groups[k]),', ',names(pop.groups[h]),', ',names(beta.groups[i]),', ',names(city.groups[m]),'.tif',sep='')
          writeRaster(af, filename=f1, format="GTiff", overwrite=TRUE)
          #==============================================================================================================================================================            
          af.iqr <- quantile(af)
          af.iqr <- as.matrix(af.iqr)
          af.iqr <- t(af.iqr)
          af.mean <- cellStats(af, 'mean')
          af.df <- cbind(af.mean, af.iqr)
          
          print(paste(names(city.groups[m]),', attributable fraction, ',names(conc.groups[j]),', ',names(beta.groups[i]),sep=''))
          print(af.df)
          
          min.af <- minValue(af)
          max.af <- maxValue(af)
          min.af.label <- round(minValue(af),2)
          max.af.label <- round(maxValue(af),2)
          mean.af <- (min.af+max.af)/2
          mean.af.label <- round(mean.af,2)
          
          af.df <- rasterToPoints(af)
          af.df <- data.frame(af.df)
          colnames(af.df) <- c('lon','lat','val')
          
          library(extrafont)
          loadfonts()
          loadfonts(device="postscript")
          
          af_plot <- ggplot()  +
            geom_polygon(data = shp.f, aes(x = long, y = lat, group = group), 
                         fill="grey50",alpha=0.5)+
            geom_tile(data=af.df,aes(lon, lat, fill = val),alpha=0.8) +
            scale_fill_gradient2("Attributable Fraction (%)",
                                 low = "#3ec267", 
                                 mid = "#fff429",  #ff7e29
                                 high = "#fc0339", ##ff1f40
                                 midpoint = mean.af,
                                 breaks=c(min.af,mean.af,max.af),
                                 labels=c(min.af.label,mean.af.label,max.af.label),
                                 limits=c(min.af, max.af),
                                 na.value = 'grey50',
                                 guide = guide_colourbar(
                                   direction = "horizontal",
                                   label=TRUE,
                                   keyheight = unit(2, units = "mm"),
                                   title.position = 'top',
                                   title.hjust = 0.5,
                                   label.hjust = 0.5,
                                   barwidth = 15,
                                   nrow = 1,
                                   byrow = T,
                                   label.position = "bottom"))+
            theme_map()+
            geom_path(data = shp.f, aes(x = long, y = lat, group = group), 
                      color = "grey60", size = 0.5)+
            labs(title='Attributable Fraction',
                 subtitle=paste0('Range: ',min.af.label,' to ',max.af.label,'%.'),sep='')
          ggsave(paste(names(city.groups[m]),' AF All-Cause and PM25 ',names(beta.groups[i]),' ',names(conc.groups[j]),' ',names(pop.groups[h]),' ',names(rate.groups[k]),'.png',sep=''),dpi=300)
          print('af')
          
          
          #==============================================================================================================================================================    
          
          # Create mortality risk	
          mr <- crop(mr, shp)
          mr <- mask(mr, shp)
          
          # Cell stats sink to csv
          mr.iqr <- quantile(mr)
          mr.iqr <- as.matrix(mr.iqr)
          mr.iqr <- t(mr.iqr)
          mr.mean <- cellStats(mr, 'mean')
          df.mr <- cbind(mr.mean, mr.iqr)
          print(paste(names(city.groups[m]),' Mortality Risk ',names(conc.groups[j]),' ',names(beta.groups[i]),' ',names(rate.groups[k]),sep=''))
          print(df.mr)
          
          f2 = paste("/Users/mdcc/Box/dc_hia/HIA_calculations_mc/results/all_doh/mr/",names(conc.groups[j]),', ',names(rate.groups[k]),', ',names(pop.groups[h]),', ',names(beta.groups[i]),', ',names(city.groups[m]),'.tif',sep='')
          writeRaster(mr, filename=f2, format="GTiff", overwrite=TRUE)
          
          # Create maps mortality risk 
          min.mr <- minValue(mr)
          min.mr.label <- round(minValue(mr),2)
          max.mr <- maxValue(mr)
          max.mr.label <- round(maxValue(mr),2)
          mean.mr <- (min.mr+max.mr)/2
          mean.mr.label <- round(mean.mr,2)
          
          mr.df <- rasterToPoints(mr)
          mr.df <- data.frame(mr.df)
          colnames(mr.df) <- c('lon','lat','val')
          
          # MAP MORTALITY RISK
          mr_plot <- ggplot()  +
            geom_polygon(data = shp.f, aes(x = long, y = lat, group = group), 
                         fill="grey50",alpha=0.5)+
            geom_tile(data=mr.df,aes(lon, lat, fill = val),alpha=0.8) +
            scale_fill_gradient2("Mortality Risk per 10,000",
                                 low = "#3ec267", 
                                 mid = "#fff429",  #ff7e29
                                 high = "#fc0339", ##ff1f40
                                 midpoint = mean.mr,
                                 breaks=c(min.mr,mean.mr,max.mr),
                                 labels=c(min.mr.label,mean.mr.label,max.mr.label),
                                 limits=c(min.mr, max.mr),
                                 na.value = 'grey50',
                                 guide = guide_colourbar(
                                   direction = "horizontal",
                                   label=TRUE,
                                   keyheight = unit(2, units = "mm"),
                                   title.position = 'top',
                                   title.hjust = 0.5,
                                   label.hjust = 0.5,
                                   barwidth = 15,
                                   nrow = 1,
                                   byrow = T,
                                   label.position = "bottom"))+
            theme_map()+
            geom_path(data = shp.f, aes(x = long, y = lat, group = group), 
                      color = "grey60", size = 0.5)+
            labs(title='Mortality Risk',
                 caption=paste0(names(beta.groups[i]),', ',names(conc.groups[j]),', ',names(rate.groups[k]),'.',sep=''),
                 subtitle=paste0('Range: ',min.mr.label,' to ',max.mr.label,' per 10,000. '),sep='')
          ggsave(paste(names(city.groups[m]),' MR All-Cause and PM25 ',names(beta.groups[i]),' ',names(conc.groups[j]),' ',names(pop.groups[h]),' ',names(rate.groups[k]),'.png',sep=''),dpi=300)
          print('mr')
          
          #===========================================================================================================================      
          # Crop and map HIA files
          
          hia <- crop(hia, shp)
          hia <- mask(hia, shp)
          
          # f3 = paste('/GWSPH/home/mdcastillo/results/dc/pm25/all_doh/adults/paf/',names(conc.groups[j]),', ',names(rate.groups[k]),', ',names(pop.groups[h]),', ',names(beta.groups[i]),', ',names(city.groups[m]),'.tif',sep='')
          f3 = paste("/Users/mdcc/Box/dc_hia/HIA_calculations_mc/results/all_doh/paf/",names(conc.groups[j]),', ',names(rate.groups[k]),', ',names(pop.groups[h]),', ',names(beta.groups[i]),', ',names(city.groups[m]),'.tif',sep='')
          writeRaster(hia, filename=f3, format="GTiff", overwrite=TRUE)
          
          #Cell stats sink to csv
          hia.iqr <- quantile(hia)
          hia.iqr <- as.matrix(hia.iqr)
          hia.iqr <- t(hia.iqr)
          hia.mean <- cellStats(hia, 'mean')
          hia.sum <- cellStats(hia,'sum')
          df.hia <- cbind(hia.sum,hia.mean, hia.iqr)
          print(paste(names(city.groups[m]),' excess mortality ',names(conc.groups[j]),' ',names(beta.groups[i]),' ',names(pop.groups[h]),' ',names(rate.groups[k]),sep=''))
          print(df.hia)
          
          hia.df <- rasterToPoints(hia)
          hia.df <- data.frame(hia.df)
          colnames(hia.df) <- c('lon','lat','val')
          
          # ----------------------------------------------
          # CREATE PAF MAPS 
          min.hia <- minValue(hia)
          max.hia <- maxValue(hia)
          min.hia.label <- round(minValue(hia),2)
          max.hia.label <- round(maxValue(hia),2)
          mean.hia <- (min.hia+max.hia)/2
          mean.hia.label <- round(mean.hia,2)
          
          paf_plot <- ggplot()  +
            geom_polygon(data = shp.f, aes(x = long, y = lat, group = group), 
                         fill="grey50",alpha=1)+
            geom_tile(data=hia.df,aes(lon, lat, fill = val),alpha=0.8) +
            scale_fill_gradient2("Excess cases (n)",
                                 low = "#3ec267", 
                                 mid = "#fff429",  #ff7e29
                                 high = "#fc0339", ##ff1f40
                                 midpoint = mean.hia,
                                 breaks=c(min.hia,mean.hia,max.hia),
                                 labels=c(min.hia.label,mean.hia.label,max.hia.label),
                                 limits=c(min.hia, max.hia),
                                 na.value = 'grey50',
                                 guide = guide_colourbar(
                                   direction = "horizontal",
                                   label=TRUE,
                                   keyheight = unit(2, units = "mm"),
                                   title.position = 'top',
                                   title.hjust = 0.5,
                                   label.hjust = 0.5,
                                   barwidth = 15,
                                   nrow = 1,
                                   byrow = T,
                                   label.position = "bottom"))+
            theme_map()+
            geom_path(data = shp.f, aes(x = long, y = lat, group = group), 
                      color = "grey60", size = 0.5)+
            labs(title=bquote('All-Cause Excess Mortality Attributable to '~PM[2.5]~''),
                 caption=paste0(names(beta.groups[i]),', ',names(conc.groups[j]),', ',names(rate.groups[k]),',\n',names(pop.groups[h]),'.',sep=''),
                 subtitle=paste0('Range: ',min.hia.label,' to ',max.hia.label,' per grid cell. '),sep='', size = 21)
          ggsave(paste(names(city.groups[m]),' PAF All-Cause and PM25 excess cases, ',names(beta.groups[i]),' ',names(conc.groups[j]),' ',names(pop.groups[h]),' ',names(rate.groups[k]),'.png',sep=''),dpi=300)
          print('hia')
          
          
          #=========================================================================================================================
          
          clip.groups <- c('png_dc','tract_dc','ward_dc')
          names(clip.groups) <- c('PNGs','Tracts','Wards')
          
          for (s in 1:length(clip.groups)){
            print(clip.groups[s])
            
            clips.shp <- readOGR(dsn=path.expand("/Users/mdcc/Box/dc_hia/HIA_calculations_mc/run/clip_dc"), layer=paste(clip.groups[s]))
            list.files("/Users/mdcc/Box/dc_hia/HIA_calculations_mc/run/clip_dc") # list all files in folder
            list.files("/Users/mdcc/Box/dc_hia/HIA_calculations_mc/run/clip_dc/", pattern='\\.shp$') # list all shapefiles in folder
            file.exists("/Users/mdcc/Box/dc_hia/HIA_calculations_mc/run/clip_dc/png_dc.shp")
            
            crs(clips.shp) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
            
            requireNamespace("rgeos") 
            
            clips.shp <- crop(clips.shp, shp)    # crop all clips (pngs, tracts and wards) to the extent of shp (boundary_dc)
            
            clips.shp.f <- fortify(clips.shp) %>%
              mutate(id = as.numeric(id))
            
            zone.in <- clips.shp
            raster.in <- hia
            
            require(rgdal)
            shp2 <- ZonalPipe(zone.in, raster.in, stat="sum")
            shp2@data <- shp2@data %>% mutate(id = row.names(.))
            
            shp_df <- fortify(shp2, region = "id")
            shp_df <- shp_df %>% left_join(shp2@data, by = c("id"="id"))
            shp_df <- as.data.frame(shp_df)
            shp_df[,ncol(shp_df)][shp_df[,ncol(shp_df)] == 0] <- NA
            r.min <- min(shp_df[,ncol(shp_df)],na.rm=TRUE)
            r.max <- max(shp_df[,ncol(shp_df)],na.rm=TRUE)
            r.med <- median(shp_df[,ncol(shp_df)],na.rm=TRUE)
            colnames(shp_df)[ncol(shp_df)] <- 'hia.val'
            
            r.mean <- (r.min+r.max)/2
            r.mean.label <- round(r.mean,2)
            r.min.label <- round(r.min,2)
            r.max.label <- round(r.max,2)
            r.med.label <- round(r.med,2)
            
            zone.in <- clips.shp
            raster.in <- c
            
            shp3 <- ZonalPipe(zone.in, raster.in, stat="sum")
            shp3@data <- shp3@data %>% mutate(id = row.names(.))
            pop_df <- fortify(shp3, region = "id")
            pop_df <- pop_df %>% left_join(shp3@data, by = c("id"="id"))
            pop_df <- as.data.frame(pop_df)
            pop_df[,ncol(pop_df)][pop_df[,ncol(pop_df)] == 0] <- NA
            colnames(pop_df)[ncol(pop_df)] <- "pop.val"
            
            rate_df <- merge(shp_df,pop_df,by='order')
            rate_df <- as.data.frame(rate_df)
            rate_df$rate <- NA
            rate_df$rate <- (rate_df$hia.val*100000)/rate_df$pop.val
            #rate_df$rate[rate_df$rate==0]<-NA
            rate.min <- min(rate_df[,ncol(rate_df)],na.rm=TRUE)
            rate.max <- max(rate_df[,ncol(rate_df)],na.rm=TRUE)
            rate.med <- median(rate_df[,ncol(rate_df)],na.rm=TRUE)
            rate.min.label <- round(rate.min,2)
            rate.max.label <- round(rate.max,2)
            rate.med.label <- round(rate.med,2)
            rate.mean <- (rate.min+rate.max)/2
            rate.mean.label <- round(rate.mean,2)
            
            
            # write results to csv and raster 
            write.csv(rate_df, paste(names(clips.shp[s]),' ',names(beta.groups[i]),' ',names(conc.groups[j]),' ',names(pop.groups[h]),' ',names(rate.groups[k]),'results.csv'))
            
            # MAP OF EXCESS MORTALITY PER PNG
            e <- ggplot()  +
              geom_polygon(data = shp_df, aes(x = long, y = lat, group = group, fill = shp_df[,ncol(shp_df)]),alpha=1)+
              scale_fill_gradient2("Count Cases (n)",
                                   low = "#3ec267",
                                   mid = "#fff429",  #ff7e29
                                   high = "#fc0339", ##ff1f40
                                   midpoint = r.mean,
                                   na.value='grey50',
                                   breaks=c(r.min,r.mean,r.max),
                                   labels=c(r.min.label,r.mean.label,r.max.label),
                                   limits=c(r.min, r.max),
                                   guide = guide_colourbar(
                                     direction = "horizontal",
                                     label=TRUE,
                                     keyheight = unit(2, units = "mm"),
                                     title.position = 'top',
                                     title.hjust = 0.5,
                                     label.hjust = 0.5,
                                     barwidth = 15,
                                     nrow = 1,
                                     byrow = T,
                                     label.position = "bottom"))+
              theme_map()+
              geom_path(data = clips.shp.f, aes(x = long, y = lat, group = group),
                        color = "grey60", size = 0.3)+
              labs(title=bquote('All-Cause Excess Mortality Attributable to '~PM[2.5]~''),
                   caption=paste0(names(beta.groups[i]),', ',names(conc.groups[j]),', \n',names(rate.groups[k]),', ',names(pop.groups[h]),'.',sep=''),
                   subtitle=paste0('Range: ',r.min.label,' to ',r.max.label),sep='')
            ggsave(paste(names(clips.shp[s]),' PAF counts, ',names(beta.groups[i]),' ',names(conc.groups[j]),' ',names(pop.groups[h]),' ',names(rate.groups[k]),'.png',sep=''),dpi=300)
            print('clips.shp')
            
            
            e_rate <- ggplot()  +
              geom_polygon(data = rate_df, aes(x = long.x, y = lat.x, group = group.x, fill = rate_df$rate),alpha=1) +
              scale_fill_gradient2("Rate per 100,000",
                                   low = "#3ec267",
                                   mid = "#fff429",  #ff7e29
                                   high = "#fc0339", ##ff1f40
                                   midpoint = rate.mean,
                                   na.value ='grey50',
                                   breaks = c(rate.min,rate.mean,rate.max),
                                   labels = c(rate.min.label,rate.mean.label,rate.max.label),
                                   limits = c(rate.min, rate.max),
                                   guide = guide_colourbar(
                                     direction = "horizontal",
                                     label=TRUE,
                                     keyheight = unit(2, units = "mm"),
                                     title.position = 'top',
                                     title.hjust = 0.5,
                                     label.hjust = 0.5,
                                     barwidth = 15,
                                     nrow = 1,
                                     byrow = T,
                                     label.position = "bottom"))+
              theme_map() +
              geom_path(data = clips.shp.f, aes(x = long, y = lat, group = group), # add clip to map
                        color = "grey60", size = 0.3) +
              labs(title = bquote('All-Cause Excess Mortality Attributable to '~PM[2.5]~' (2011)', plot.title = element_text (hjust = 0.5)),
                   caption = paste0(names(beta.groups[i]),', ',names(conc.groups[j]),', \n',names(rate.groups[k]),', ',names(pop.groups[h]),'.',sep=''),
                   subtitle = paste0('Range: ',rate.min.label,' to ',rate.max.label,' per 100,000. '),sep='')
            ggsave(paste(names(clips.shp[s]),' PAF rates, ',names(beta.groups[i]),' ',names(conc.groups[j]),' ',names(pop.groups[h]),' ',names(rate.groups[k]),'.png',sep=''),dpi=300)
            print('rate.png')
            
            
          
            # ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            
            # composite 2 with color scale 
            clips.png <- readOGR(dsn=path.expand("/Users/mdcc/Box/dc_hia/HIA_calculations_mc/run/clip_dc"), layer="png_dc")
            crs(clips.png) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
            requireNamespace("rgeos")
            clips.png <- crop(clips.png, shp)    # crop all clips (pngs, tracts and wards) to the extent of shp (boundary_dc)
            
            require(ggplot2)
            e_rate3 <- ggplot() +
              geom_polygon(data = rate_df, aes(x = long.x, y = lat.x, group = group.x, fill = rate_df$rate), alpha=4)+
              scale_fill_gradient2( "Rate per 100,000", 
                                    # low = "#3ec267",
                                    low = "#3ec267",   #ff7e29    
                                    mid = "#fff429",
                                    high = "#fc0339", ##ff1f40
                                    midpoint = rate.mean,
                                    na.value='grey50',
                                    guide = guide_colourbar(
                                      direction = "horizontal",
                                      label=TRUE,
                                      keyheight = unit(2, units = "mm"),
                                      title.position = 'top',
                                      title.hjust = 0.5,
                                      label.hjust = 0.5,
                                      barwidth = 15,
                                      nrow = 1,
                                      byrow = T,
                                      label.position = "bottom")) +
              theme_map() + theme(panel.border = element_blank()) +    # override rectangular box for map to none
              geom_path(data = clips.shp.f, aes(x = long, y = lat, group = group),
                        color = "grey60", size = 0.5) +
              labs(title=paste0('All (2015)'),
                   subtitle=paste0('Range: ',rate.min.label,' - ',rate.max.label,' per 100,000'),sep='') + 
              theme (plot.title = element_text(face="bold", hjust = 0.5, size = 23),
                     plot.subtitle = element_text(face = "italic", hjust = 0.5, size = 21),
                     legend.title = element_text(size=18), legend.text = element_text(size=18))
            print(e_rate3)
            ggsave(paste('All_comp_color ',names(clips.shp[s]),', ',names(conc.groups[j]),', ',names(beta.groups[i]),', ',names(rate.groups[k]),'.png',sep=''),
                   width = 8, height = 10, dpi = 300, units = "in")
            print("ALL")
            
            
            
            # create a plot for PNGs and save it to construct composite figure
            clips.png <- readOGR(dsn=path.expand("/Users/mdcc/Box/dc_hia/HIA_calculations_mc/run/clip_dc"), layer="png_dc")
            crs(clips.png) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
            requireNamespace("rgeos")
            clips.png <- crop(clips.png, shp)    # crop all clips (pngs, tracts and wards) to the extent of shp (boundary_dc)
            
            require(ggplot2)
            e_rate2 <- ggplot()  +
              geom_polygon(data = rate_df, aes(x = long.x, y = lat.x, group = group.x, fill = rate_df$rate), alpha=4)+
              scale_fill_gradient2( "Rate per 100,000",
                                    # low = "#3ec267",
                                    low = "slategray1",   #ff7e29    #  3ec267 is green, fff429 is yellow, fc0339 is red
                                    mid = "slategray3",
                                    high = "slategray4", ##ff1f40
                                    midpoint = rate.mean,
                                    na.value='grey50',
                                    guide = FALSE
              )+
              theme_map() + theme(panel.border = element_blank()) +    # override rectangular box for map to none
              geom_path(data = clips.shp.f, aes(x = long, y = lat, group = group),
                        color = "grey60", size = 0.5) +
              # labs(title=bquote('ALL')) #,
              labs(title=paste0('ALL (2011)'),
                   subtitle=paste0('Range: ',rate.min.label,' - ',rate.max.label,' per 100,000'),sep='') + 
              theme (plot.title = element_text(face="bold", hjust = 0.5, size = 23),
                     plot.subtitle = element_text(face = "italic", hjust = 0.5, size = 21))
            print(e_rate2)
            ggsave(paste('ALL_composite_',names(clips.shp[s]),names(beta.groups[i]),'_',names(rate.groups[k]),'.png',sep=''),dpi=300)
            print("lc")
          }
          
          rm(clip.groups)
          
        }
        rm(shp)
        rm(c)
        rm(hia)
        rm(b)
        rm(mr)
        rm(af)
        rm(af2)
        rm(a)
      }
      
    }
    
  }
}


