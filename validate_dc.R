# code to validate results from HIA from Veronica
# adapted to DC PM25 02252020
# example: DC 

library(raster)
library(rgdal)

# setwd('/home/vtinney/run/results/no2/all.cause/paf/bay/')
setwd("C:/Users/mdcastillo/Desktop/MC_R/HIA/results/dc/pm25/ihd/adults/paf/")
getwd()

# Select results to validate. In this case: All ages, png code rate, larkin concentrations and atkinson beta
test <- raster("PM25 Global Estimates, BDR, IHD, Adults, 2010, GPWv4, Lepuele 2012, Point Estimate, PNG clip.tif")

# Read in all input files
setwd("C:/Users/mdcastillo/Desktop/MC_R/HIA/run/")
pop_gpw <- raster('pop_dc/gpw_2010_dc.tif')
rates_png <- raster('rates_dc/mortality.ihd_25_2010_dc.tif')
conc_vdk <- raster('conc_dc/conc.pm25_2010_dc.tif')
popfrac <- raster('pop/popfrac.bay.day.tif')
beta <- c(0.023111172,0.013102826,0.033647224)

# beta <- raster('betas/beta.atkin.ug.tif')
# beta.groups <- c(0.023111172,0.013102826,0.033647224) # Lepuele 2012 (25-99)
# names(beta.groups) <- c('Lepuele 2012, Point Estimate','Lepuele 2012, Lower CI','Lepuele 2012, Upper CI')

# Read in rate and weighted rate results to test as well
setwd("C:/Users/mdcastillo/Desktop/MC_R/HIA/results/dc/pm25/ihd/adults/")
rate <- raster('rate.atkin.ug.png.all.bay.day.larkin.ug.tif')
wrate <- raster('wrate.atkin.ug.png.all.bay.day.larkin.ug.tif')


#Random sample of the results. Randomly choose a grid cell > 0 in order to test results for HIA and rates.
z <- sampleRandom(test, size=10, na.rm=TRUE, xy=TRUE, rowcol=TRUE)
# > z
# row col         x        y PM25_Global_Estimates._BDR._IHD._Adults._2010._GPWv4._Lepuele_2012._Point_Estimate._PNG_clip
# [1,]  10  18 -76.94426 38.89891                                                                                   11.4806938
# [2,]   7   2 -77.10475 38.92956                                                                                    3.6199024
# [3,]   5   6 -77.06463 38.95000                                                                                    8.8854504
# [4,]  18  10 -77.02451 38.81718                                                                                    0.7080541
# [5,]   9  17 -76.95429 38.90913                                                                                    3.4813495
# [6,]   7   9 -77.03454 38.92956                                                                                   55.1798935
# [7,]   3   7 -77.05460 38.97043                                                                                    5.4230652
# [8,]   9   7 -77.05460 38.90913                                                                                    8.1407413
# [9,]  12  18 -76.94426 38.87848                                                                                   11.1683207
# [10,]   9   9 -77.03454 38.90913                                                                                   21.3333912

# Test - use the row and column from results to test
# [6,]   7   9 -77.03454 38.92956                                                                                   55.1798935


# ***************************************************************************************
# z
# row  col         x        y bay.atkin.ug.png.all.bay.day.larkin.ug
# [1,]  812 1408 -122.4596 38.18875                              0.0000000
# [2,]  288 1016 -122.7862 38.62542                              0.0000000
# [3,] 2107 2129 -121.8587 37.10958                              0.0000000
# [4,]  971 2407 -121.6271 38.05625                              0.0000000
# [5,] 1886 2565 -121.4954 37.29375                              0.0000000
# [6,] 1097 1917 -122.0354 37.95125                              0.0262124
# [7,] 1154 2318 -121.7012 37.90375                              0.0000000
# [8,] 1712 2245 -121.7621 37.43875                              0.0000000
# [9,] 1554 2403 -121.6304 37.57042                              0.0000000
#[10,]  396  602 -123.1313 38.53542                              0.0000000


# Test - use the row and column from results to test
# [6,] 1097 1917 -122.0354 37.95125                              0.0262124
# ***************************************************************************************


# CREATE RASTER STACK
r <- stack(pop, png, beta, larkin, popfrac, wrate, rate)
r <- stack(pop_gpw, rates_png, beta, conc_vdk)


# GET RATE VALUES FROM THE SAME GRIDCELL 
wrate[1097, 1917] 
  #0.000244152
rate[1097, 1917]
  #45.98666


# SUBSET BASED ON THE ABOVE RANDOM SAMPLE TO GET THE INDIVIDUAL RASTER LAYERS 
r[1097, 1917] # row, column of the grid
    #     pop.bay.day mortality.png.all beta.atkin.ug conc.larkin.ug popfrac.bay.day
    #[1,]          57          86.39152   0.002078254          26.32    5.309192e-06

r[7,9]
    #     gpw_2010_dc mortality.ihd_25_2010_dc conc.pm25_2010_dc
    # [1,]    9245.815                   242.05             12.25



# MANUAL HIA CALCULATION USING GRID CELL VALUES (CAN DID THIS IN EXCEL TOO) 
t <- 57*86.39152*(10^-4)*(1-exp(-0.002078254*26.32))
    # [1] 0.0262124
    # t - matches results from the loop 

t <- 9245.815*242.05*(10^-4)*(1-exp(-0.023111172*12.25))
    # [1] 55.17989
 


# TEST RATE
q <- (0.0262124*100000)/57
q
    #[1] 45.98667
    # matches results from the loop - rate per 100,000

q <- (55.1798935*100000)/9245.815
q
    # [1] 596.8094


#  TEST RATE WITH POPULATION FRACTION
p <- (0.0262124*100000)/57*5.309192e-06
p
    # [1] 0.000244152
    # matches results from the loop - rate weighted by population fraction (pop in cell/pop in total)


