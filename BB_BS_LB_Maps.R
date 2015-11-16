##############################################################################
# title         : BB_BS_LB_Maps.R;
# purpose       : Generate .png files for display of predicted disease severity;
# producer      : prepared by A. Sparks;
# last update   : in Los Baños, Laguna, November 2015;
# inputs        : EPIRICE output from 2001-2008 for BB, BS and LB;
# outputs       : maps of BB, BS and LB for BGD, IND and NPL ;
# remarks 1     :;
# Licence:      : GPL2;
##############################################################################

#### Load libraries #####
library(raster)
library(rgdal)
library(maptools)
library(gpclib)
library(ggplot2)
library(plyr)
library(RColorBrewer)
library(cartography)
library(classInt)
#### End load libraries ####

#### Load data ####
# GAUL Level 1 country layer (FAO)
IND <- readOGR(dsn = "/Users/asparks/Google Drive/Data/gaul/g2015_2014_1/India",
                layer = "India")

# EPIRICE Output
diseases <- list(stack(list.files(path = "Data/EPIRICE 25deg 01-08 PK1/",
                                  pattern = "bblight_audpc.tif$", full.names = TRUE)),
                 stack(list.files(path = "Data/EPIRICE 25deg 01-08 PK1/",
                                  pattern = "bspot_audpc.tif$", full.names = TRUE)),
                 stack(list.files(path = "Data/EPIRICE 25deg 01-08 PK1/",
                                  pattern = "_blast_audpc.tif$", full.names = TRUE)))

names(diseases) <- c("bb", "bs", "lb")

for (i in 1:3) {
  # extract the values for each state
  j <- extract(mean(diseases[[i]], na.rm = TRUE), IND, coordinates(IND),
               method = "bilinear", fun = mean)

  # unlist and generate mean values for each polygon
  j <- data.frame(unlist(lapply(j, FUN = mean)))

  row.names(j) <- row.names(IND)
  names(j) <- names(diseases[i])

  row.names(IND) <- row.names(IND)

  assign(paste("IND", names(diseases[i]), sep = "."),
         spCbind(IND, j))
}

rm("i", "j", "diseases")

IND.bb.df <- fortify(IND.bb, id = "bb", region = "bb")
IND.bs.df <- fortify(IND.bs, id = "bs", region = "bs")
IND.lb.df <- fortify(IND.lb, id = "lb", region = "lb")
names(IND.bb.df) <- names(IND.bs.df) <- names(IND.lb.df) <- c("Longitude", "Latitude", "order", "hole", "piece", "group", "id")

IND.bb.breaks <- round(classIntervals(as.numeric(IND.bb.df$id), 5,
                                      style = "equal", labels = FALSE)$brks, 0)
IND.bs.breaks <- round(classIntervals(as.numeric(IND.bs.df$id), 5,
                                      style = "equal", labels = FALSE)$brks, 0)
IND.lb.breaks <- round(classIntervals(as.numeric(IND.lb.df$id), 5,
                                      style = "equal", labels = FALSE)$brks, 0)

labs <- c("Low", "Moderately Low", "Moderate", "Moderately Severe", "Severe")

IND.bb@data$bb <- cut(IND.bb@data$bb, breaks = IND.bb.breaks,
                      include.lowest = TRUE,
                      labels = labs)
IND.bs@data$bs <- cut(IND.bs@data$bs, breaks = IND.bs.breaks,
                      include.lowest = TRUE,
                      labels = labs)
IND.lb@data$lb <- cut(IND.lb@data$lb, breaks = IND.lb.breaks,
                      include.lowest = TRUE,
                      labels = labs)

writeOGR(IND.bb, dsn = "Cache", layer = "India_Bacterial_Blight_2001-2008",
         driver = "ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(IND.bs, dsn = "Cache", layer = "India_Brown_Spot_2001-2008",
         driver = "ESRI Shapefile", overwrite_layer = TRUE)
writeOGR(IND.lb, dsn = "Cache", layer = "India_Leaf_Blast_2001-2008",
        driver = "ESRI Shapefile", overwrite_layer = TRUE)
#eos
