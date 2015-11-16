##############################################################################
# title         : BB_BS_LB_Maps.R;
# purpose       : Generate .png files for display of predicted disease severity;
# producer      : prepared by A. Sparks;
# last update   : in Los Baños, Laguna, September 2015;
# inputs        : EPIRICE output from 2001-2008 for BB, BS and LB;
# outputs       : maps of BB, BS and LB for BGD, IND and NPL ;
# remarks 1     : requires BB_BS_LB_csv_for_GIS.R to have been run and generated
#               : obligatory .csv files;
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
#### End load libraries ####

#### Load data ####
# GAUL Level 1 country layer (FAO)
gaul <- readOGR(dsn = "/Users/asparks/Google Drive/Data/gaul/g2015_2014_1/India",
                layer = "India")

# EPIRICE Output
bb <- stack(list.files(path = "Data/EPIRICE 25deg 01-08 PK1/",
                       pattern = "bblight_audpc.tif$",
                       full.names = TRUE))
bs <- stack(list.files(path = "Data/EPIRICE 25deg 01-08 PK1/",
                       pattern = "bspot_audpc.tif$",
                       full.names = TRUE))
lb <- stack(list.files(path = "Data/EPIRICE 25deg 01-08 PK1/",
                       pattern = "_blast_audpc.tif$",
                       full.names = TRUE))

#### Extract values for disease by state level

bb <- extract(bb, gaul, method = "bilinear", df = TRUE)
bs <- data.frame(values(bs))
lb <- data.frame(values(lb))

IND.BB$BB <- factor(IND.BB$BB, levels(IND.BB$BB)[c(1, 3, 2, 4, 5)])
IND.BS$BS <- factor(IND.BS$BS, levels(IND.BS$BS)[c(1, 3, 2, 4, 5)])
IND.LB$LB <- factor(IND.LB$LB, levels(IND.LB$LB)[c(1, 3, 2, 4, 5)])

IND <- SA[SA@data$ADM0_NAME == "India", ]

IND@data$id <- rownames(IND@data)

IND.df <- fortify(IND, region = "id")

IND.df <- join(IND.df, IND@data, by = "id")

IND.BB.df <- join(IND.df, IND.BB, by = "ADM1_CODE")
IND.BS.df <- join(IND.df, IND.BS, by = "ADM1_CODE")
IND.LB.df <- join(IND.df, IND.LB, by = "ADM1_CODE")

# IND
# BB
ggplot(data = IND.BB.df, aes(long, lat, group = group, fill = BB)) +
  geom_polygon(color = "white", size = 0.2) +
  scale_fill_brewer(palette = "GnBu", name = "Relative Risk") +
  scale_y_continuous(name = "Latitude") +
  scale_x_continuous(name = "Longitude") +
  theme(axis.title = element_text(face = "bold", size = 6),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 5),
        strip.text.x = element_text(size = 6),
        legend.title = element_blank()) +
  ggtitle("Relative Risk of Bacterial Blight for India") +
  coord_map("lambert", lat0 = 6.755997, lat1 = 33.17194)
ggsave("Maps/IND_BB.png", width = 6, height = 6, units = "in")

# BS
ggplot(data = IND.BS.df, aes(long, lat, group = group, fill = BS)) +
  geom_polygon(color = "white", size = 0.2) +
  scale_fill_brewer(palette = "GnBu", name = "Relative Risk") +
  scale_y_continuous(name = "Latitude") +
  scale_x_continuous(name = "Longitude") +
  theme(axis.title = element_text(face = "bold", size = 6),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 5),
        strip.text.x = element_text(size = 6),
        legend.title = element_blank()) +
  ggtitle("Relative Risk of Brown Spot for India") +
  coord_map("lambert", lat0 = 6.755997, lat1 = 33.17194)
ggsave("Maps/IND_BS.png", width = 6, height = 6, units = "in")

# LB
ggplot(data = IND.LB.df, aes(long, lat, group = group, fill = LB)) +
  geom_polygon(color = "white", size = 0.2) +
  scale_fill_brewer(palette = "GnBu", name = "Relative Risk") +
  scale_y_continuous(name = "Latitude") +
  scale_x_continuous(name = "Longitude") +
  theme(axis.title = element_text(face = "bold", size = 6),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6),
        plot.title = element_text(face = "bold", size = 8),
        legend.text = element_text(size = 5),
        strip.text.x = element_text(size = 6),
        legend.title = element_blank()) +
  ggtitle("Relative Risk of Leaf Blast for India") +
  coord_map("lambert", lat0 = 6.755997, lat1 = 33.17194)
ggsave("Maps/IND_LB.png", width = 6, height = 6, units = "in")


#eos
