#.................................
# Dependencies
#.................................
library(tidyverse)
library(sf)
library(raster)
library(rgeos)
library(rgdal)

#..............................
# Pull down the border cntrs w/ raster
#..............................
brdrcnt <- lapply(c("UGA", "SSD", "CAF", "COG", "AGO", "ZMB", "TZA", "RWA", "BDI", "GAB", "CMR", "GNQ",
                    "NGA", "TCD", "SDN", "BEN", "TGO", "GHA", "NER", "BFA", "CIV", "COD", "ETH",
                    "KEN", "MWI", "MOZ", "ZWE"),
                  function(x){
                    ret <- sf::st_as_sf(rgeos::gSimplify(
                      raster::getData(name = "GADM", country = x, level = 0, path = "data/map_bases/"), tol = 1e-2))
                    return(ret)

                  })


basemap <- list(
  geom_sf(data = brdrcnt[[1]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[2]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[3]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[4]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[5]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[6]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[7]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[8]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[9]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[10]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[11]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[12]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[13]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[14]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[15]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[16]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[17]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[18]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[19]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[20]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[21]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[22]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[23]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[24]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[25]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[26]], fill = "#f0f0f0", lwd = 0.5),
  geom_sf(data = brdrcnt[[27]], fill = "#f0f0f0", lwd = 0.5),
  coord_sf(xlim = c(-1, 38.5), ylim = c(-16, 10.8), datum = sf::st_crs(4326)),
  theme(panel.background = element_rect(fill = "#9ecae1"),
        panel.grid = element_line(colour="transparent"),
        axis.text = element_blank())
)

saveRDS(basemap, file = "data/map_bases/basemap.rds")
