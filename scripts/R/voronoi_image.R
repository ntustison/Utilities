doInstall <- TRUE  # Change to FALSE if you don't want packages installed.
toInstall <- c("ReadImages", "reshape", "ggplot2")
if(doInstall){install.packages(toInstall, repos = "http://cran.r-project.org")}
lapply(toInstall, library, character.only = TRUE)

# Image URL:
allImageURLs <- c("http://media.charlesleifer.com/blog/photos/thumbnails/akira_940x700.jpg",
                  "http://upload.wikimedia.org/wikipedia/commons/thumb/e/ec/Mona_Lisa%2C_by_Leonardo_da_Vinci%2C_from_C2RMF_retouched.jpg/402px-Mona_Lisa%2C_by_Leonardo_da_Vinci%2C_from_C2RMF_retouched.jpg",
                  "http://upload.wikimedia.org/wikipedia/commons/thumb/e/e9/Official_portrait_of_Barack_Obama.jpg/441px-Official_portrait_of_Barack_Obama.jpg",
                  "http://cache.boston.com/universal/site_graphics/blogs/bigpicture/obama_11_05/obama22_16604051.jpg",
                  "http://upload.wikimedia.org/wikipedia/commons/thumb/e/ea/Van_Gogh_-_Starry_Night_-_Google_Art_Project.jpg/758px-Van_Gogh_-_Starry_Night_-_Google_Art_Project.jpg",
                  "http://www.10mfh.com/wp-content/uploads/2011/09/dino_riders.jpg",
                  "http://images3.alphacoders.com/855/8557.jpg",
                  "http://inapcache.boston.com/universal/site_graphics/blogs/bigpicture/ngm_101912/bp19.jpg",
                  "http://inapcache.boston.com/universal/site_graphics/blogs/bigpicture/ngm_101912/bp26.jpg",
                  "http://inapcache.boston.com/universal/site_graphics/blogs/bigpicture/ngm_101912/bp35.jpg",
                  "http://inapcache.boston.com/universal/site_graphics/blogs/bigpicture/balloon/bp6.jpg",
                  "http://inapcache.boston.com/universal/site_graphics/blogs/bigpicture/smithsonian_030512/bp14.jpg",
                  "http://inapcache.boston.com/universal/site_graphics/blogs/bigpicture/smithsonian_030512/bp15.jpg",
                  "http://inapcache.boston.com/universal/site_graphics/blogs/bigpicture/earth_day_2012/bp6.jpg",
                  "http://inapcache.boston.com/universal/site_graphics/blogs/bigpicture/2011part2/bp1.jpg",
                  "http://inapcache.boston.com/universal/site_graphics/blogs/bigpicture/2011part2/bp4.jpg",
                  "http://inapcache.boston.com/universal/site_graphics/blogs/bigpicture/2011part2/bp15.jpg",
                  "http://inapcache.boston.com/universal/site_graphics/blogs/bigpicture/2011part2/bp27.jpg",
                  "http://inapcache.boston.com/universal/site_graphics/blogs/bigpicture/natural_world_2011/bp40.jpg",
                  "http://inapcache.boston.com/universal/site_graphics/blogs/bigpicture/ngmphotocontest_111811/bp10.jpg",
                  "http://inapcache.boston.com/universal/site_graphics/blogs/bigpicture/ngmphotocontest_111811/bp54.jpg")

imageLoader <- function(url){  # This function takes a URL, and generates a data.frame with pixel locations and colors
  # Download to disk, load
  download.file(url, "tempPicture.jpg", mode = "wb")  # Stash image locally
  readImage <- read.jpeg("tempPicture.jpg")

  longImage <- melt(readImage)
  rgbImage <- reshape(longImage, timevar = "X3",
                      idvar = c("X1", "X2"), direction = "wide")
  rgbImage$X1 <- -rgbImage$X1
  return(rgbImage)
  }

##########
# Part 3 # Completely, unbelievably awesome Voronoi clusters:
##########

rgbImage <- imageLoader(allImageURLs[5])  # Pick one, or use your own URL.

nRegions <- 500  # Number of voronoi regions, large numbers are slow.
voronoiMeans <- kmeans(rgbImage, centers = nRegions, iter.max = 50)
voronoiColor <- voronoiMeans$centers[voronoiMeans$cluster, 3:5]
with(rgbImage, plot(X2, X1, col = rgb(voronoiColor), asp = 1, pch = "."))

rgbImage <- imageLoader(allImageURLs[5])  # Pick one, or use your own URL.
rgbImage[, 1:2] <- sweep(rgbImage[, 1:2], 2, apply(abs(rgbImage[, 1:2]), 2, max), "/")
rgbImage[, 1:2] <- rgbImage[, 1:2] * 1

nRegions <- 100  # Number of voronoi regions, large numbers are slow.
voronoiMeans <- kmeans(rgbImage, centers = nRegions, iter.max = 50)
voronoiColor <- voronoiMeans$centers[voronoiMeans$cluster, 3:5]
with(rgbImage, plot(X2, X1, col = rgb(voronoiColor), asp = 1, pch = "."))
