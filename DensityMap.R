
### Libraeies  ####
library(readr)
library(ggplot2)
library(scales)
library('png')

### Load Data ####

rawData <-  read_table2("~/Desktop/Projects/2018_Block_Typologies/Data/AllClusterLocations1MItr_5.csv")
rawData$x = rawData$x +1
rawData$y = rawData$y +1

somDensity <- data.frame(matrix(0,nrow = max(rawData$y, na.rm=TRUE),ncol = max(rawData$x, na.rm=TRUE)));
clusterData <- data.frame(matrix(0,nrow = max(rawData$y, na.rm=TRUE),ncol = max(rawData$x, na.rm=TRUE)));
img <- readPNG("images/cubediagonal.png")
img1 <- readPNG('images/bremm.png')
img2 <- readPNG('images/teulingfig2.png')
img3 <- readPNG('images/cubediagonal.png')

### function ####
### get colour from image

getColour2D3D = function (i,j){
  pixel = img[1,512,]
  #x=1
  #i = as.integer(i*5.12)
  #j = as.integer(j*5.12)
  pixel = img[j*5.12,i*5.12,]
  #pixel[1]*255
  x = paste0(pixel[1]*255, sep=",",pixel[2]*255,sep=",",pixel[3]*255)
  pixel = img2[j*5.12,i*5.12,]
  x = paste0(x,sep=",",pixel[1]*255, sep=",",pixel[2]*255,sep=",",pixel[3]*255)
  #pixel[1]*255
  #x = paste0(x,sep=",",pixel[1]*255, sep=",",pixel[2]*255,sep=",",pixel[3]*255)
  return(x)
}
rawData$'colourR,colourG,colourB,TeuColourR,TeuColourG,TeuColourB' = 1
rawData$'colourR,colourG,colourB,TeuColourR,TeuColourG,TeuColourB' = mapply(getColour2D3D,rawData$x,rawData$y)
colnames(rawData)[1] <- "image1,image2"
write.csv(rawData, "~/Desktop/Projects/2018_Block_Typologies/Data/trained_Som_classified_images_melsyd1MItr_I.csv", row.names = TRUE, quote = FALSE)
write.csv(clusterData, "~/Desktop/Projects/2018_Block_Typologies/Data/clusterXY.csv", row.names = FALSE, quote = FALSE)

### Data Prep ####
for (i in 1:nrow(rawData)){
  #print(i)
  somDensity[rawData[i,]$x,rawData[i,]$y]=somDensity[rawData[i,]$x,rawData[i,]$y]+1;
  clusterData[rawData[i,]$x,rawData[i,]$y]=rawData[i,]$cluster;
}  
View(somDensity);

finalData = matrix(1,nrow = 1,ncol = 4)

for (i in 1:nrow(somDensity)){
  for(j in 1:ncol(somDensity)){
    newRow = c(i,j,somDensity[i,j],clusterData[i,j])
    finalData =rbind(finalData,newRow)
  }
}  
finalData = as.data.frame(finalData)
names(finalData) = c("x","y","weight","cluster")
finalData = finalData[-1,] 
finalData = subset(finalData, finalData$weight>0)
finalData$x = finalData$x-1
finalData$y = finalData$y-1

### Plotting ####

ggplot(finalData, aes(cluster)) + geom_histogram()
ggsave("som_histogram.png", width = 6, height = 5)

ggplot(finalData, aes(x,y)) + geom_raster(aes(fill = weight), interpolate = FALSE) 
ggsave("som_Density.png", width = 6, height = 5)

ggplot(finalData, aes(x, y, fill = log(weight))) + geom_raster() + scale_fill_gradient(low = "red",high = "black") 
ggsave("som_Density_log.png", width = 6, height = 5)

ggplot(finalData, aes(x, y, fill = cluster)) + geom_raster() + scale_fill_gradientn(colours=c("#0000FFFF", "#f4d03f","#2e4053","#7d6608","#FF0000FF","#0000FFFF"))
ggsave("som_clusters_.png", width = 6, height = 5)

ggplot(rawData, aes(x=x)) + geom_density()
ggplot(rawData, aes(x=y)) + geom_density()


ggplot(finalData, aes(x,y)) + 
  geom_raster(aes(fill = weight), interpolate = FALSE) +
  scale_fill_gradientn( colours=c("#0000FFFF", "#f4d03f","#2e4053","#7d6608","#FF0000FF"), limits=c(0,4))


### get cluster centres ####

### urlStringImputBase = "http://10.100.8.30:5002/~thud/projects/city_typology/download_new_style/"

for(i in 0:max(rawData$cluster)){
  tempData = subset(rawData, rawData$cluster == i)
  tempDataOne = tempData[which.min(tempData$distance),]
  print(i)
  print(tempDataOne$image)
  #page = read_html(urlStringImput)
}



### hotspots in SOM ###


minClipNumber = 10000
tempSomDensity = subset(finalData, finalData$weight>minClipNumber)

for(i in 1 :nrow(tempSomDensity)){
  highDensityPoint = which(rawData$x == tempSomDensity[i,]$x & rawData$y ==  tempSomDensity[i,]$y)
  highDensityPoint = highDensityPoint[1]
  rawData[highDensityPoint,]$image
  #highDensityPoint = c(tempSomDensity[i,]$x,tempSomDensity[i,]$y)
  print(tempSomDensity[i,]$cluster)
  print(rawData[highDensityPoint,]$image)
}