setwd("/users/ryan/Desktop/test")
library(jpeg)
img = readJPEG("./lenna.jpg")  ##input a grayscale image as a matrix

##if(length(dim(img)>2)) img = img[,,1]  ##some picture need this.
num_rdn = 0.4*ncol(img)*nrow(img)  ##number of random missing

set.seed(2016000100)
x_rdn = sample.int(n=nrow(img), size=num_rdn,replace = TRUE)
y_rdn = sample.int(n=ncol(img), size=num_rdn,replace = TRUE)

img_rdn = img
for(i in 1:num_rdn){
  img_rdn[x_rdn[i],y_rdn[i]] = NA
}     ##randomly miss

P_Omega = matrix(0,nrow(img_rdn),ncol(img_rdn))
for(i in 1:nrow(img_rdn)){
  for(j in 1:ncol(img_rdn)){
    if(!is.na(img_rdn[i,j])) P_Omega[i,j] = 1
  }
}     ##indicator matrix of missing coordinate



#plot(0,0, type="n")
#rasterImage(img_rdn, -1,-1,1,1)


library("softimpute")

zz = softimpute(img_rdn, P_Omega, 0.5)
zz$num_iteration   ##number of iteration

plot(0,0, type="n")
rasterImage(zz$Z, -1,-1,1,1)

