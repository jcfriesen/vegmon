###########################################################################
## Vegetation metric provcessing
# The following code calculates a greenness index and 
# canopy cover from regular photos. As greenness index 
# the VIgreen index is calculated (St. Peter et al., 2018) 
# and the canopy cover is calulated based on Rosin method 
# (Chianucci et al. 2018, Rosin 2001). The input data are 
# two photos taken by the ODK VegMon app. Photo 'p45' is 
# taken at an angle of 45 deg and photo 'p90' is taken at 
# a 90 deg angle (bird's eye view). Both photos are taken 
# with smartphones held approximately at breast height.

# Photo 45 and 90 will both be used to caclulate indivudual VIgreen values
# Photo 90 will, in addition, be used to calculate canopy cover

## Authors: # Shakir AHMED # Jan FRIESEN
## References: 
# St. Peter et al. (2018) Remote Sensing, 10.
# Chianucci et al. (2018) Biosystems Engineering, 169.
# Rosin (2001) Pattern Recognition, 34.
###########################################################################
# Required libraries
library(imagerExtra)

# image location and paths

####### CHANGE IMAGE PATH TO p45.jpg and p90.jpg local image folder
####### CHANGE OUTPUT PATH TO local folder
image_path  <- "...";
out_path  <- "...";

image45_file <- paste(image_path, "p45.jpg",sep = "")
image90_file <- paste(image_path, "p90.jpg",sep = "")

# parameters  
darkvalue=25

## Photo 45
#################################################
# read image and separate into channels
img <- load.image(image45_file) # read RGB file
r<-R(img)
g<-G(img)
b<-B(img)

## Vegetation index (VIgreen)
VIgreen <- ((g-r)/(g+r))

####### Results photo 45 (VIgreen)
vig45_value <- mean(VIgreen,trim=0, na.rm =TRUE) # average VIgreen value for whole p45 image

remove(img, r, g, b, VIgreen)

## Photo 90
#################################################
# read image and separate into channels
img <- load.image(image90_file) # read RGB file
r<-R(img)
g<-G(img)
b<-B(img)

## Vegetation index (VIgreen)
VIgreen <- ((g-r)/(g+r))

####### Results photo 90 (VIgreen)
vig90_value <- mean(VIgreen,trim=0, na.rm =TRUE) # average VIgreen value for whole p45 image

remove(VIgreen)

## Canopy cover (CC_Ros)
#get img dimensions
img_cols <- ncol(img)
img_rows <- nrow(img)

#convert to integer and scale to 0 to 255 
r<-as.integer(as.double(r)*255)
g<-as.integer(as.double(g)*255)
b<-as.integer(as.double(b)*255)

##Identify pixel groups
g_1= ((r+b-g-g)>=0)             #group 1 
g_2=(g<=darkvalue)              #group 2 
g_3=((g>r)&(g>b)&(g>darkvalue)) #group 3
g_4=(!g_3 &!g_1 & !g_2)         #group 4

##Visible-vegetation index (GLA) based on RBG channels
GLA=(g+g-r-b)/(g+g+r+b); #apply algorithm
GLA=as.double((GLA+1)/2)  # set accuracy
GLA_img <-as.cimg(GLA,img_rows,img_cols,1,1)  # convert vector to image 
GLA_bs <- BalanceSimplest(inpaint(GLA_img,1), 1, 1, range = c(0, 255)) # stretch contrast
gla <- round(GLA_bs, digits = 0) # convert to integer

# get histogram of g_4 GLA pixels
hist_gla_g_4=hist(gla[g_4], breaks = seq(0,255,by=1))  #group 4
h_gg4 <- hist_gla_g_4$counts # extract histogram count vector
h_start <- min(which(h_gg4 != 0)) # first non-empty bin position
h_end <- max(which(h_gg4 != 0)) # last bin position
h_maxpos <- which(h_gg4 == max(h_gg4)) # max bin position
h_max <- max(h_gg4) # max bin count

slope <- h_max/h_end-h_maxpos
alpha <- atan(1/slope)

# Rosin (2001) corner detection
H <- matrix(0,1,256)
O <- matrix(0,1,256)
O[,] <- NA
for (i in h_maxpos:h_end) {
  H[1,i] <- slope*(0-i)-h_gg4[i]
  O[1,i]=H[1,i]*sin(alpha);
} 

Ros_thr <- which(O == max(O, na.rm = TRUE)) - 1 # Rosin threshold
Ros <- threshold(gla, thr = Ros_thr) # Binary GLA pixset using Rosin threshold


####### Results photo 90 (CC_Ros)
CC_Ros <- sum(Ros, na.rm=TRUE) / length(Ros) # canopy cover ratio

# canopy cover image export
CC_img <- as.cimg(Ros)
save.image(CC_img, paste(out_path, "cc_rosin_90.png",sep = ""))

remove(image_path, out_path, darkvalue,
       img, r, g, b, img_cols, img_rows,
       g_1, g_2, g_3, g_4, gla, GLA, GLA_bs, GLA_img,
       hist_gla_g_4, h_gg4, h_start, h_end, h_maxpos, h_max,
       slope, alpha, H, O, i, Ros_thr, Ros)


