########################
###
### Leonard Mada
### [the one and only]
###
### Image Processing
### Lab 6: Integral Image
###
### draft v.0.1a

# - only the Lab part;
# - code ported to R;
# - efficient computation of the integral image,
#   using 2 orthogonal cummulative sums;

###################

# install.packages("EBImage")

library(EBImage)


setwd(".../Img")


#######################

# Load image
img.str = "Circuit_BlurAdaptive.jpg"

img = readImage(img.str)

########################



######################
### Integral Image ###

### Method:
# uses cumsum() twice: on each orthogonal direction;

integral.img = function(img, ch=1) {
	img.int = apply(img[,, ch], 2, cumsum)
	tmp = apply(img.int, 1, cumsum)
	tmp = rbind(0, tmp[ - dim(img)[2], ])
	# tmp = rbind(rep(0, dim(img)[1]), tmp[ - dim(img)[2], ])
	img.int = img.int + t(tmp)
	return(img.int)
}
shift.m = function(m, nrow=0, ncol=0, val=0) {
	# shifts a matrix (e.g. the integral image);
	m = m[1:(nrow(m) - nrow), 1:(ncol(m) - ncol)]
	m = rbind(matrix(val, nrow=nrow, ncol=ncol(m)), m)
	m = cbind(matrix(val, ncol=ncol, nrow=nrow(m)), m)
	m
}

### Green channel
img.int = integral.img(img, ch=2)


### Test
### Average: 4x7 window
# - computed over the entire image;
wd.h = 4; wd.w = 7;
img.mean = img.int +
	- cbind(matrix(0, ncol=wd.w, nrow=nrow(img.int)), img.int[, seq(ncol(img.int) - wd.w)]) +
	- rbind(matrix(0, nrow=wd.h, ncol=ncol(img.int)), img.int[seq(nrow(img.int) - wd.h), ]) +
	+ shift.m(img.int, wd.h, wd.w)
display(img.mean / max(img.mean)) # a 4x7 blur!


### Other Test
x.min = 5
y.min = 4
img.int[10,10] - img.int[10, y.min] - img.int[x.min,10] + img.int[x.min, y.min]
sum(img[(x.min+1):10, (y.min+1):10, 1])


###################

### initial Method:

### 2D array/matrix
img.int = array(0, dim(img)[1:2])

for(y in seq(dim(img)[2])) {
	sum = cumsum(as.vector(img[ , y, 1]))
	
	if(y == 1) {
		img.int[ , 1] = sum;
	} else {
		img.int[ , y] = sum + img.int[ , y-1]
	}
	if(y %% 20 == 1) print(y)
}

img.int[10,10] - img.int[10,4] - img.int[4,10] + img.int[4,4]
sum(img[5:10, 5:10, 1])
