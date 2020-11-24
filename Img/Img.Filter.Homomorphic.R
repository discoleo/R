#########################
###
### Image Processing
### Homomorphic Filtering
###
### Leonard Mada
### draft v.0.1


# based on:
# 1.) https://en.wikipedia.org/wiki/Homomorphic_filtering
# 2.) http://www.worldcolleges.info/sites/default/files/enggnotes/homomorphic-processing-and-its-application-to-image-enhancement-.pdf


##################

### R Code

# install.packages("EBImage")

library(EBImage)


setwd("C:\\Users\\Leo Mada\\Desktop\\Practica/Img")

###################

### Load image
# img.str = "Brugia_40x.jpg"
img.str = "Circuit_BlurAdaptive.jpg"
img.str = "Houses-2085752-0F6D963700000578-924_964x642.jpg"
# Tested on:
# - 1st image with houses in London:
# https://www.dailymail.co.uk/news/article-2085752/Eye-sky-The-pretty-patterns-housing-developments-world.html


img = readImage(img.str)


######################

######################
### Img Processing ###

### Log-Transform
img.log = log(1 + img)


img.dim = dim(img.log)
img.dim


### FFT
# native fft seems equivalent to fftw2d
# & is fully automatic;
x.fft = fft(img.log)

##############

##############
### Filter ###

### generate Filter
filter.f = function(dim, D0=40, c=0.5, H=2, L=0.5) {
	# c is redundant;
	diff = H - L;
	x0.d = dim[1] / 2; y0.d = dim[2] / 2;
	dist.gauss = function(x, y) {
		# L + diff * (1 - exp(- c * ((x-x0.d)^2 + (y-y0.d)^2) / D0^2))
		L + diff * (1 - exp(- c * (x^2 + y^2) / D0^2))
	}
	x0 = dim[1] %/% 2; y0 = dim[2] %/% 2;
	id.x = if(dim[1] %% 2 == 0) c(x0:1, 1:x0) else c(x0:1, 1, 1:x0);
	id.y = if(dim[2] %% 2 == 0) c(y0:1, 1:y0) else c(y0:1, 1, 1:y0);
	outer(id.x, id.y, dist.gauss)
}
### apply Filter
filter.himg = function(img.fft, flt) {
	img.dim = dim(img.fft)
	if(length(img.dim) > 2) {
		# multiple channels
		img.flt = img.fft;
		sapply(seq(img.dim[3]), function(ch) {
			img.flt[,, ch] <<- flt * img.fft[,, ch];
			ch; # print processed channel;
		})
	} else {
	# 1 channel
		img.flt = flt * img.fft
	}
	img.flt;
}

#################

### Create Filter
flt = filter.f(img.dim[1:2], D0=10)

# display Filter:
flt[1:10, 1:10]

display(flt)


### Apply Filter
img.flt = filter.himg(x.fft, flt)


### Convert back to Image

### Inverse FFT:
img2 = Re(fft(img.flt, inverse=1)/prod(img.dim[1:3]))

### Reverse Log-Transform
img2 = exp(img2) - 1
img2[img2 < 0] = 0

### Corrected image:
# - looks better (even better than the original)
#   & seems to be the improved original;
# - normalizing the initial filter by 1/H (= 1/2)
#   does NOT achive this effect!
# - filtered image: img2 = img2 - (L - 1)*img;
img2 = img2 - img;
# alternative:
# img2 = img2 - L*img; # may be useful as well;


###############
### Display ###

# img  = original;
# img2 = Homomorphic filtered;
# img2 - img = difference;
# either img2 or the difference (img2 - img) can be useful;
display(combine(img, img2, img2-img))
# variant: tested with the houses image;
display(combine(img, img2-img, img2))


### Save image
# png(file = "FFT.Homomorphic.Circuits.png", width=360, height=600)
# png(file = "FFT.Homomorphic.Houses.png", width=360, height=600)
old.par = par(mfrow = c(3,1))
plot(img)
plot(img2)
plot(img2 - img)

par(old.par)

# dev.off()


