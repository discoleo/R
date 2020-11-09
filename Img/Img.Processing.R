########################
###
### Leonard Mada
### [the one and only]
###
### Image Processing: Tools
###
### draft v.0.1c


### History

### draft v.0.1c:
# - "Adaptive" Blur / Quasi-Non-Local Means;
### draft v.0.1b:
# - added the fft-based convolution;
### draft v.0.1a: [09-11-2020]
# - classical pixel-wise convolution;

# Note:
# - the specific algorithms were mostly developed over the previous years,
#   but never oficially published;
# - TODO:
#  -- add all algorithms;
#  -- native R impleemntations are much slower than the Java or C# variants;


####################

### helper Functions

### Kernel Operations

### classical (pixelwise) convolution
conv.grey = function(img, m, center = dim(m) %/% 2) {
	if(length(dim(img)) > 3) {
		stop("Currently only Greyscale images supported!")
	}
	# Kernel
	k = decompose.kernel(m)
	k.width = max(k$x); k.height = max(k$y)
	k.xstart = min(k$x); k.ystart = min(k$y)
	# New image: assumes base pixels remain the same;
	# TODO: implement processing of margins;
	new.img = img; # copy
	pix.line = rep(0, dim(new.img)[1])

	y.max = dim(new.img)[2] - k.height + k.ystart;
	x.max = dim(new.img)[1] - k.width + k.xstart;
	
	for(y in 1:y.max) {
		pix_line = pix.line; # copy vs instantiate inline?
		sapply(1:nrow(k), function(id) {
			x = 1:x.max + k$x[id] - 1
			y.tot = y - 1 + k$y[id]
			# scan image line & apply convolution
			pix = img[x, y.tot] * k$w[id]
			# add modified line to "cache"
			# TODO: evaluate efficiency of:
			# return(pix) + addition outside of sapply;
			pix_line[center[1]:(x.max + center[1] - 1)] <<-
				pix_line[center[1]:(x.max + center[1] - 1)] + pix;
		} )
		new.img[ , y + center[2]] = pix_line;
		# TODO: Progress
		if(y %% 20 == 1) {
			# cat(y)
		}
	}
	return(new.img)
}

### Quasi-Non-Local Means/Convolution
### aka Similar-Means / "Adaptive" Convolution
conv.adaptive = function(img, m, thresh=1, center = dim(m) %/% 2) {
	if(length(dim(img)) != 3 && dim(img)[3] != 3) {
		stop("Currently only Color images supported!")
		# TODO: grey-scale;
	}
	# Kernel
	k = decompose.kernel(m)
	k.width = max(k$x); k.height = max(k$y)
	k.xstart = min(k$x); k.ystart = min(k$y)
	# New image: assumes base pixels remain the same;
	# TODO: implement processing of margins;
	new.img = img; # copy
	pix.line = matrix(0, nrow=dim(new.img)[1], ncol=3)
	w.line = matrix(0, nrow=dim(new.img)[1], ncol=3)

	y.max = dim(new.img)[2] - k.height + k.ystart;
	x.max = dim(new.img)[1] - k.width + k.xstart;

	for(y in 1:y.max) {
		pix_line = pix.line; # copy vs instantiate inline?
		w_line = w.line; # copy
		pixel = img2[center[1]:(x.max + center[1] - 1), y + center[2], ]
		dim(pixel) = c(x.max, 3)
		sapply(1:nrow(k), function(id) {
			x = 1:x.max + k$x[id] - 1
			y.tot = y - 1 + k$y[id]
			# scan image line
			pix = img2[x, y.tot,]
			dim(pix) = c(x.max, 3)
			# print(dim(pix)); print(dim(pixel))
			diff = abs(pix - pixel)
			diff = diff[,1] + diff[,2] + diff[,3]
			isThresh = (diff <= thresh); # similar pixels
			### Convolve
			pix[,1] = ifelse(isThresh, pix[,1] * k$w[id], 0)
			pix[,2] = ifelse(isThresh, pix[,2] * k$w[id], 0)
			pix[,3] = ifelse(isThresh, pix[,3] * k$w[id], 0)
			# update pixel-values
			pix_line[center[1]:(x.max + center[1] - 1), ] <<-
				pix_line[center[1]:(x.max + center[1] - 1), ] + pix;
			# update weights: the sum changes!
			w_line[center[1]:(x.max + center[1] - 1), ][isThresh, ] <<-
				w_line[center[1]:(x.max + center[1] - 1), ][isThresh, ] + k$w[id];
		} )
		new.img[ , y  + center[2], ] = pix_line / w_line;
	
		# TODO: Progress
		if(y %% 20 == 1) {
			cat(paste0(y, ", "))
		}
	}
	return(new.img)
}

### Kernel Tools
decompose.kernel = function(m) {
	dim = dim(m)
	pos.x = integer()
	pos.y = integer()
	w = numeric()
	for(y in seq(dim[2])) {
		for(x in seq(dim[1])) {
			if(m[x, y] != 0) {
				pos.x = c(pos.x, x)
				pos.y = c(pos.y, y)
				w = c(w, m[x, y])
			}
		}
	}
	return(data.frame("y"=pos.y, "x"=pos.x, "w"=w))
}

####################


# install.packages("EBImage")

library(EBImage)


setwd(".../Img")


#######################

# Load image
img.str = "Brugia_40x.jpg"
img.str = "Circuit_BlurAdaptive.jpg"

img = readImage(img.str)

########################

# Kernel: can have even rows/columns;
m = matrix(c(
	0, 1, 1, 1, 1, 1, 0,
	0, 1, 1, 1, 1, 1, 0,
	1, 1, 1, 3, 1, 1, 1,
	0, 1, 1, 1, 1, 1, 0,
	0, 1, 1, 1, 1, 1, 0
	), ncol=7, byrow=T)
m = m / sum(m)

##############

### Test
# Center of the kernel:
# - can be asymmetric;
k.center = dim(m) %/% 2;

### gray-scale
grey.img = Image(img, colormode='Grayscale')
grey.img = grey.img[,,1]

t.time = Sys.time()
new.img = conv.grey(grey.img, m, k.center);
t.time = Sys.time() - t.time;
t.time
# 8 s with 5x5 kernel;
# 13.4 - 15 s with 7x7 kernel!

display(new.img)

display(new.img - grey.img)

####################

### "Adaptive" Blur
### [Color]
# threshhold:
# = 0.2; # very strict;
# = 1; 1.5; # medium;
# = 2; # extensive blur;

thresh = 2; # 0.75; # 1
t.time = Sys.time()
new.img = conv.adaptive(img, m, thresh);
t.time = Sys.time() - t.time;
t.time
# ... with 5x5 kernel;
# ~2 mins with 7x7 kernel!

display(new.img)

display(new.img - img)


###############

### Brush:
flt = makeBrush(5, shape='disc', step=FALSE)
# display(flt, title='Disc filter')

flt = flt/sum(flt)


###############

###############
### FFT-Way ###

img2 = filter2(img, m, boundary=2)

display(img2)

### explicit
# native fft seems equivalent to fftw2d
# & is fully automatic;
x.fft = fft(img)
dim.len = length(dim(img));

# but needs a transformed kernel!
# img2 = Re(fft(x.fft * fft(_transformed_m_), inverse=1)/prod(dim(img)[1:dim.len]))


