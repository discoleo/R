########################
###
### Leonard Mada
### [the one and only]
###
### Image Processing: Tools
###
### draft v.0.1a


### History

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
conv.grey = function(img, m, center) {
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
		pix_line = pix.line; # copy
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
			# print(y)
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

### Test
# Center of the kernel:
# - can be asymmetric;
k.center = dim(m) %/% 2;

# gray-scale
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


