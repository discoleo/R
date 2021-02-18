########################
###
### Leonard Mada
### [the one and only]
###
### Image Processing: Tools
###
### draft v.0.2c


### History

### draft v.0.2c:
# - wapply: apply to window;
#   e.g. MaxPool done correctly;
### draft v.0.2a - v.0.2b:
# - Texture analysis:
#   initial experiments;
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
img.str = "Houses-2085752-0F6D963700000578-924_964x642.jpg"
img.str = "Circuit_BlurAdaptive.jpg"
img.str = "Corn_Maize_Corncobs.jpg"

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


####################
####################

### Integral Image


### Method with cumsum():
# img.int = array(0, dim(img)[1:2])
# for(y in seq(dim(img)[2])) {
#	sum = cumsum(as.vector(img[ , y, 1]))
#	img.int[ , y] = sum;
# }


integral.img = function(img, ch=1) {
	img.int = apply(img[,, ch], 2, cumsum)
	tmp = apply(img.int, 1, cumsum)
	tmp = rbind(0, tmp[ - dim(img)[2], ])
	# tmp = rbind(rep(0, dim(img)[1]), tmp[ - dim(img)[2], ])
	img.int = img.int + t(tmp)
	return(img.int)
}
shift.m = function(m, nrow=0, ncol=0, val=0) {
	m = m[1:(nrow(m) - nrow), 1:(ncol(m) - ncol)]
	m = rbind(matrix(val, nrow=nrow, ncol=ncol(m)), m)
	m = cbind(matrix(val, ncol=ncol, nrow=nrow(m)), m)
	m
}

###
img.int = integral.img(img, ch=2)


### Average: 4x7 window
wd.h = 4; wd.w = 7;
img.mean = img.int +
	- cbind(matrix(0, ncol=wd.w, nrow=nrow(img.int)), img.int[, seq(ncol(img.int) - wd.w)]) +
	- rbind(matrix(0, nrow=wd.h, ncol=ncol(img.int)), img.int[seq(nrow(img.int) - wd.h), ]) +
	+ shift.m(img.int, wd.h, wd.w)
display(img.mean / max(img.mean))


### Other Test
x.min = 5
y.min = 4
img.int[10,10] - img.int[10, y.min] - img.int[x.min,10] + img.int[x.min, y.min]
sum(img[(x.min+1):10, (y.min+1):10, 1])


### initial Method:
### 2D
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


#######################
#######################


### Texture Analysis

diff.offset = function(img, off) {
	r.dim = dim(img)[1:2] - off;
	r.m = matrix(0, nrow=r.dim[1], ncol=r.dim[2]);

	for(y in seq(r.dim[2])) {
		for(ch in seq(dim(img)[3])) {
			id = seq(r.dim[1])
			r.m[id,y] = r.m[id,y] + abs(img[id, y, ch] - img[id + off[1], y + off[2], ch])
		}
	}
	r.m = round(255 * r.m / 3)
	invisible(r.m)
}
pad = function(img, dim, val=0) {
	img = rbind(img, matrix(val, nrow=dim[1], ncol=dim(img)[2]));
	img = cbind(img, matrix(val, ncol=dim[2], nrow=dim(img)[1]));
	invisible(img)
}
table.img = function(img, img.diff, ch=1, div=16) {
	isImg1 = missing(img.diff);
	d.dim = if(isImg1) dim(img)[1:2] else dim(img.diff)[1:2];
	img.part = if(isImg1 && length(dim(img)) == 2) img else
		img[seq(d.dim[1]), seq(d.dim[2]), ch];
	if(isImg1) return(table(round(img.part/div)));
	table(round(img.part*255/div), round(img.diff/div))
}


off = c(5, 3)
img.diff = diff.offset(img, off)
image(img.diff)


d.m = table.img(img, img.diff, ch=1)
median(d.m)
median(d.m[d.m > 0])
mean(d.m)

table.img(img.diff)

image(img[,,2] - pad(img.diff, off)/255)


###################

### WAPPLY

wapply = function(d, window, FUN, diagonal, ...) {
	wnd = window;
	if(length(wnd) < 2) wnd = c(wnd, wnd);
	noDiagonal = missing(diagonal);
	hasDiagonal = (! noDiagonal) && any(diagonal != 0);
	f = function(i, j) {
		r = rep(0, length(i))
		if(hasDiagonal) r[i == j] = diagonal;
		isUpper = if(noDiagonal) TRUE else (i > j);
		i2 = i[isUpper]; j2 = j[isUpper];
		r[isUpper] = sapply(seq(length(i2)), function(id)
			FUN(d[i2[id]:(i2[id] + wnd[1]), j2[id]:(j2[id] + wnd[2])], ...))
		return(r)
	}
	d.dim = dim(d)
	outer(seq(d.dim[1] - wnd[1]), seq(d.dim[2] - wnd[2]), f)
}

wnd = c(5, 5)
# "MaxPool" done correctly
tx.max = wapply(img.diff, wnd, FUN=max)
summary(as.vector(tx.max))

image(tx.max)
display(tx.max/255)

