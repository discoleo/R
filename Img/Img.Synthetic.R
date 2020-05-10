
################
###
### Generator: Synthetic Images
### e.g. to test Sobel Filters
###
### Leonard Mada
###
### draft v 0.2

# Generate synthetic images
# to test various graphical algorithms:
# - Rectangles;
# - Diagonal polygons;
# - Circles;


# https://github.com/discoleo/R/tree/master/Img/Img.Synthetic.R

################

# install.packages("plotrix")

library(plotrix)


# Save images:
setwd("...")


###############

### helper Functions

# Colors
col.gen = function(lumi) {
	col.base.gen = function(lumi) {
		rgb(lumi, lumi, lumi)
	}
	col.max = max(lumi)
	if(col.max > 1) {
		lumi = lumi / 255
	}
	cols = sapply(lumi, col.base.gen)
	return(cols)
}

################

### Rectangles

synth.rect = function(total, sensitivity=50, width=10, col, save=FALSE) {
	id = 1:total
	if(missing(col)) {
		cols = col.gen(sensitivity*(id - 1))
	} else {
		cols=col
	}
	if(save) {
		file = paste("Test.Sobel.Rect.Grey.", total, ".png", sep="", collapse="")
		# TODO: dimensions of png
		png(file=file, bg = "white")
	}
	par.old = par(mar=c(0.5,0.5,0,0))
	plot.new()
	lim.max = (total + 2) * width
	# Init window
	plot.window(xlim=c(0,lim.max), ylim=c(0,lim.max))
	
	coord = width * (total - id)
	coord_next = coord + width
	# Rect
	# x.L, y.B, x.R, y.T
	rect(0, coord, 20, coord_next, col=cols, border=NA)
	shift.y = 3
	rect(coord, lim.max-10-shift.y, coord_next, lim.max - shift.y, col=cols, border=NA)

	shift.x = total * width
	shift.y = -width
	for(p.id.base in id) {
		p.id = p.id.base
		polygon(shift.x + width*c(0, 1, 1, 0),
			lim.max - shift.y - width*(p.id + c(2,2,1,2)), col=cols[p.id.base], border=NA)
		p.id = p.id + 1
		polygon(shift.x + width*c(0, 0, 1, 0),
			lim.max - shift.y - width*(p.id + c(1,2,1,1)), col=cols[p.id.base], border=NA)
	}

	par(par.old)
	if(save) {
		dev.off()
	}
}

# Diagonal Rectangles/Polygons
synth.diag.rect = function(total, sensitivity=50, width=10, col, save=FALSE) {
	id = 1:total
	if(missing(col)) {
		cols = col.gen(sensitivity*(id - 1))
	} else {
		cols=col
	}
	if(save) {
		file = paste("Test.Sobel.Rect.D.Grey.", total, ".png", sep="", collapse="")
		# TODO: dimensions of png
		png(file=file, bg = "white")
	}
	par.old = par(mar=c(0.5,0.5,0,0))
	plot.new()
	lim.max = (total + 2) * width
	# Init window
	plot.window(xlim=c(0,lim.max), ylim=c(0,lim.max))
	
	# Rect
	# x.L, y.B, x.R, y.T
	draw.diag = function(shift.x, shift.y, pattern) {
		for(p.id.base in id) {
			p.id = p.id.base
			polygon(shift.x + width * pattern[,1],
				shift.y + width*(p.id - 1 + pattern[,2]), col=cols[p.id.base], border=NA)
			p.id = p.id + 1
			polygon(shift.x + width * pattern[,3],
				shift.y + width*(p.id -1 + pattern[,4]), col=cols[p.id.base], border=NA)
		}
	}
	
	pattern = matrix(
	c(0,1,1,0, # X
	  1,1,0,1, # Y
	  0,0,1,0, # X
	  0,1,0,0) # Y
	, nrow=4)
	shift.x = 1
	shift.y = 0
	draw.diag(shift.x, shift.y, pattern)
	shift.x = shift.x + width + 5
	shift.y = 0
	draw.diag(shift.x, shift.y, pattern[,c(3,1,1,3)])
	
	pattern = matrix(
	c(0,1,1,0, # X
	  2,0,1,2, # Y
	  0,0,1,0, # X
	  2,1,0,2) # Y
	, nrow=4)
	shift.x = shift.x + width + 5
	shift.y = 0
	draw.diag(shift.x, shift.y, pattern)
	shift.x = shift.x + width + 5
	shift.y = 0
	pattern.inv = pattern[,c(3,4,1,2)]
	pattern.inv[,c(2,4)] = 2 - pattern.inv[,c(2,4)]
	# print(pattern.inv)
	draw.diag(shift.x, shift.y, pattern.inv)
	

	par(par.old)
	if(save) {
		dev.off()
	}
}

################

### Parameters
total = 5
id = 1:total
# Colors
sesitivity = 50
cols = col.gen(sesitivity*(id - 1))

# Rectangles
synth.rect(total, col=cols)
# synth.rect(total, col=cols, save=TRUE)


### Diagonals
# TODO:
# various diagonals;

###
total = 4

synth.diag.rect(total)
# synth.diag.rect(total, save=TRUE)

###
total = 5

synth.diag.rect(total)
# synth.diag.rect(total, save=TRUE)


############

### Circles

synth.circle = function(total, r=5, seinsitivity=50, col, save=FALSE) {
	id = 1:total
	if(missing(col)) {
		cols = col.gen(sesitivity*(id - 1))
	} else {
		cols=col
	}
	if(save) {
		file = paste("Test.Sobel.Circles.Grey.", total, ".png", sep="", collapse="")
		# TODO: proper dimensions of png
		png(file=file, bg = "white", width=600)
	}
	par.old = par(mar=c(0.5,0.5,0,0))
	plot.new()
	
	lim.max = 2 * total * r + r + 2
	plot.window(xlim=c(0, 2*lim.max), ylim=c(0, lim.max))
	draw.circle(  r + lim.max/3, lim.max /2, r*rev(id), col=cols, border=NA)
	draw.circle(2*r + lim.max*4/3, lim.max /2, r*rev(id), col=rev(cols), border=NA)

	par(par.old)
	if(save) {
		dev.off()
	}
}

### Parameters
total = 4
# Circles
synth.circle(total)
# synth.circle(total, save=TRUE)

### Parameters
total = 5
# Circles
synth.circle(total)
# synth.circle(total, save=TRUE)

