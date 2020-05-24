
################
###
### Generator: Synthetic Images
### e.g. to test Sobel Filters
###
### Leonard Mada
###
### draft v 0.3f

# Generate synthetic images
# to test various graphical algorithms:
# - Rectangles;
# - Diagonal polygons;
# - Circles;
# - Cells;

# https://github.com/discoleo/R/tree/master/Img/Img.Synthetic.R


##################

# install.packages("plotrix")

library(plotrix)


# Save images:
setwd("\\img")


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
# only on Gray-levels
add.col = function(fill, base.col, n=length(base.col)) {
	# cols = original colors;
	if(is.list(fill)) {
		if(length(fill$rel) != 0) {
			# TODO: How should rel work?
			fill = rgb(t(col2rgb(base.col)/255 * fill))
		} else if(length(fill$add) != 0) {
			fill = rgb(t(col2rgb(base.col)/255 + fill))
		}
	} else if(length(fill) == 1) {
		if(fill == TRUE) {
			fill = base.col
		} else if(fill < 0) {
			# fill = base.col + fill
			fill = rgb(t(col2rgb(base.col)/255 + fill))
		} else {
			fill = rep(fill, n)
		}
	}
	return(fill)
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

synth.circle = function(total, r=5, sensitivity=50, col, save=FALSE) {
	id = 1:total
	if(missing(col)) {
		cols = col.gen(sensitivity*(id - 1))
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



#######################
#######################

#######################
### Test Cells (Biology)

### Generator Function
synth.cells = function(grid, r, width, r.jitter = r,
		col.opt = list(sens=100, var=c(0.5, 1.5)), cols, fill=NA, save=FALSE, img.id=NA,
		margin=FALSE, alternating=FALSE) {
	#
	n = grid$n * grid$n
	if(missing(cols)) {
		cols.palette = col.gen(col.opt$sens*runif(n, col.opt$var[1], col.opt$var[2]))
		cols = sample(cols.palette, n, replace=TRUE)
	}
	isFill = TRUE
	if( ! missing(fill) ) {
		fill = add.col(fill, cols)
	} else {
		isFill = FALSE
	}
	# Centers
	id.x = runif(n, -r.jitter, r.jitter)
	id.y = runif(n, -r.jitter, r.jitter)
	grid.all = rep(1:grid$n, grid$n)
	c.x = (grid.all) * grid$d + id.x
	c.y = (grid.all) * grid$d + id.y
	c.y = c(t(matrix(c.y, nrow=grid$n)))
	if(alternating) {
		id.between = (n %/% grid$n) %% 2 == 1
		c.x[id.between] = c.x[id.between] - grid$d / 2
	}

	if(save) {
		img.id = if(missing(img.id)) "" else paste(".ID", img.id, sep="")
		file = paste("Test.Cells.Grey", img.id, ".N_", n, ".R_", r, ".png", sep="", collapse="")
		# TODO: proper dimensions of png
		png(file=file, bg = "white") # width = ?
	}
	### Draw Circles
	if(margin) {
		par.old = par(mar=c(0.5,0.5,0,0))
	} else {
		par.old = par(mar=c(0,0,0,0))
	}
	plot.new()
	
	lim.max = grid$n * grid$d + r
	plot.window(xlim=c(0, lim.max), ylim=c(0, lim.max))
	for(id in 1:n) {
		draw.circle(  c.x[id], c.y[id], r, col=ifelse(isFill, fill[id], NA), border=cols[id], lwd=width)
	}

	par(par.old)
	if(save) {
		dev.off()
	}
}

# TODO:
# high lumi fill, alternating grid rows;
# image background;
# populations of cells: fill, R;

test.cells.gen = function(samples, grid, r=6, width=4, save=FALSE) {
	id.series = 1:samples
	base.lumi = 100
	SAVE = save
	
	# Baseline
	series.id = id.series
	for(id in id.series) {
		synth.cells(grid, r, width, img.id=series.id[id], save=SAVE)
	}
	# Luminosity Drift
	series.id = series.id + samples
	diff = seq(0, by=0.8/samples, length.out=samples)
	for(id in id.series) {
		col.opt = list(sens=base.lumi, var=c(1 - diff[id], 1+ diff[id]))
		synth.cells(grid, r, width, col.opt=col.opt, img.id=series.id[id], save=SAVE)
	}
	# Luminosity Base
	series.id = series.id + samples
	diff = (base.lumi - 20)/samples
	lumi = id.series * diff
	for(id in id.series) {
		col.opt = list(sens=lumi[id], var=c(0.5, 1.5))
		synth.cells(grid, r, width, col.opt=col.opt, img.id=series.id[id], save=SAVE)
	}
	# Radius
	series.id = series.id + samples
	diff = 2*r/samples
	r.seq = seq(diff, by=diff, length.out=samples)
	for(id in id.series) {
		synth.cells(grid, r=r.seq[id], width, img.id=series.id[id], save=SAVE)
	}
}

test.cells.gen(5, dim.grid, save=TRUE)


######################

grid.n = 12 # 10
grid.d = 15 # 18
#
dim.grid = list(n = grid.n, d = grid.d)
#
r = 5 # 6
width = 4
# hole = r - width

synth.cells(dim.grid, r, width=width)

# synth.cells(dim.grid, r, width=width, save=T)

synth.cells(dim.grid, r, width=width, fill=-0.1)


synth.cells(dim.grid, r, width=width, alternating=TRUE)

