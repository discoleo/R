
### Leonard Mada
###
### Fractals
### draft v 0.1b


### TODO:
# Rachel Skipper: Finiteness properties for simple groups
# https://www.youtube.com/watch?v=AWgoI5jspYc
# - is it meaningful to apply various formalisms?

inner.poly = function(center, width, isH=TRUE, alpha=1/3) {
	if(isH) {
		dy = width * alpha
		x = c(center[1] - width, center[1], center[1] + width, center[1], center[1] - width)
		y = c(center[2], center[2] + dy, center[2], center[2] - dy, center[2])
	} else {
		dx = width * alpha
		x = c(center[1], center[1] - dx, center[1], center[1] + dx, center[1])
		y = c(center[2] - width, center[2], center[2] + width, center[2], center[2] - width)
	}
	
	polygon(x, y, col="white", border=NA)
}

###############

### Fractals

# Init plot
par.old = par(mar=c(0.5,0.5,1.5,1.5))
plot.new()
plot.window(xlim=c(0.5,5.5), ylim=c(.5,3.5))
polygon(c(1,1,5,5), c(1,3,3,1), col="black")

# Problems:
# Q: What is the Fractal Dimension?
# small axis:
# - extension possible only for 1 generation;
# - OR different scaling for elements on small axis: 1/4 x 1/2 of original scaling;


# Fractal

# TODO:
# - correct some of the positions (on L2, L3);
# - see also the problems with the small axis;

# L1
inner.poly(c(3, 2), 1)
# L2
inner.poly(c((1+2)/2, 2), 0.5, isH=F)
inner.poly(c((4+5)/2, 2), 0.5, isH=F)
inner.poly(c(3, (2.5+3)/2), 0.5)
inner.poly(c(3, (1.5+1)/2), 0.5)
# L3
inner.poly(c((1+2)/2, 1+0.5/2), 0.25)
inner.poly(c((1+2)/2, 3-0.5/2), 0.25)
inner.poly(c((4+5)/2, 1+0.5/2), 0.25)
inner.poly(c((4+5)/2, 3-0.5/2), 0.25)
inner.poly(c((2+2.5)/2, 3-0.5/2), 0.25, isH=F)
inner.poly(c((3.5+4)/2, 3-0.5/2), 0.25, isH=F)
inner.poly(c((2+2.5)/2, 1+0.5/2), 0.25, isH=F)
inner.poly(c((3.5+4)/2, 1+0.5/2), 0.25, isH=F)
# + another 4 pieces for L3;


par(par.old)

