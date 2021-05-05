###
### Vest University:
### Group Project 2021
###
### Team: Modeling Infection Spread
### Supervisor: Leonard Mada
###
### draft v0.2a


########################
### Models: Diagrams ###
########################

# install.packages("shape")

library(shape)


####################

### helper functions

### line-arrow
lnarrow = function(xy1, xy2, r=r0, head.jt=0.01, lwd=1.5, ...) {
	tol = 0.002;
	if(abs(xy1[2] - xy2[2]) <= tol) {
		if(xy2[1] >= xy1[1]) {
			Arrows(xy1[1] + r, xy1[2], xy2[1] - r - head.jt, xy2[2], lwd=lwd, ...)
		} else {
			Arrows(xy1[1] - r - head.jt, xy1[2], xy2[1] + r, xy2[2], lwd=lwd, ...)
		}
	} else if(abs(xy1[1] - xy2[1]) <= tol) {
		if(xy2[2] >= xy1[2]) {
			Arrows(xy1[1], xy1[2] + r, xy2[1], xy2[2] - r - head.jt, lwd=lwd, ...)
		} else {
			Arrows(xy1[1], xy1[2] - r, xy2[1], xy2[2] + r + head.jt, lwd=lwd, ...)
		}
	} else {
		slope = (xy2[2] - xy1[2]) / (xy2[1] - xy1[1])
		slope.sc = r / sqrt(1 + slope^2);
		rx = slope.sc; ry = slope * slope.sc;
		if(xy1[1] < xy2[1]) {
			Arrows(xy1[1] + rx, xy1[2] + ry,
				xy2[1] - rx - head.jt, xy2[2] - ry - sign(slope)*head.jt, lwd=lwd, ...)
		} else {
			# TODO
			ry = - ry;
			Arrows(xy1[1] - rx, xy1[2] + ry,
				xy2[1] + rx + head.jt, xy2[2] - ry + sign(slope)*head.jt, lwd=lwd, ...)
		}
	}
}
### curved-arrow
cvarrow = function(xy, col, w=w0, r=r0, mid.off = c(0.005, 0.035)) {
	xy = xy - c(r, 0);
	filledcircle(r/2 + w, r/2, mid=xy + mid.off, from=pi/4, to=-pi, col=col)
	Arrowhead(xy[1] - r/2, xy[2] + 2*w, angle=-90, arr.type="triangle", lcol=col)
	# TODO: other orientations;
}

# TODO:
# - other orientations;
# - better parametrization;
cvlarrow = function(xy1, xy2, col, w=w0, r=r0, mid.off = c(-0.035, 0), arr.length = 0.4) {
	user = par("usr")
    pcm  = par("pin") * 2.54
    sy = (user[4] - user[3])/pcm[2]; sx = 
	#
	rsq = r * sqrt(2); dtheta = 0.1; # may need dynamic adaptation;
	mid.off = mid.off - c(r + r/3, 0);
	dy = xy1[2] - xy2[2];
	if(dy > r) {
		xy = xy1 - c(rsq, rsq);
		r2 = sqrt(sum((xy2 - xy)^2)) + r/2; print(r2)
		filledcircle(r2 + w, r2, mid=xy2 + c(r2,0) + mid.off  + c(0.0015, 0),
			from=pi/2, to= -pi - dtheta, col=col)
		Arrowhead(xy2[1] + mid.off[1], xy2[2] + 2*w + mid.off[2], angle=-90,
			arr.type="triangle", lcol=col, arr.length=arr.length)
	} else if(dy < -r) {
		xy = xy1 # - c(rsq, rsq);
		r2 = sqrt(sum((xy2 - xy)^2)) - r/2; print(r2)
		filledcircle(r2 + w, r2, mid=xy2 + c(r2,0) + mid.off  + c(0.0015, 0),
			from= - pi + dtheta, to=-pi/2, col=col)
		Arrowhead(xy2[1] + mid.off[1], xy2[2] + 2*w + mid.off[2] - arr.length*sy*1.4, angle=90,
			arr.type="triangle", lcol=col, arr.length=arr.length)
	} else {
	}
}
circle = function(txt, xy, r=r0, col) {
	filledcircle(r, mid=xy, col=col)
	text(xy[1], xy[2], txt)
}
### double arrows
dlarrows = function(xy1, xy2, x.jt=0.005) {
	if(length(x.jt) < 2) x.jt = c(x.jt, 0);
	# TODO: all directions;
	# - see cases in lnarrow();
	lnarrow(xy1 + x.jt, xy2 + x.jt)
	lnarrow(xy2 - x.jt, xy1 - x.jt)
}

r0 = 0.1;
w0 = 0.01;


################

################
### Diagrams ###
################

x0 = 0.05;
y0 = 0.6;

### SIR + Hospital + Death

dev.new(width = 11.7, height = 8.3)
# png(file="Diagram.Model.H.png", width = 11.7, height = 8.3, units="in", res=100) # run this to save as png;
par.old = par(mar=c(0,0,2,0) + 0.01)
emptyplot(main = "SIR Model")

### Compartments
xyS = c(x0, y0);
circle("S", xyS, col="yellow");

xyI = c(x0 + 0.3, y0);
circle("I", xyI, col="red");

xyHs = c(x0, y0 - 3*r0);
circle("Hs", xyHs, col="grey");

xyHi = c(x0 + 0.3, y0 - 3*r0);
circle("Hi", xyHi, col="indianred1");

### R & D
xyR = c(x0 + 0.7, y0);
circle("R", xyR, col="green");

xyD = c(x0 + 0.7, y0 - 3*r0);
circle("D", xyD, col="indianred1");


### Arrows
lnarrow(xyS, xyI)
lnarrow(xyHs, xyHi)
# R & D
lnarrow(xyI, xyR)
lnarrow(xyI, xyD)
lnarrow(xyHi, xyR)
lnarrow(xyHi, xyD)

dlarrows(xyS, xyHs)
lnarrow(xyI, xyHi)


cvarrow(xyI, col="red")
cvarrow(xyHi, col="indianred1")


#######################
#######################

### SIR: Old age
# including Hospital + Death

### TODO: Extensions
# - H => (S -> I);
# - H => Hs + Hi;

Diagram2 = function(x0=0.02, y0=0.6, file="Diagram.Model.oldAge.png", save.png=FALSE) {
	if(save.png) {
		# run this to save as png;
		png(file=file, width = 11.7, height = 8.3, units="in", res=100)
	} else {
		dev.new(width = 11.7, height = 8.3)
	}
	par.old = par(mar=c(0,0,2,0) + 0.01)
	emptyplot(main = "SIR Model + Hospitalization")

	### Compartments

	### S & Old
	xyS = c(x0, y0);
	circle(expression(S[Y]), xyS, col="yellow");

	xyO = c(x0, y0 - 3*r0);
	circle(expression(S[Old]), xyO, col="grey");

	### Is & Iv
	xyIs = c(x0 + 0.3, y0);
	circle(expression(I[Y]), xyIs, col="red");

	xyIv = c(x0 + 0.3, y0 - 3*r0);
	circle(expression(I[Old]), xyIv, col="indianred1");

	### H
	# - only H due to Infection;
	xyH = c(x0 + 0.65, y0 + 0.25);
	circle("H", xyH, col="orange");

	### R & D
	xyR = c(x0 + 1, y0);
	circle("R", xyR, col="green");

	xyD = c(x0 + 1, y0 - 3*r0);
	circle("D", xyD, col="indianred1");


	### Arrows
	lnarrow(xyS, xyIs)
	lnarrow(xyO, xyIv)

	cvarrow(xyIs, col="red")
	cvarrow(xyIv, col="indianred1")
	cvlarrow(xyIs, xyIv, col="red")
	cvlarrow(xyIv, xyIs, col="indianred1")

	lnarrow(xyIs, xyR)
	lnarrow(xyIs, xyD)
	lnarrow(xyIv, xyR)
	lnarrow(xyIv, xyD)

	# Hospitalization
	lnarrow(xyIv, xyH)
	lnarrow(xyIs, xyH)
	# Outcome from H:
	lnarrow(xyH, xyR)
	lnarrow(xyH, xyD)
	if(save.png) {
		# close png
		Sys.sleep(0.1); # may be needed
		dev.off();
	}
}

x0 = 0.02;
y0 = 0.6;

Diagram2();


# png(file="Diagram.Model.oldAge.png", width = 11.7, height = 8.3, units="in", res=100)
	# ... run specific code: without dev.new()!
# dev.off() # close file



#######################
#######################

