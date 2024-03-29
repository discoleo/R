###
### West University:
### Team Project 2021
###
### Team: Modeling Infection Spread
### Supervisor: Leonard Mada / Syonic
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
circle = function(txt, xy, r=r0, col, cex) {
	filledcircle(r, mid=xy, col=col)
	if(missing(cex)) cex = par("cex");
	text(xy[1], xy[2], txt, cex=cex)
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

### Simple SIR

diagram.SIR = function(file=NULL, cex=2, dim=c(11.7, 8.3), xy=c(0.05, 0.6), main="SIR Model") {
	if(is.null(file)) {
		dev.new(width = dim[1], height = dim[2])
	} else {
		png(file=file, width = dim[1], height = dim[2], units="in", res=100);
	}
	
	par.old = par(mar=c(0,0,2,0) + 0.01);
	emptyplot(main = main);
	x0 = xy[1]; y0 = xy[2];

	### Compartments
	
	# S
	xyS = c(x0, y0);
	circle("S", xyS, col="yellow", cex=cex);
	
	# I
	xyI = c(x0 + 0.35, y0);
	circle("I", xyI, col="red", cex=cex);
	
	# R
	xyR = c(x0 + 0.7, y0);
	circle("R", xyR, col="green", cex=cex);


	### Arrows
	lnarrow(xyS, xyI)
	# R
	lnarrow(xyI, xyR)
	# I -> S
	cvarrow(xyI, col="red")
	
	if( ! is.null(file)) dev.off();
	return(invisible());
}

diagram.SIR()


# diagram.SIR(file = "Diagram.Model.Simple.png", dim=c(8, 4))


################
################

### Extended SIR
### + Hospital + Death

diagram.ExtHD = function(xy=c(0.05, 0.6), cex=2, file=NULL, dim=c(11.7, 8.3),
		main="SIR Model\n+ Hospitalization & Death", save.png=FALSE) {
	if(is.null(file) && ! save.png) {
		dev.new(width = dim[1], height = dim[2]);
	} else {
		if(is.null(file)) file = "Diagram.Model.H.png";
		png(file=file, width = dim[1], height = dim[2], units="in", res=100);
	}
	
	par.old = par(mar=c(0,0,2,0) + 0.01)
	emptyplot(main = main)
	x0 = xy[1]; y0 = xy[2];
	
	### Compartments
	xyS = c(x0, y0);
	circle("S", xyS, col="yellow", cex=cex);
	
	xyI = c(x0 + 0.35, y0);
	circle("I", xyI, col="red", cex=cex);
	
	xyHs = c(x0, y0 - 3*r0);
	circle(expression(H[S]), xyHs, col="grey", cex=cex);
	
	xyHi = c(x0 + 0.35, y0 - 3*r0);
	circle(expression(H[I]), xyHi, col="indianred1", cex=cex);

	### R & D
	xyR = c(x0 + 0.7, y0);
	circle("R", xyR, col="green", cex=cex);

	xyD = c(x0 + 0.7, y0 - 3*r0);
	circle("D", xyD, col="indianred1", cex=cex);


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
	# Cross-Category
	cvlarrow(xyHi, xyI, col="indianred1")
	cvlarrow(xyI, xyHi, col="red")
	
	if( ! is.null(file)) dev.off();
	return(invisible());
}

diagram.ExtHD()

# diagram.ExtHD(save.png=TRUE)


#######################
#######################

### SIR: Old age
# including Hospital + Death

### TODO: Extensions
# - H => (S -> I);
# - H => Hs + Hi;

diagram.ExtAge = function(xy = c(0.02, 0.6), cex=2, file="Diagram.Model.oldAge.png",
		main="SIR Model\nOld Age", save.png=FALSE) {
	if(save.png) {
		png(file=file, width = 11.7, height = 8.3, units="in", res=100)
	} else {
		dev.new(width = 11.7, height = 8.3)
	}
	
	x0 = xy[1]; y0 = xy[2];
	par.old = par(mar=c(0,0,2,0) + 0.01)
	emptyplot(main=main)

	### Compartments

	### S & Old
	xyS = c(x0, y0);
	circle(expression(S[Y]), xyS, col="yellow", cex=cex);

	xyO = c(x0, y0 - 3*r0);
	circle(expression(S[Old]), xyO, col="grey", cex=cex);

	### Is & Iv
	xyIs = c(x0 + 0.35, y0);
	circle(expression(I[Y]), xyIs, col="red", cex=cex);

	xyIo = c(x0 + 0.35, y0 - 3*r0);
	circle(expression(I[Old]), xyIo, col="indianred1", cex=cex);

	### H
	# - only H due to Infection;
	xyH = c(x0 + 0.65, y0 + 0.25);
	circle("H", xyH, col="orange", cex=cex);

	### R & D
	xyR = c(x0 + 1, y0);
	circle("R", xyR, col="green", cex=cex);

	xyD = c(x0 + 1, y0 - 3*r0);
	circle("D", xyD, col="indianred1", cex=cex);


	### Arrows
	lnarrow(xyS, xyIs)
	lnarrow(xyO, xyIo)

	cvarrow(xyIs, col="red")
	cvarrow(xyIo, col="indianred1")
	# Cross-Category
	cvlarrow(xyIs, xyIo, col="red")
	cvlarrow(xyIo, xyIs, col="indianred1")

	lnarrow(xyIs, xyR)
	lnarrow(xyIs, xyD)
	lnarrow(xyIo, xyR)
	lnarrow(xyIo, xyD)

	# Hospitalization
	lnarrow(xyIo, xyH)
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


diagram.ExtAge();


# diagram.ExtAge(save.png=TRUE);


#######################
#######################

