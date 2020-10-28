########################
###
### Leonard Mada
### [the one and only]
###
### TSP Models
### Data Generators
###
### draft v.0.3b


###############
### History ###

### draft v.0.3a - v.0.3b:
# - moved Analysis & Tours to TSP.Solver.R;
#   [basically all v.0.2x-code]
# - better formatting & comments (v.0.3b);
### draft v.0.2c - v.0.2d:
# - basic annealing algorithm;
# - initialisation using global medians as regularisation terms (v.0.2d);
# - TODO: explore also local neighbourhood medians;
### draft v.0.2b:
# - concept for algorithm:
#   add penalty for points with short median distance;
### draft v.0.2a:
# - started work on Analysis tools:
#  -- distance between locations on tour;
#  -- distance between all locations;
### draft v.0.1e - v.0.1f:
# - added concentric circles;
# - added variable number of points per circle (v.0.1f);
### draft v.0.1d:
# - added 2 overlapping densities (mixture of 2 gaussians);
### draft v.0.1c:
# - added overlapping circles;


####################
####################

### Generate Special Graphs for TSP
# Orbits resembling:
# 1.) Ellipses spanning the same Chord / Watermelon
# 2.) Intersecting Circles
# 3.) Concentric Circles
# 4.) Non-Orbit Densities

###############################

### various Generator functions

### Watermelon / Ellipses
watermelon.gen = function(n, m, l, random=TRUE, v.scale=3) {
	# n = number of points in a half-ellipse;
	# m = number of half-ellipses / half-orbits;
	# l = length/radius;
	r2 = l^2 / 4;
	ellipse.f = function(id, x) {
		a = 2*id - m - 1
		if(a == 0) {
			return (rep(0,n));
		}
		a = a/m
		x.x = abs(x[,id] - l/2)
		return(sign(a)*sqrt(r2 - x.x^2 + v.scale*(l/2 - x.x)*l*a^2))
	}
	if(random) {
		p.x = sapply(1:m, function(id) runif(n, 0, l))
	} else {
		# TODO: fix repeating (0, 0);
		p.x = matrix(rep( (1/2 + 0:(n-1)) * l/n, m), ncol=m)
		p.x[1, round(m/2)] = 0 # first Element;
	}
	p.y = sapply(1:m, ellipse.f, x=p.x)
	return(list("x"=p.x, "y"=p.y))
}
### Intersecting Circles
circles.int.gen = function(n, n.c, r=1) {
	# intersecting/overlapping circles
	central.th = (0:(n.c-1)) * 2*pi / n.c;
	circle.th = (0:(n-1)) * 2*pi / n;
	circle.f = function(theta, r, centre=c(0,0)) {
		x = r * cos(theta) + centre[1]	
		y = r * sin(theta) + centre[2]
		return(cbind(x, y))
	}
	circle.shift = function(id, base.circle, centre) {
		x = base.circle[1,] + centre[1, id];
		y = base.circle[2,] + centre[2, id];
		return(cbind(x, y));
	}
	# TODO: nicer code;
	circles.c = sapply(central.th, circle.f, r=r)
	circle.base = sapply(circle.th, circle.f, r=r)
	circles.pxy = sapply(1:ncol(circles.c), circle.shift, base.circle=circle.base, centre=circles.c)
	circles.pxy = list(
		x = as.vector(circles.pxy[1:n,]),
		y = as.vector(circles.pxy[-(1:n),]))
	# print(circles.pxy)
	return(circles.pxy)
}
### Concentric Circles
radial.gen = function(n, n.c, r=1, phi=0, d, addCenter=TRUE) {
	# concentric centers:
	# n = number of points/circle;
	# n.c = number of circles;
	# d = distance between points (variable points / circle);
	if(missing(n.c)) {
		if(length(r) > 1) {
			n.c = length(r)
		} else {
			stop("Stop: Must specify number of circles!")
		}
	}
	if(length(r) == 1) {
		r = r * 1:n.c;
	}
	if(missing(n)) {
		if(missing(d)) {
			stop("Stop: Number of points n or distance between points d must be provided!")
		}
		n = 2*pi*r / d;
		th = 2*pi/n
		n = round(n) # check if useful?
	} else {
		th = rep(2*pi/n, n.c)
		n = rep(n, n.c)
	}
	if(addCenter) {
		xy.df = data.frame(x=0, y=0)
	} else {
		xy.df = data.frame(x=numeric(0), y=numeric(0))
	}
	for(i in 1:n.c) {
		id = 0:(n[i]-1)
		x = r[i] * cos(th[i] * id + phi*i)
		y = r[i] * sin(th[i] * id + phi*i)
		xy.df = rbind(xy.df, cbind(x, y))
	}
	return(xy.df)
}
### 2D: N(0, sd1) x N(0, sd2)
rnorm2d.gen = function(n1, n2=n1, sd1=1, sd2=1) {
	len = n1*n2
	p.x = rnorm(len, sd=sd1)
	p.y = rnorm(len, sd=sd2)
	return(list("x"=p.x, "y"=p.y))
}
### 2D: U(0, max) x U(0, max)
runif2d.gen = function(n1, n2=n1, max=2) {
	# max = 2 is more close to rnorm
	# alternative: max = 4*sd
	len = n1*n2
	p.x = runif(len, max=max)
	p.y = runif(len, max=max)
	return(list("x"=p.x, "y"=p.y))
}
### 2D: 2 bands * epochs * N(0, sd)
r2d.gen = function(n, epochs, sd=1, sep.scale=1, x.jitter=NA) {
	# n = number of points per epoch
	s = rnorm(10*n, sd=sd)
	ep.id = if(length(epochs) == 1) 1:epochs else epochs;
	#
	y1 = as.vector(sapply(ep.id, function(id) sample(s, n))) + sd*sep.scale
	y2 = as.vector(sapply(ep.id, function(id) sample(s, n))) - sd*sep.scale
	ep.all = rep(ep.id, each=n)
	ep.all = c(ep.all, ep.all)
	if( ! is.na(x.jitter)) {
		ep.all = jitter(ep.all)
	}
	return(list("x"=ep.all, "y"=c(y1, y2)))
}

### Tools

### find specific location
find.base = function(m, y=0, middle=FALSE, type=c("gaussian", "regular")) {
	# m = matrix with coordinates of cities
	# y = y-limits to filter the city;
	# - currently works best/only with non-random data;
	# TODO:
	# - more options;
	# - regular vs gaussian vs uniform;
	
	# y-Coordinate
	if(length(y) == 2) {
		isZero = (m[,2] >= y[1]) & (m[,2] <= y[2])
	} else {
		isZero = m[,2] == y
	}
	#
	if( ! middle) {
		x.min = min(m[isZero, 1])
		isZero = isZero & m[,1] == x.min
	} else {
		m0 = m[isZero, 1]
		x.mid = m0[rank(m0) == round((length(m0) + 1)/2)]
		isZero = isZero & m[,1] == x.mid
	}
	id = match(TRUE, isZero)
	return(id)
}

write.tsp = function(x, file, asInt=TRUE, scale=1000) {
	if(asInt) {
		tsp.int = round(x * scale)
		if(inherits(x, "ETSP")) {
			tsp.int = ETSP(tsp.int)
			# still does NOT save as INTEGER!
			write_TSPLIB(tsp.int, file=file)
		} else {
			print("Nothing written: class NOT supported!")
		}
	} else {
		write_TSPLIB(x, file=file)
	}
}

### Analysis

### Distance between tour-locations
# best for regular lattices;
dist.tour = function(cities, tour) {
	dist.f = function(c1, c2) {
		return(sqrt( sum((c1 - c2)^2) ))
	}
	dist.byid = function(id) {
		dist.f(cities[tour[tour.id[id]], ], cities[tour[tour.id[id + 1]], ])
	}
	tour.id = c(seq_along(tour), 1)
	return(sapply(seq_along(tour), dist.byid))
}

dist.all = function(cities) {
	dist.byid = function(id) {
		return(sqrt( ((cities[id,1] - cities[,1])^2 + (cities[id,2] - cities[,2])^2) ))
	}
	id = 1:nrow(cities)
	return(sapply(id, dist.byid))
}
### Median distance to all other points
dist.med = function(cities) {
	d = dist.all(cities)
	return(sapply(1:nrow(d), function(id) median(d[id, ])))
}

###########################

# install.packages("TSP")

library(TSP)

### 3D
library(rgl)


setwd("/Math")


###########################
###########################

#######################
### Test Generators ###

###
ell.n = 7
p = watermelon.gen(20, m=ell.n, 10, random=F)
plot(p$x, p$y)

cities = matrix(c(as.vector(p$x), as.vector(p$y)), ncol=2)
id = find.base(cities, middle=T)
id

### TSP data
etsp <- ETSP(cities)
etsp

### calculate a tour
# tour <- solve_TSP(etsp, method = "nn"
tour <- solve_TSP(etsp, method = "nn", control=list(start=id))
tour

tour_length(tour)
plot(etsp, tour, tour_col = "red")
points(cities[id,1], cities[id,2], col="green")


### Save tour as image
SAVE_PNG = TRUE
if(SAVE_PNG) {
png(file="TSP.Watermelon.png")
	plot(etsp, tour, tour_col = "red", xlab="X-Coord", ylab="Y-Coord")
	points(cities[id,1], cities[id,2], col="green")
dev.off()
}

### 3D
library(rgl)
#
plot3d(p$x, p$y, p$x^2 + p$y^2, type="s", col=1/2*abs(p$y)+1)
val.r = (p$x - max(p$x)/2)^2 + (abs(p$y) - max(p$y)/2)^2
plot3d(p$x, p$y, val.r, type="s", col=1/2*abs(p$y)+1)
plot3d(p$x, p$y, sqrt(val.r), type="s", col=1/2*abs(p$y)+1)

### Save coordinates as TSP file
# as Integer:
write.tsp(etsp, file="Watermelon.tsp", asInt=TRUE)

# as initial data
write_TSPLIB(etsp, file="Watermelon.tsp")


### Other: for ATSP
# image(atsp, tour)


###########################

###
p = rnorm2d.gen(10)
# p = runif2d.gen(10)
plot(p$x, p$y)


cities = matrix(c(as.vector(p$x), as.vector(p$y)), ncol=2)

etsp <- ETSP(cities)
etsp

### calculate a tour
tour <- solve_TSP(etsp, method = "nn")
tour

tour_length(tour)
plot(etsp, tour, tour_col = "red")


#########################

### Overlapping Circles
p = circles.int.gen(17, 5, 2)
plot(p$x, p$y)


cities = matrix(c(as.vector(p$x), as.vector(p$y)), ncol=2)

etsp <- ETSP(cities)
etsp

### calculate a tour
tour <- solve_TSP(etsp, method = "nn")
tour

tour_length(tour)
plot(etsp, tour, tour_col = "red")


#######################

### Specific Densities
p = r2d.gen(5, 10, sd=2, sep.scale=3, x.jitter=T)
plot(p$x, p$y)

### Q:
# Are there any phase transitions determined by sd & sep.scale?
# Are there various phase transitions when sep.scale changes:
# between ~1.5 and ~2 and between ~2 and ~3?
# How can we measure and describe phase transitions?


cities = matrix(c(as.vector(p$x), as.vector(p$y)), ncol=2)

etsp <- ETSP(cities)
etsp

### calculate a tour
tour <- solve_TSP(etsp, method = "nn")
tour

tour_length(tour)
plot(etsp, tour, tour_col = "red")


### with Base-city:
id = find.base(cities, y=c(0, 2.5), middle=T)
id
#
tour <- solve_TSP(etsp, method = "nn", control=list(start=id))
tour

tour_length(tour)
plot(etsp, tour, tour_col = "red")
points(cities[id,1], cities[id,2], col="green")


##########################
##########################

##########################
### Concentric Circles ###

phi = 0.3 # 0.2 # 0;
p = radial.gen(17, 5, phi=phi)
plot(p$x, p$y)


cities = matrix(c(as.vector(p$x), as.vector(p$y)), ncol=2)

### with Base-city:
# id = find.base(cities, y=0, middle=F)
# id = find.base(cities, y=c(2, 5), middle=T)
id = find.base(cities, y=c(0.5, 2), middle=T)
id

etsp <- ETSP(cities)
etsp

### calculate a tour
tour <- solve_TSP(etsp, method = "nn", control=list(start=id))
tour

tour_length(tour)
plot(etsp, tour, tour_col = "red", xlab="X-Coord", ylab="Y-Coord")
points(cities[id,1], cities[id,2], col="green")


### 3D
plot3d(p$x, p$y, p$x^2 + p$y^2, type="s", col=1/2*abs(p$y)+1)
#
val.r = p$x^2 + p$y^2
plot3d(p$x, p$y, val.r, type="s", col = val.r/2 + 1)
plot3d(p$x, p$y, sqrt(val.r), type="s", col = sqrt(val.r)+1)


###########

### variable number of points / circle
phi = -0.01 # 0.2 # 0;
p = radial.gen(d=1 + (1:7)/17, n.c=7, phi=phi)
plot(p$x, p$y)


cities = matrix(c(as.vector(p$x), as.vector(p$y)), ncol=2)

### with Base-city:
# id = find.base(cities, y=0, middle=F)
# id = find.base(cities, y=c(2, 5), middle=T)
id = find.base(cities, y=c(0.5, 2), middle=T)
id

etsp <- ETSP(cities)
etsp

### calculate a tour
tour <- solve_TSP(etsp, method = "nn", control=list(start=id))
tour

tour_length(tour)
plot(etsp, tour, tour_col = "red", xlab="X-Coord", ylab="Y-Coord")
points(cities[id,1], cities[id,2], col="green")


# write.tsp(etsp, file="TSP.Circular.tsp", asInt=TRUE)

#######################
#######################
