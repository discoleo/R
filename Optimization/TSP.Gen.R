########################
###
### Leonard Mada
### [the one and only]
###
### TSP Models
### Data Generators
###
### draft v.0.1d


###############
### History ###

### draft v.0.1d:
# - added 2 overlapping densities (mixture of 2 gaussians);
### draft v.0.1c:
# - added overlapping circles;


####################
####################

### Generate Special Graphs for TSP

watermelon.gen = function(n, m, l, random=TRUE, d=3) {
	# n = number of points in a half-ellipse;
	# m = number of half-ellipses;
	r2 = l^2 / 4 # * m
	ellipse.f = function(id, x) {
		a = 2*id - m - 1
		if(a == 0) {
			return (rep(0,n));
		}
		a = a/m
		x.x = abs(x[,id] - l/2)
		return(sign(a)*sqrt(r2 - x.x^2 + d*(l/2 - x.x)*l*a^2))
	}
	if(random) {
		p.x = sapply(1:m, function(id) runif(n, 0, l))
	} else {
		# TODO: fix repeating (0, 0);
		p.x = matrix(rep( (0:(n-1)) * l/n, m), ncol=m)
	}
	p.y = sapply(1:m, ellipse.f, x=p.x)
	return(list("x"=p.x, "y"=p.y))
}

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

rnorm2d.gen = function(n1, n2=n1, sd1=1, sd2=1) {
	len = n1*n2
	p.x = rnorm(len, sd=sd1)
	p.y = rnorm(len, sd=sd2)
	return(list("x"=p.x, "y"=p.y))
}

runif2d.gen = function(n1, n2=n1, max=2) {
	# max = 2 is more close to rnorm
	# alternative: max = 4*sd
	len = n1*n2
	p.x = runif(len, max=max)
	p.y = runif(len, max=max)
	return(list("x"=p.x, "y"=p.y))
}

r2d.gen = function(n, ep, sd=1, sep.scale=1) {
	# ep = epochs
	s = rnorm(10*n, sd=sd)
	ep.id = 1:ep
	y1 = as.vector(sapply(ep.id, function(id) sample(s, n))) + sd*sep.scale
	y2 = as.vector(sapply(ep.id, function(id) sample(s, n))) - sd*sep.scale
	ep.all = rep(ep.id, each=n)
	return(list("x"=c(ep.all, ep.all), "y"=c(y1, y2)))
}

find.base = function(m, y=0, middle=FALSE) {
	# m = matrix with coordinates of cities
	# only with non-random
	# TODO: more options
	if( ! middle) {
		isZero = m[,1] == y & m[,2] == y
	} else {
		isZero = m[,2] == y
		m0 = m[isZero, 1]
		mid = m0[rank(m0) == round(length(m0)/2)]
		isZero = isZero & m[,1] == mid
	}
	id = match(TRUE, isZero)
	return(id)
}

###########################

# install.packages("TSP")

library(TSP)


setwd("C:/Users/Leo Mada/Desktop/DB/Math")


###########################

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


SAVE_PNG = TRUE
if(SAVE_PNG) {
png(file="TSP.Watermelon.png")
	plot(etsp, tour, tour_col = "red")
dev.off()
}

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
p = r2d.gen(5, 10, sd=2, sep.scale=1.5)
plot(p$x, p$y)

# Q: Are there any phase transitions determined by sd & sep.scale?


cities = matrix(c(as.vector(p$x), as.vector(p$y)), ncol=2)

etsp <- ETSP(cities)
etsp

### calculate a tour
tour <- solve_TSP(etsp, method = "nn")
tour

tour_length(tour)
plot(etsp, tour, tour_col = "red")


#######################

### Analysis
# TODO

# 1.) Phase Transitions in the data
# - How to measure phase transitions?
# 2.) Invariants, pseudo-Invariants
# - Higher Moments;
# - Higher "Moments" of Correlation;
# - other pseudo-invariants;

