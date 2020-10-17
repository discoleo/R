########################
###
### Leonard Mada
### [the one and only]
###
### TSP Models
### Data Generators
###
### draft v.0.1b


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

find.base = function(m, middle=FALSE) {
	# m = matrix with coordinates of cities
	# only with non-random
	# TODO: more options
	if( ! middle) {
		isZero = m[,1] == 0 & m[,2] == 0
	} else {
		m0 = m[m[,2] == 0, 1]
		mid = m0[rank(m0) == round(length(m0)/2)]
		isZero = m[,2] == 0 & m[,1] == mid
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

