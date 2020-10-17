


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
		p.x = matrix(rep( (0:(n-1)) * l/n, m), ncol=m)
	}
	p.y = sapply(1:m, ellipse.f, x=p.x)
	return(list("x"=p.x, "y"=p.y))
}

###########################

# install.packages("TSP")

library(TSP)


setwd("/Math")


###########################

###
ell.n = 7
p = watermelon.gen(20, m=ell.n, 10, random=F)
plot(p$x, p$y)

cities = matrix(c(as.vector(p$x), as.vector(p$y)), ncol=2)

etsp <- ETSP(cities)
etsp

### calculate a tour
tour <- solve_TSP(etsp, method = "nn")
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


