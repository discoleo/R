########################
###
### Leonard Mada
### [the one and only]
###
### TSP Solver
###
### draft v.0.1g


###############
### History ###

### draft v.0.1e - v.0.1g:
# - early experiments with multiple (generalized) polynomial terms (0.1e & 0.1f & 0.1g);
# - early experiments with bands of data split into "epochs" (0.1g);
# - new idea: density of locations in a specified neighbourhood (v.0.1f);
### draft v.0.1b - v.0.1d:
# - some experiments with median-based regularisation;
# - more experiments & explorations: concentric circles need local optimizations;
# - bug fix: in tour length (v.0.1d);
### draft v.0.1a:
# - moved from TSP.Gen.R;
### from TSP.Gen.R [v.0.2c - v.0.2d]:
# - basic annealing algorithm;
# - initialisation using global medians as regularisation terms (v.0.2d);
# - TODO: explore also local neighbourhood medians;


##########################

### helper functions

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

########
### Tour
init.tour = function(cities, start=1, alpha=1/3, p=1, dist.m) {
	len = nrow(cities)
	remain.id = seq(len)
	if(missing(dist.m)) {
		d = dist.all(cities)
		d.med = sapply(1:nrow(d), function(id) median(d[id, ]))
	} else if( ! is.null(dim(dist.m))) {
		d = dist.m
		d.med = sapply(1:nrow(d), function(id) median(d[id, ]))
	} else {
		d.med = dist.m
	}
	d = d.med
	#
	dist.reg = function(id, current.city, remain.id, d, current.max) {
		sqrt(sum((cities[remain.id[id],] - current.city)^2)) +
		alpha * (current.max^p - d[id]^p)
	}
	# multiple coefficients alpha & p
	dist.mreg = function(id, current.city, remain.id, d, current.max) {
		sqrt(sum((cities[remain.id[id],] - current.city)^2)) +
		sum(alpha * (current.max^p - d[id]^p))
	}
	if(length(alpha) > 1 || length(p) > 1) {
		dist.reg = dist.mreg;
	}
	# else d = d^p; p = 1; # optimization
	#
	tour = remain.id[start]
	current.id = start
	iter = rep(0, len-2)
	for(i in iter) {
		current.city = cities[remain.id[current.id], ]
		# remove current city
		remain.id = remain.id[ - current.id]
		d = d[ - current.id]
		# cost of next city
		current.max = max(d);
		cost = sapply(seq(remain.id), dist.reg, current.city, remain.id, d, current.max)
		cost.min = min(cost)
		next.id = match(cost.min, cost)
		# update
		current.id = next.id
		tour = c(tour, remain.id[current.id])
	}
	tour = c(tour, remain.id[-current.id])
	return(tour)
}

### basic Annealing Algorithm:
optim.tour = function(cities, tour=NA, T_Max=1000, T_Min=1E-4) {
	accept = function(current_cost, new_cost, T) {
		if (new_cost <= current_cost) {
			return(TRUE);
		}
		if (runif(1, 0, 1) < exp(-(new_cost - current_cost)/T)) {
			return(TRUE);
		} else return(FALSE);
	}
	update.Temp = function(T, k) {
		return(T * 0.99995)
	}
	eval.tour = function(cities, tour) {
        return (sum(dist.tour(cities, tour)))
	}
	perturb.tour = function(ids, S) {
		id = sample(ids, 2) # i, j
        if (id[1] > id[2]){ id = id[c(2, 1)]; }
		# S is already a copy of the original;
        for(k in seq((id[2] - id[1]) %/% 2)) {
			S[c(id[1]+k, id[2]-k)] = S[c(id[2]-k, id[1]+k)]
		}
        return(S);
	}
	init.tour = function(cities) {
		# naive initialisation
		len = nrow(cities)
		return(sample(seq(len), len));
	}
	
	##########
	# SA(prob, T_Max, T_Min):
	if(is.na(tour)) {
		S = init.tour(cities);
	} else {
		S = tour;
	}
    S_cost = eval.tour(cities, S);
    
    S_best = S
    S_best_cost = S_cost
    
    T = T_Max
    k = 0
	idAll = 1:nrow(cities)
    while (T > T_Min) {
        k = k+1
        S_prim = perturb.tour(idAll, S)
        S_prim_cost = eval.tour(cities, S_prim);
        
        if (accept(S_cost, S_prim_cost, T)) {
            S = S_prim;
            S_cost = S_prim_cost
		}
        if (S_cost < S_best_cost) {
            S_best = S;
            S_best_cost = S_cost
            print(paste(S_best_cost, T))
		}
            
        T = update.Temp(T, k);
	}
    
    return (list("cost"=S_best_cost, "S"=S_best))
}

scan.tour = function(cities, start.id=1, alpha=(-5:5)/3, p=(-5:5)/3, param.grid, dist.m, ignore.start=FALSE) {
	### scan parameter-space
	if(missing(dist.m)) {
		dist.m = dist.med(cities)
	}
	if(missing(param.grid)) {
		param.grid = expand.grid(alpha, p)
	}
	#
	tours = sapply(seq(nrow(param.grid)),
		function(id) init.tour(cities, start=if(ignore.start) 1 else start.id, dist.m=dist.m,
			alpha=param.grid[id,1], p=param.grid[id,2]))
	tours.d = sapply(seq(ncol(tours)), function(id) sum(dist.tour(cities, tours[,id])))
	d.min = min(tours.d)
	tour.id = seq_along(tours.d) [tours.d == d.min]
	rez = list(min=d.min, id=tour.id, alpha=param.grid[tour.id,1], p=param.grid[tour.id,2])
	plot.tour(tours[ , tour.id[1]], cities, etsp, start.id)
	return(rez)
}
scan2.tour = function(cities, start.id=1, alpha=2/3, p=2/3, param.grid, dist.m, ignore.start=FALSE) {
	### scan parameter-space: Level 2
	if(missing(dist.m)) {
		dist.m = dist.med(cities)
	}
	if(missing(param.grid)) {
		param.grid = expand.grid((-5:5)/3, (-5:5)/3)
	}
	#
	tours = sapply(seq(nrow(param.grid)),
		function(id) init.tour(cities, start=if(ignore.start) 1 else start.id, dist.m=dist.m,
			alpha=c(alpha, param.grid[id,1]), p=c(p, param.grid[id,2])))
	tours.d = sapply(seq(ncol(tours)), function(id) sum(dist.tour(cities, tours[,id])))
	d.min = min(tours.d)
	tour.id = seq_along(tours.d) [tours.d == d.min]
	rez = list(min=d.min, id=tour.id, alpha=param.grid[tour.id,1], p=param.grid[tour.id,2])
	plot.tour(tours[ , tour.id[1]], cities, etsp, start.id)
	return(rez)
}

#############################
#############################

################
### Analysis ###

# 1.) Phase Transitions in the data
# - How to measure phase transitions?
# 2.) Invariants, pseudo-Invariants
# - Higher Moments;
# - Higher "Moments" of Correlation;
# - non-linear correlation, x-"autocorrelation" or y-"autocorrelation";
# - "divergence", "curl";
# - separation in higher dimensions;
# - other pseudo-invariants;
# 2.) Data
# - compressible/incompressible vectors/regions/clusters;

### TODO:
# - proper concepts of analysis;


#################
### Algorithm ###

### TODO:
# - penalty term based on Median distance:
#   Cost = Distance() +
#         (Max(median_distances of remaining points) - median_distance(current_point));
#   Cost = Distance() + alpha * Diff(Max(remaining_median)^p, Median(current)^p),
#   where alpha = coefficient, p = some power (e.g. p = 1/2; may use an L2 norm as well);
# - penalize points that have short median distance;
# - force to include first points that are poorly connected (large median distance);
# - once a point is added to the tour, remove its median value
#   from the list of remaining points;
# - it is possible to use regularisation terms both for the global & the local medians;
#  -- global = all remaining locations, not yet included in the tour;
#  -- local = all locations in a local neighbourhood, not yet included in the tour;


################
################

# uses functions from TSP.Gen.R;

test.Tour.gen = function(n=4, phi=-0.01, shift.scale=1/17,
			start.city=list(y=c(0.5, 2), middle=TRUE), printXY=TRUE) {
	# phi = -0.01 # 0.2 # 0;
	p = radial.gen(d = 1 + (1:n)*shift.scale, r = 5 + (1:n)/3, phi=phi, addCenter=F)
	if(printXY) {
		plot(p$x, p$y)
	}

	cities = matrix(c(as.vector(p$x), as.vector(p$y)), ncol=2)

	### with Base-city:
	id = find.base(cities, y=start.city$y, middle=start.city$middle)
	print(id)

	etsp = ETSP(cities)
	return(list(cities=cities, etsp=etsp, id=id))
}

plot.tour = function(tour, cities, etsp, id, col="red", start.col="green", digits=3) {
	plot(etsp, tour, tour_col = col, xlab="X-Coord", ylab="Y-Coord")
	points(cities[id,1], cities[id,2], col=start.col)
	d = sum(dist.tour(cities, tour))
	mtext(paste0("D = ", round(d, digits)), 3, line=1, adj=0.5)
	return(d)
}

################

################
### Examples ###


##################
### Examples: Tour
data = test.Tour.gen(4)
cities = data$cities; etsp = data$etsp; id = data$id;

### Tour
# alpha = 15; p = 1/5
alpha = 0; # 95.65
alpha = -1/2; p = 0.9; # 82.2
alpha = 1/5; p = 1; # 78.55
alpha = 2/3; p = 1; # 77.747
tour = init.tour(cities, alpha=alpha, p=p)
plot.tour(tour, cities, etsp, id)


# Note: takes long!!!
# and may NOT improve the tour!
Temp = 200; # 1000 # lets start lower
tour = optim.tour(cities, tour=tour, Temp, 1E-4)
tour; tour = tour$S
plot.tour(tour, cities, etsp, id)
# starting with lower Temp and the 82.2 tour:
# => 72.64! [but takes quit a while!]
# [another tour with 71.95 is also known]
tour = c(1, 33, 65, 97, 32, 64, 96, 128, 31, 63, 95, 127, 30, 62, 94, 126, 29, 61, 93, 125, 28, 60, 92, 124, 27, 59,
	91, 123, 26, 58, 90, 122, 25, 57, 89, 121, 24, 56, 88, 120, 23, 55, 87, 119, 22, 54, 86, 118, 117, 85, 53,
	21, 20, 52, 84, 116, 115, 83, 51, 19, 18, 50, 82, 114, 113, 81, 49, 17, 16, 48, 80, 112, 15, 47, 79, 111,
	110, 78, 46, 14, 13, 45, 77, 109, 108, 76, 44, 12, 11, 43, 75, 107, 106, 74, 42, 10, 9, 41, 73, 105, 104,
	72, 40, 8, 7, 39, 71, 103, 6, 38, 70, 102, 101, 69, 37, 5, 4, 36, 68, 100, 99, 67, 35, 3, 2, 34, 66, 98)

# more analysis:
# see Section at line ~ 400;


##################

##################
### Examples: Tour
data = test.Tour.gen(5)
cities = data$cities; etsp = data$etsp; id = data$id;

### Tour
alpha = 1/5; p = 1;
alpha = 1.1; p = 1/4;
alpha = 0; # 110.73
alpha = -0.5; p =1/4; # 88.96
alpha = -2/3; p =1/5; # 88.96
tour = init.tour(cities, alpha=alpha, p=p)
plot.tour(tour, cities, etsp, id)



#############################
#############################

### Examples: Analysis

### Example 1:
### variable number of points / circle
phi = -0.01 # 0.2 # 0;
cities.p = radial.gen(d=1 + (1:7)/17, n.c=7, phi=phi)
plot(cities.p$x, cities.p$y)

cities = matrix(c(as.vector(cities.p$x), as.vector(cities.p$y)), ncol=2)

### with Base-city:
id = find.base(cities, y=c(0.5, 2), middle=T)
id

etsp <- ETSP(cities)
etsp

### calculate a tour
tour <- solve_TSP(etsp, method = "nn", control=list(start=id))
tour
tour_length(tour)

plot.tour(tour, cities, etsp, id) # 186.986, 177.465

alpha = 0; p = 1; # 178.773
alpha = -2/3; p = 1/5; # 175.9
alpha = 1/3; p = 1; # 169.798
tour = init.tour(cities, alpha=alpha, p=p)
plot.tour(tour, cities, etsp, id)


### SA
# Note: takes long!!!
# and may NOT improve the tour!
Temp = 200; # 180; # 1000 # lets start lower
tour = optim.tour(cities, tour=tour, Temp, 1E-4)
tour; tour = tour$S
plot.tour(tour, cities, etsp, id)
# starting with lower Temp and the 169.8 tour:
# => 163.898 (with T=200) & 165.04 (with T=180)! [but takes quit a while!]
# [Note: this is a random process!]


### scan parameter-space
dist.m = dist.med(cities)
param.grid = expand.grid((-5:5)/3, (-5:5)/3)
tours = sapply(seq(nrow(param.grid)),
	function(id) init.tour(cities, dist.m=dist.m, alpha=c(param.grid[id,1], -1/5), p=c(param.grid[id,2], 1)))
tours.d = sapply(seq(ncol(tours)), function(id) sum(dist.tour(cities, tours[,id])))
d.min = min(tours.d)
tour.id = seq_along(tours.d) [tours.d == d.min]
list(min=d.min, id=tour.id, alpha=param.grid[tour.id,1], p=param.grid[tour.id,2])
plot.tour(tours[ , tour.id[1]], cities, etsp, id)
# 166.65
tour = init.tour(cities, alpha=c(5/3, -1/5), p=c(2/3, 1))
plot.tour(tour, cities, etsp, id)


### Analysis ###

### all Locations
d = dist.all(cities)
# d.v = d.m[d.m != 0]
# summary(d.v)

d.min = sapply(1:nrow(d), function(id) min(d[id, d[id,] != 0]))
d.med = sapply(1:nrow(d), function(id) median(d[id, ]))
d.len = sapply(1:nrow(d), function(id) sum(d[id, ] < 2.5))
summary(d.min)
summary(d.med)
summary(d.len)
plot3d(cities.p$x, cities.p$y, d.min, type="s", col = round(10*d.min) + 1)
# nice 3D decomposition!
plot3d(cities.p$x, cities.p$y, d.med, type="s", col = round(10*d.min) + 1)
# less distortion
plot3d(cities.p$x, cities.p$y, (d.med)^(1/3), type="s", col = round(d.min + d.med) + 1)
# density of local points
plot3d(cities.p$x, cities.p$y, d.len, type="s", col = round(d.min + d.med) + 1)
#
plot3d(cities.p$x, cities.p$y, d.med, type="s", col = d.len - min(d.len) + 1)


### Tour
d.tr = dist.tour(cities, tour)
summary(d.tr)

len = tour_length(tour)
id.col = (1:len) / len
plot3d(p$x[tour], p$y[tour], d.med[tour], type="s", col = rgb(id.col, 1-id.col, 0))


#######################
#######################

##############
### Example 2:
### variable number of points / circle
phi = -0.01 # 0.2 # 0;
cities.p = radial.gen(d=1 + (1:4)/17, r= 5 + (1:4)/3, phi=phi, addCenter=F)
plot(cities.p$x, cities.p$y)

cities = matrix(c(as.vector(cities.p$x), as.vector(cities.p$y)), ncol=2)

### with Base-city:
id = find.base(cities, y=c(0.5, 2), middle=T)
id

etsp <- ETSP(cities)
etsp

### calculate a tour
tour <- solve_TSP(etsp, method = "nn", control=list(start=id))
tour # 102.14

tour_length(tour)
plot.tour(tour, cities, etsp, id)

################
### Analysis ###

### all Locations
d = dist.all(cities)
# d.v = d.m[d.m != 0]
# summary(d.v)

d.min = sapply(1:nrow(d), function(id) min(d[id, d[id,] != 0]))
d.med = sapply(1:nrow(d), function(id) median(d[id, ]))
d.len = sapply(1:nrow(d), function(id) sum(d[id, ] < 2.5))
summary(d.min)
summary(d.med)
summary(d.len)
plot3d(cities.p$x, cities.p$y, d.min, type="s", col = round(10*d.min) + 1)
# nice 3D decomposition!
plot3d(cities.p$x, cities.p$y, d.med, type="s", col = round(10*d.min) + 1)
# less distortion
plot3d(cities.p$x, cities.p$y, (d.med)^(1/3), type="s", col = round(d.min + d.med) + 1)
# density of local points
plot3d(cities.p$x, cities.p$y, d.len, type="s", col = round(d.min + d.med) + 1)
#
plot3d(cities.p$x, cities.p$y, d.med, type="s", col = d.len - min(d.len) + 1)


### Tour ###
### scan parameter-space
dist.m = dist.med(cities)
param.grid = expand.grid((-5:5)/3, (-5:5)/3)
tours = sapply(seq(nrow(param.grid)),
	function(id) init.tour(cities, dist.m=dist.m, alpha=param.grid[id,1], p=param.grid[id,2]))
tours.d = sapply(seq(ncol(tours)), function(id) sum(dist.tour(cities, tours[,id])))
d.min = min(tours.d)
tour.id = seq_along(tours.d) [tours.d == d.min]
list(min=d.min, id=tour.id, alpha=param.grid[tour.id,1], p=param.grid[tour.id,2])
plot.tour(tours[ , tour.id[1]], cities, etsp, id)
# 77.747
tour = init.tour(cities, alpha=2/3, p=1)
plot.tour(tour, cities, etsp, id)

###
tours = sapply(seq(nrow(param.grid)),
	function(id) init.tour(cities, dist.m=dist.m, alpha=c(2/3, param.grid[id,1]), p=c(1, param.grid[id,2])))
tours.d = sapply(seq(ncol(tours)), function(id) sum(dist.tour(cities, tours[,id])))
d.min = min(tours.d)
tour.id = seq_along(tours.d) [tours.d == d.min]
list(min=d.min, id=tour.id, alpha=param.grid[tour.id,1], p=param.grid[tour.id,2])
plot.tour(tours[ , tour.id[1]], cities, etsp, id)
# 77.55 ### minimal improvement
tour = init.tour(cities, alpha=c(2/3, -1), p=c(1, -5/3))
plot.tour(tour, cities, etsp, id)

### SA
# Note: takes long!!!
# and may NOT improve the tour!
Temp = 100; # 180; # lets start lower
tour = optim.tour(cities, tour=tour, Temp, 1E-4)
tour; tour = tour$S
plot.tour(tour, cities, etsp, id)
# starting with lower Temp and the 77.55 tour:
# => 71.9509! [but takes quit a while!]
tour = c(1, 33, 65, 97, 98, 66, 34, 2, 3, 35, 67, 99, 100, 68, 36, 4, 5, 37, 69, 101, 102,
	70, 38, 6, 7, 39, 71, 103, 104, 72, 40, 8, 9, 41, 73, 105, 106, 74, 42, 10, 107, 75, 43,
	11, 12, 44, 76, 108, 109, 77, 45, 13, 14, 46, 78, 110, 111, 79, 47, 15, 16, 48, 80, 112,
	113, 81, 49, 17, 114, 82, 50, 18, 19, 51, 83, 115, 116, 84, 52, 20, 21, 53, 85, 117, 118,
	86, 54, 22, 119, 87, 55, 23, 120, 88, 56, 24, 121, 89, 57, 25, 122, 90, 58, 26, 123, 91,
	59, 27, 124, 92, 60, 28, 125, 93, 61, 29, 126, 94, 62, 30, 127, 95, 63, 31, 128, 96, 64, 32)
plot.tour(tour, cities, etsp, id)


##################
##################

##################
### Watermelon ###

ell.n = 7
cities.p = watermelon.gen(20, m=ell.n, 10, random=F)
plot(cities.p$x, cities.p$y)

cities = matrix(c(as.vector(cities.p$x), as.vector(cities.p$y)), ncol=2)
id = find.base(cities, middle=T)
id

### TSP data
etsp <- ETSP(cities)
etsp

### calculate a tour
# tour <- solve_TSP(etsp, method = "nn")
tour <- solve_TSP(etsp, method = "nn", control=list(start=id))
tour

tour_length(tour) # 132.2
plot.tour(tour, cities, etsp, id)

### Tour ###
### scan parameter-space
dist.m = dist.med(cities)
param.grid = expand.grid((-5:5)/3, (-5:5)/3)
#
scan.tour(cities, id, param.grid=param.grid, ignore.start=T)
# 112.57
alpha = 2/3; p = 2/3;
tour = init.tour(cities, alpha=alpha, p=p)
plot.tour(tour, cities, etsp, id)

### Level 2:
tours = sapply(seq(nrow(param.grid)),
	function(id) init.tour(cities, dist.m=dist.m, alpha=c(2/3, param.grid[id,1]), p=c(2/3, param.grid[id,2])))
tours.d = sapply(seq(ncol(tours)), function(id) sum(dist.tour(cities, tours[,id])))
d.min = min(tours.d)
tour.id = seq_along(tours.d) [tours.d == d.min]
list(min=d.min, id=tour.id, alpha=param.grid[tour.id,1], p=param.grid[tour.id,2])
plot.tour(tours[ , tour.id[1]], cities, etsp, id)
# 110.01
alpha = c(2/3, -5/3) ; p = c(2/3, -1);
tour = init.tour(cities, alpha=alpha, p=p)
plot.tour(tour, cities, etsp, id)

### SA
# Note: takes long!!!
# and may NOT improve the tour!
Temp = 120; # 180; # lets start lower
tour = optim.tour(cities, tour=tour, Temp, 1E-4)
tour; tour = tour$S
plot.tour(tour, cities, etsp, id)
# starting with lower Temp and the 112.57 tour:
# => ! [but takes quit a while!]
# only 111.887; [but should be possible to get lower]


##################
##################

##############################
### 2D: 2 bands ~ N(0, sd) ###

epochs.n = 10
cities.p = r2d.gen(6, epochs.n, sep.scale=3)
plot(cities.p$x, cities.p$y)

cities = matrix(c(as.vector(cities.p$x), as.vector(cities.p$y)), ncol=2)
id = find.base(cities, y=c(0, 2))
id

### TSP data
etsp <- ETSP(cities)
etsp

### calculate a tour
# tour <- solve_TSP(etsp, method = "nn")
tour <- solve_TSP(etsp, method = "nn", control=list(start=id))
tour

tour_length(tour) # 104.91 ### Note: random data!
plot.tour(tour, cities, etsp, id)

### Tour ###
### scan parameter-space
dist.m = dist.med(cities)
param.grid = expand.grid((-5:5)/3, (-5:5)/3)
#
tour = scan.tour(cities, id, dist.m=dist.m, param.grid=param.grid, ignore.start=T)
tour
# examples: 84.8, or 90.27 # but data is random!

### Level 2:
alpha = -2/3; p = -5/3;
alpha = 4/3; p = 2/3;
scan2.tour(cities, id, alpha=alpha, p=p, dist.m=dist.m, param.grid=param.grid, ignore.start=T)
# Level 2 does NOT usually help much!

