########################
###
### Leonard Mada
### [the one and only]
###
### TSP Solver
###
### draft v.0.1a


###############
### History ###

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
		dist.f(cities[tour[id], ], cities[tour[id + 1], ])
	}
	id = c(seq_along(tour), 1)
	return(sapply(seq(length(tour)-1), dist.byid))
}

dist.all = function(cities) {
	dist.byid = function(id) {
		return(sqrt( ((cities[id,1] - cities[,1])^2 + (cities[id,2] - cities[,2])^2) ))
	}
	id = 1:nrow(cities)
	return(sapply(id, dist.byid))
}

########
### Tour
init.tour = function(cities, start=1, alpha=1/3, p=1) {
	len = nrow(cities)
	remain.id = seq(len)
	d = dist.all(cities)
	d.med = sapply(1:nrow(d), function(id) median(d[id, ]))
	d = d.med
	#
	dist.reg = function(id, current.city, remain.id, d, current.max) {
		sqrt(sum((cities[remain.id[id],] - current.city)^2)) +
		alpha * (current.max^p - d[id]^p)
	}
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
    
    return (list("cost"=S_best_cost, "S"=S))
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
# - once a point is added to the toor, remove its median value
#   from the list of remaining points;
# - it is possible to use regularisation terms both for the global & the local medians;
#  -- global = all remaining locations, not yet included in the tour;
#  -- local = all locations in a local neighbourhood, not yet included in the tour;


################
################

################
### Examples ###

### Examples: Tour
phi = -0.01 # 0.2 # 0;
p = radial.gen(d=1 + (1:4)/17, r= 5 + (1:4)/3, phi=phi, addCenter=F)
plot(p$x, p$y)

cities = matrix(c(as.vector(p$x), as.vector(p$y)), ncol=2)

### with Base-city:
id = find.base(cities, y=c(0.5, 2), middle=T)
id

etsp = ETSP(cities)
etsp

### alpha = 1/5; p = 1;
### alpha = 15; p = 1/5

tour = init.tour(cities, alpha=1/5)
plot(etsp, tour, tour_col = "red", xlab="X-Coord", ylab="Y-Coord")
sum(dist.tour(cities, tour))

# Note: takes long!!!
# and may increase the cost: NO elitism in current implementation;
tour = optim.tour(cities, tour=tour, 1000, 1E-4)
sum(dist.tour(cities, tour))

plot(etsp, tour$S, tour_col = "red", xlab="X-Coord", ylab="Y-Coord")



#############################
#############################

### Examples: Analysis

### Example 1:
### variable number of points / circle
phi = -0.01 # 0.2 # 0;
p = radial.gen(d=1 + (1:7)/17, n.c=7, phi=phi)
plot(p$x, p$y)

cities = matrix(c(as.vector(p$x), as.vector(p$y)), ncol=2)

### with Base-city:
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

################
### Analysis ###

### all Locations
d = dist.all(cities)
# d.v = d.m[d.m != 0]
# summary(d.v)

d.min = sapply(1:nrow(d), function(id) min(d[id, d[id,] != 0]))
d.med = sapply(1:nrow(d), function(id) median(d[id, ]))
summary(d.min)
summary(d.med)
plot3d(p$x, p$y, d.min, type="s", col = round(10*d.min) + 1)
# nice 3D decomposition!
plot3d(p$x, p$y, d.med, type="s", col = round(10*d.min) + 1)
# less distortion
plot3d(p$x, p$y, (d.med)^(1/3), type="s", col = round(d.min + d.med) + 1)

### Tour
d.tr = dist.tour(cities, tour)
summary(d.tr)

len = tour_length(tour)
id.col = (1:len) / len
plot3d(p$x[tour], p$y[tour], d.med[tour], type="s", col = rgb(id.col, 1-id.col, 0))


#######################
#######################

### Example 2:
### variable number of points / circle
phi = -0.01 # 0.2 # 0;
p = radial.gen(d=1 + (1:4)/17, r= 5 + (1:4)/3, phi=phi, addCenter=F)
plot(p$x, p$y)

cities = matrix(c(as.vector(p$x), as.vector(p$y)), ncol=2)

### with Base-city:
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

################
### Analysis ###

### all Locations
d = dist.all(cities)
# d.v = d.m[d.m != 0]
# summary(d.v)

d.min = sapply(1:nrow(d), function(id) min(d[id, d[id,] != 0]))
d.med = sapply(1:nrow(d), function(id) median(d[id, ]))
summary(d.min)
summary(d.med)
plot3d(p$x, p$y, d.min, type="s", col = round(10*d.min) + 1)
# nice 3D decomposition!
plot3d(p$x, p$y, d.med, type="s", col = round(10*d.min) + 1)
# less distortion
plot3d(p$x, p$y, (d.med)^(1/3), type="s", col = round(d.min + d.med) + 1)

