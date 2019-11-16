
### Lab 2: Optimizations
### Leonard Mada

### v1.0
### [includes solution for Problem 2]

###############

# install.packages("lpSolve")

library(lpSolve)

#### Problems ####

### Problem 1
### Telephone Calls

# minimize objective
objective.coeff = c(1, 1.6)
objective.coeff = list("Morning"=1, "Evening"=1.6)
# constraints
constr.mat = matrix(c(0.3, 0.3, 0.1, 0.1,  0.3, 0.2, 0.3, 0.15), nrow=2, byrow=TRUE)
constr.mat = t(constr.mat)
constr.val = c(150, 110, 120, 100)
constr.dir = c(">=", ">=", ">=", ">=") # direction

# NON-integer solution
optimum = lp(direction="min", objective.coeff, constr.mat, constr.dir, constr.val)
optimum
optimum$solution

# Integer solution
# we actually need the integer solution
# all.int=TRUE
optimum = lp(direction="min", objective.coeff, constr.mat, constr.dir, constr.val, all.int=TRUE)
optimum
optimum$solution
optimum$objval

#################

### Problem 2
### Number of Nurses

days.of.week <- weekdays(as.Date(4,"1970-01-01",tz="GMT")+0:6)

# constraints
x.row = rep(1, 7)
x.row[6] = 0
x.row[7] = 0
# a lot of code to create "automatically" the matrix
constr.mat = matrix(x.row, nrow=7, ncol=7, byrow=TRUE)
rownames(constr.mat) = days.of.week

m.shift = function(by) {
	r = constr.mat[by , ((1:7 - by) %% 7) + 1, drop=FALSE]
	# rownames(r) = rownames(constr.mat)[by]
	return(r)
}
# preserve name
l = lapply(1:7, m.shift)
constr.mat = t(simplify2array(l, higher=FALSE))
rownames(constr.mat) = sapply(l, rownames)
constr.mat


# minimize objective
objective.coeff = rep(1, 7)
objective.coeff = lapply(1:7, function(x) objective.coeff[[x]])
names(objective.coeff) = days.of.week
# constraints
constr.mat = t(constr.mat)
constr.val = c(17, 13, 15, 19, 14, 16, 11)
constr.dir = rep(">=", 7) # direction

# Integer solution
optimum = lp(direction="min", objective.coeff, constr.mat, constr.dir, constr.val, all.int=TRUE)
optimum
optimum$solution
optimum$objval


### Experiments with feasability of system
N_EXACT = 3
constr.dir = sample(c(rep(">=", 7 - N_EXACT), rep("==", N_EXACT)))
constr.dir # direction

optimum = lp(direction="min", objective.coeff, constr.mat, constr.dir, constr.val, all.int=TRUE)
optimum
optimum$solution
optimum$objval

### All Permutations
permutations = function( x, prefix = c() ) {
	if(length(x) == 0 ) return(prefix)
	do.call(rbind, sapply(1:length(x),
		FUN = function(idx) permutations( x[-idx], c( prefix, x[idx])), simplify = FALSE))
}
constr.dir.perm = unique(permutations(c(rep(">=", 7 - N_EXACT), rep("==", N_EXACT))))
constr.dir.perm

apply(constr.dir.perm, 1, function(constr.dir) {
	lp(direction="min", objective.coeff, constr.mat, constr.dir, constr.val, all.int=TRUE)} )
