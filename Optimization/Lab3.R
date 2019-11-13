
### Lab 3: Optimizations
### Leonard Mada

### v1.0

###############

# install.packages("lpSolve")

library(lpSolve)

#### Problems ####

### Problem 1
### Cat Food
### from the Pulp tutorial


### P1.1: Simplified version

# minimize objective
objective.coeff = c(0.013, 0.008)
objective.coeff = list("Chicken"=0.013, "Beef"=0.008)
# constraints
constr.mat = matrix(
	c(
	0.100, 0.080, 0.001, 0.002,
	0.200, 0.100, 0.005, 0.005), nrow=2, byrow=TRUE)

constr.mat = rbind( c(1, 1), t(constr.mat) )
constr.val = c(100, 8, 6, 2, 0.4)
constr.dir = c("=", ">=", ">=", "<=", "<=") # direction

# NON-integer solution
optimum = lp(direction="min", objective.coeff, constr.mat, constr.dir, constr.val)
optimum
optimum$solution


### P1.2: "Complete" version

# minimize objective
objective.coeff = c(0.013, 0.008, 0.002, 0.005)
objective.coeff = list("Chicken"=0.013, "Beef"=0.008, "Rice"=0.002, "Wheat_bran"=0.005)
# constraints
constr.mat = matrix(
	c(
	0.100, 0.080, 0.001, 0.002,
	0.200, 0.100, 0.005, 0.005,
	0.000, 0.010, 0.100, 0.002,
	0.040, 0.010, 0.150, 0.008), nrow=4, byrow=TRUE)

constr.mat = rbind( rep(1, 4), t(constr.mat) )
constr.val = c(100, 8, 6, 2, 0.4)
constr.dir = c("==", ">=", ">=", "<=", "<=") # direction

# NON-integer solution
optimum = lp(direction="min", objective.coeff, constr.mat, constr.dir, constr.val)
optimum
optimum$solution


#################

### Problem 2
### Gas Stations

gen.matrix = function(len, min.quantity=0, huge.quantity=0, val.default=1) {
	len_A = len[1]
	len_B = len[2]
	
	# Min constraints
	len_tot = len_A * len_B * len_B
	if(min.quantity > 0) {
		# will add dummy binary variables
		len_tot = 2 * len_tot
	}
	constr.mat = matrix(
		rep(0, len_tot),
		nrow=len_B, byrow=FALSE)
	for(i in 0:(len_B - 1)) {
		constr.mat[i + 1, 1:len_A + i*len_A] = val.default
	}
	
	# Max constraints
	tmp = matrix(rep(diag(len_A), len_B), nrow=len_A)
	if(min.quantity > 0) {
		tmp = cbind(tmp, matrix(0, nrow=len_A, ncol=len_A * len_B))
	}
	constr.mat = rbind(constr.mat, tmp)
	
	# Min quantity constraints
	if(min.quantity > 0) {
		tmp  = diag(len_A * len_B)
		tmp2 = cbind(tmp, tmp * (-min.quantity) )
		tmp2 = rbind(tmp2, cbind(-tmp, tmp * (huge.quantity) ))
		constr.mat = rbind(constr.mat, tmp2)
	}
	
	return(constr.mat)
}

# generate some distances
rdist = function(len, option=1) {
	if(option == 1) {
		val = runif(len, 10, 20)
	} else if(option == 2) {
		val = runif(len, 5, 20)
		val = val * val
	} else if(option == 3) {
		val = rnorm(len, 10, 8)
		val = val * val
	} else {
		val = rpois(len, 15)
	}
	return(val)
}

### Example

# Max quantity per supplier
supplier.q = c(200, 500, 300)
len = length(supplier.q)
# number of end users
len_gas = 9
# minimize objective
# Distance to gas stations
objective.coeff = rdist(len * len_gas, option=2)
boxplot(objective.coeff)

constr.mat = gen.matrix(c(len, len_gas))

min_val = 85 # 100
constr.val = c(rep(min_val, len_gas), supplier.q)
constr.dir = c(rep(">=", len_gas), rep("<=", len)) # direction

# NON-integer solution
optimum = lp(direction="min", objective.coeff, constr.mat, constr.dir, constr.val)
optimum
optimum$solution
t(matrix(optimum$solution, nrow=len))


#########
### Using Minimal-Quantity Constraint
min_val = 85 # 100
min_transported = 10
constr.mat = gen.matrix(c(len, len_gas), min.quantity = min_transported, huge.quantity=max(supplier.q))
# Constraints
len_tot = len * len_gas
constr.bin.val = c(constr.val, rep(0, 2*len_tot))
constr.bin.dir = c(constr.dir, rep(">=", 2*len_tot))
# Binary variables
bin.var = len_tot + 1:len_tot
objective.bin.coeff = c(objective.coeff, rep(0, len_tot))

optimum = lp(direction="min", objective.bin.coeff, constr.mat, constr.bin.dir, constr.bin.val, binary.vec = bin.var)
optimum
optimum$solution
t(matrix(optimum$solution, nrow=len))
