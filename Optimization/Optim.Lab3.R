
### Lab 3: Optimizations
### Leonard Mada

### v1.0

###############

# install.packages("lpSolve")

library(lpSolve)

#### Problems ####

### Problem 1
### Cat Food
# from the Pulp tutorial
# https://pythonhosted.org/PuLP/CaseStudies/a_blending_problem.html


### P1.1: Simplified version

# minimize objective
food.names = c("Chicken", "Beef")
objective.coeff = c(0.013, 0.008)
objective.coeff = as.list(objective.coeff)
names(objective.coeff) = food.names

# constraints
constr.mat = matrix(
	c(
	0.100, 0.080, 0.001, 0.002,
	0.200, 0.100, 0.005, 0.005), nrow=2, byrow=TRUE)
colnames(constr.mat) = c("Protein", "Fat", "Fiber", "Salt")
rownames(constr.mat) = food.names

constr.mat = rbind( c(1, 1), t(constr.mat) )
constr.mat
constr.val = c(100, 8, 6, 2, 0.4)
constr.dir = c("=", ">=", ">=", "<=", "<=") # direction

# NON-integer solution
optimum = lp(direction="min", objective.coeff, constr.mat, constr.dir, constr.val)
optimum
optimum$solution

apply(constr.mat, 1, function(row) sum(row * optimum$solution))


### P1.2: Complete version

# minimize objective
food.names = c("Chicken", "Beef", "Mutton", "Rice", "Wheat_bran", "Gel")
objective.coeff = c(0.013, 0.008, 0.010, 0.002, 0.005, 0.001)
objective.coeff = as.list(objective.coeff)
names(objective.coeff) = food.names
len = length(food.names)

# constraints
constr.mat = matrix(
	c(
	0.100, 0.080, 0.001, 0.002,
	0.200, 0.100, 0.005, 0.005,
	0.150, 0.110, 0.003, 0.007,
	0.000, 0.010, 0.100, 0.002,
	0.040, 0.010, 0.150, 0.008,
	0.000, 0.000, 0.000, 0.000
	), nrow=len, byrow=TRUE)
colnames(constr.mat) = c("Protein", "Fat", "Fiber", "Salt")
rownames(constr.mat) = food.names
constr.mat

constr.mat = rbind( rep(1, len), t(constr.mat) )
constr.val = c(100, 8, 6, 2, 0.4)
constr.dir = c("==", ">=", ">=", "<=", "<=") # direction

# NON-integer solution
optimum = lp(direction="min", objective.coeff, constr.mat, constr.dir, constr.val)
optimum
optimum$solution

# Chicken    Beef   Mutton   Rice  Wheat_bran   Gel 
#       0      60        0      0           0    40

rez = apply(constr.mat, 1, function(row) sum(row * optimum$solution))
round(rez, 2)

# Total  Protein   Fat   Fiber    Salt 
# 100.0     12.0   6.0     0.3     0.3


#######################
#######################

### Problem 2
### Petrol Stations

### generic function to construct the constraints matrix
gen.matrix = function(len, min.quantity=0, huge.quantity=0, val.default=1) {
	len_A = len[1]
	len_B = len[2]
	len_tot = len_A * len_B * len_B
	if(min.quantity > 0) {
		# will add dummy binary variables
		len_tot = 2 * len_tot
	}
	
	# Max constraints
	constr.mat = matrix(rep(diag(len_A), len_B), nrow=len_A)
	if(min.quantity > 0) {
		constr.mat = cbind(constr.mat, matrix(0, nrow=len_A, ncol=len_A * len_B))
	}
	
	# Min quantity constraints
	temp.mat = matrix(
		rep(0, len_tot),
		nrow=len_B, byrow=FALSE)
	for(i in 0:(len_B - 1)) {
		temp.mat[i + 1, 1:len_A + i*len_A] = val.default
	}
	constr.mat = rbind(constr.mat, temp.mat)
	rm(temp.mat)
	
	# Min shipped quantity constraints
	if(min.quantity > 0) {
		tmp  = diag(len_A * len_B)
		tmp2 = cbind(tmp, tmp * (-min.quantity) )
		tmp2 = rbind(tmp2, cbind(-tmp, tmp * (huge.quantity) ))
		constr.mat = rbind(constr.mat, tmp2)
	}
	
	return(constr.mat)
}

# helper function to generate some random "distances"/costs
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
supplier.names = paste("S", 1:len, sep="")

# number of end users
len_gas = 9
# minimize objective
# Distance to gas stations
objective.coeff = rdist(len * len_gas, option=2)
objective.coeff = matrix(objective.coeff, ncol=len, byrow=TRUE)
colnames(objective.coeff) = supplier.names
rownames(objective.coeff) = paste("E", 1:len_gas, sep="")
objective.coeff
boxplot(objective.coeff)

constr.mat = gen.matrix(c(len, len_gas))

min_val = 85 # 100
constr.val = c(supplier.q, rep(min_val, len_gas))

constr.dir = c(rep("<=", len), rep(">=", len_gas)) # direction

# NON-integer solution
optimum = lp(direction="min", objective.coeff, constr.mat, constr.dir, constr.val)
optimum
optimum$solution
# for the unmodified (un-hacked) version of lp
t(matrix(optimum$solution, nrow=len))


#########
### Using Minimal shipped-Quantity Constraint

# 3 Suppliers
# 9 Consumers
# the binary variabiles (dummy variables)
# are used to code a variable as either being 0 or non-zero!

# Constraints
min_val = 90 # 100
# lets vary the min quantity
constr.val = c(supplier.q, rep(min_val, len_gas) + 1:len_gas )
constr.dir = c(rep("<=", len), rep(">=", len_gas)) # direction

min_transported = 40
constr.mat = gen.matrix(c(len, len_gas), min.quantity = min_transported, huge.quantity=max(supplier.q))
len_tot = len * len_gas
constr.bin.val = c(constr.val, rep(0, 2*len_tot))
constr.bin.dir = c(constr.dir, rep(">=", 2*len_tot))

# Binary variables
bin.var = len_tot + 1:len_tot
objective.bin.coeff = c(objective.coeff, rep(0, len_tot))
objective.bin.coeff = matrix(objective.bin.coeff, ncol=len, byrow=TRUE)
colnames(objective.bin.coeff) = supplier.names
rownames(objective.bin.coeff) = c(paste("E", 1:len_gas, sep=""), paste("dummy", 1:len_gas, sep=""))
objective.bin.coeff

optimum = lp(direction="min", objective.bin.coeff, constr.mat, constr.bin.dir, constr.bin.val, binary.vec = bin.var)
optimum
optimum$solution
# for the unmodified (un-hacked) version of lp
t(matrix(optimum$solution, nrow=len))
