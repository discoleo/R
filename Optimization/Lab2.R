
### Lab 2: Optimizations
### Leonard Mada

### initial draft
### TODO: solve Problem 2 as well

###############

# install.packages("lpSolve")

library(lpSolve)

#### Problems ####

### Problem 1
### Telephone Calls

# minimize objective
objective.coeff = c(1, 1.6)
# constraints
constr.mat = matrix(c(0.3, 0.3, 0.1, 0.1,  0.3, 0.2, 0.3, 0.15), nrow=2, byrow=TRUE)
constr.mat = t(constr.mat)
constr.val = c(150, 110, 120, 100)
constr.dir <-c(">=", ">=", ">=", ">=") # direction

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

# TODO

