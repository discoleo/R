

### Subset-Sum Problem

# Note:
# - the product of all permutations is a Symmetric polynomial in (x[i])^2;


####################

### Helper Functions

source("Polynomials.Helper.R")


####################

### 4 Variables:

p = as.pm("(x1+x2+x3+x4)*(x1+x2+x3-x4)*(x1+x2-x3+x4)*(x1+x2-x3-x4) *
	(x1-x2+x3+x4)*(x1-x2+x3-x4)*(x1-x2-x3+x4)*(x1-x2-x3-x4)")

p = sort.pm.vars(p, xn = paste0("x", 1:4))
p = sort.pm(p, xn = paste0("x", 1:4), sort.coeff = c(1,2,4,9))

print.pm(p, do.sort=F)

# TODO:
# - decompose into Elementary polynomials;
