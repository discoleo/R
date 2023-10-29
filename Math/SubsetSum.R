

### Subset-Sum Problem

# Note:
# - the product of all permutations is a Symmetric polynomial in (x[i])^2;


####################

### Helper Functions

source("Polynomials.Helper.R")
source("Polynomials.Helper.EP.R")


####################

### 4 Variables:

p = as.pm("(x1+x2+x3+x4)*(x1+x2+x3-x4)*(x1+x2-x3+x4)*(x1+x2-x3-x4) *
	(x1-x2+x3+x4)*(x1-x2+x3-x4)*(x1-x2-x3+x4)*(x1-x2-x3-x4)")

p = sort.pm.vars(p, xn = paste0("x", 1:4))
p = sort.pm(p, xn = paste0("x", 1:4), sort.coeff = c(1,2,4,9))

print.pm(p, do.sort=F)

# Test:

# x[i] = sqrt(x) = (...)^(1/4);
x = sqrt(c(2,3,5,7))
S = sum(x)
E2 = (S^2 - sum(x^2)) / 2;
E3 = (S^3 - sum(x^3) - 3*S*E2) / 3;
E4 = prod(x); E4x = sqrt(E4);

###
S^4 - 8*S^2*E2 + 16*E2^2 - 64*E4;
(S^2 - 4*E2)^2 - 64*E4;
(S^2 - 4*E2 - 8*E4x) * (S^2 - 4*E2 + 8*E4x);
(S^4 - 4*S^2*E2 + 2*E2^2 + 4*S*E3 - 4*E4) + 6*(E2^2 - 2*S*E3 + 2*E4) +
	- 4*(E2*S^2 - E3*S - 2*E2^2 + 4*E4) + 4*(E3*S - 4*E4) - 40*E4;
eval.pm(p, sqrt(x))

