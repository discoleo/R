########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S5: Hetero-Symmetric
### Basic Solution
###
### draft v.0.1a


### Derivations
# - Basic solution;
# - Numerical derivation;
# - moved to this file from file:
#   Poly.System.S5.Ht.Formulas.Derivation.R;


####################

### Helper Functions

### Solver Tools
source("Polynomials.Helper.Solvers.Num.R")

# Function solve.S5HtMixed.Num():
# TODO: move to new file;


####################
####################

### Basic solution:
R = c(0,1,0,0,1)

### Set 1:
x0 = c(-0.70449+0.64i, 0.8913, -0.70449-0.64i, 0.2589 + 1.08i, 0.2589 - 1.08i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
print(poly.calc(x)[4], 12)
x.all = c(x)


### Set 2:
x0 = c(-1.70449-0.64i, -0.8913, -1.70449+0.64i, -0.2589 + 1.08i, -0.2589 - 1.08i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
print(poly.calc(x)[4], 12)
x.all = c(x.all, x)


### Set 3:
x0 = c(0.05-0.94i, -2.8913, 0.05+0.04i, -0.2589 + 1.08i, -0.2589 - 1.08i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
print(poly.calc(x)[4], 12)
x.all = c(x.all, x)


### Set 4:
x0 = c(-1.274+0.729i, 1.23-0.47i, 0.58-0.67i, 0.178+0.725i, -0.71 - 0.31i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
# print(poly.calc(x)[4], 12)
x.all = c(x.all, x)


### Set 5:
x0 = c(-1.274-0.729i, 1.23+0.47i, 0.58+0.67i, 0.178-0.725i, -0.71 + 0.31i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
# print(poly.calc(x)[4], 12)
x.all = c(x.all, x)


### Set 6:
x0 = c(-0.8108+1.5014i, 0.7763-1.6039i, 0.6605-0.1778i, -0.5008-0.4342i, -0.1252+0.7145i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
# print(poly.calc(x)[4], 12)
x.all = c(x.all, x)


### Set 7:
x0 = c(-0.8108-1.5014i, 0.7763+1.6039i, 0.6605+0.1778i, -0.5008+0.4342i, -0.1252-0.7145i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
# print(poly.calc(x)[4], 12)
x.all = c(x.all, x)
x.all = matrix(x.all, nc=5, byrow=T)

round0(poly.calc(x.all)) * 27
poly.calc(apply(x.all, 1, function(x) sum(x * x[c(3,4,5,1,2)]))) * 27

