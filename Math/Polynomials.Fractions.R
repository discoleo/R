########################
###
### Leonard Mada
### [the one and only]
###
### Polynomials: Fractions
### Derivations
###
### draft v.0.1a


### Polynomial Fractions
# - various Decompositions;


######################

### Helper Functions

source("Polynomials.Helper.R")


#######################
#######################


b = c(1,1,2,3,4)
p = toPoly.pm("x^5 + b[5]*x^4 + b[4]*x^3 + b[3]*x^2 + b[2]*x + b[1]");
r = roots.pm(p)
#
dp1 = dp.pm(p, by="x")
dp2 = dp.pm(dp1, by="x")


###############
### Order 1 ###
###############

### sum( 1/(x - r) )
# = d(Q(x)) / Q(x);

### Ex 1:
x = 3
sum(1/(x - r))
eval.pm(dp1, x) / eval.pm(p, x)

####################

### Decomposition of:
### x * d(Q(x))/Q(x)

### sum( r[i] / (x - r[i]) )
# = x * d(Q(x)) / Q(x) - n;

x = 3
n = length(r)
#
sum( (r - x) / (x - r) ) + sum( x / (x - r) ) # =
- n + sum( x / (x - r) ) # =
# x * d(Q(x)) / Q(x) - n;

### Ex 1:
x = 3
n = length(r)
sum( r / (x - r) )
x * eval.pm(dp1, x) / eval.pm(p, x) - n;


###############

###############
### Order 2 ###
###############

### sum( 1/((x - r[i])*(x - r[j])) )
# = d2(Q(x)) / (2 * Q(x))

### Ex 1:
x  = 3
xd = (x - r);
xF2 = combn(xd, 2);
#
sum( apply(1/xF2, 2, prod) )
eval.pm(dp2, x) / (2 * eval.pm(p, x))


#####################

### Decomposition of:
### x * d2(Q(x))/Q(x)

### sum( (r[i] + r[j]) / ((x - r[i])*(x - r[j])) )
# = (x*d2(Q(x)) - 4*d(Q(x))) / Q(x);

# sum( (r[i] + r[j] - 2*x + 2*x) / ((x - r[i])*(x - r[j])) ) # =
# - 4*sum(1 / (x - r[i])) + 2*x*sum( 1/((x - r[i])*(x - r[j])) ) # =
# x * d2(Q(x)) / Q(x) - 4*d(Q(x)) / Q(x) # =
# (x*d2(Q(x)) - 4*d(Q(x))) / Q(x);

### Ex 1:
x  = 3
xd = (x - r);
xF2 = combn(xd, 2);
#
sum( apply(xF2, 2, function(tx) (2*x - sum(tx)) / prod(tx)) )
(x*eval.pm(dp2, x) - 4*eval.pm(dp1, x)) / eval.pm(p, x)


#####################

### Decomposition of:
### x^2 * d2(Q(x))/Q(x)

### sum( (r[i]*r[j]) / ((x - r[i])*(x - r[j])) )
# = choose(n, 2) + (1/2*x^2*d2(Q(x)) - 4*x*d(Q(x))) / Q(x);

# sum( ((x-r[i])*(x-r[j]) - x^2 + x*(r[i]+r[j])) / ((x - r[i])*(x - r[j])) ) # =
# choose(n, 2) - sum( (x^2 - x*(r[i]+r[j])) / ((x - r[i])*(x - r[j])) ) # =
# choose(n, 2) - 1/2*x^2 * d2(Q(x)) / Q(x) + x*(x*d2(Q(x)) - 4*d(Q(x))) / Q(x)
# choose(n, 2) + 1/2*x^2 * d2(Q(x)) / Q(x) - 4*x*d(Q(x)) / Q(x)

### Ex 1:
x  = 3
n  = length(r)
xd = (x - r);
xF2 = combn(xd, 2);
#
sum( apply(xF2, 2, function(tx) prod(x - tx) / prod(tx)) )
choose(n, 2) + (1/2*x^2*eval.pm(dp2, x) - 4*x*eval.pm(dp1, x)) / eval.pm(p, x)


