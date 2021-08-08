########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Asymmetric S2:
### Binomial Expansions
###
### draft v.0.1a


### Asymmetric Polynomial Systems: 2 Variables
### Binomial Expansions


### Base: Class 1 Polynomials


###############
### History ###
###############


### draft v.0.1a:
# - Order 3;


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R;
# e.g. round0(), round0.p;

solve.Cardano = function(c, d, n=3) {
	m = unity(n=n, all=TRUE);
	xdet = rootn(d^2 - c^n, 2);
	p = d + xdet; p = rootn(p, n);
	q = d - xdet; q = rootn(q, n);
	if(n %% 2 == 0 && c < 0) q = -q;
	sol = p*m + q/m;
	return(sol);
}


##########################

##########################
### Polynomial Systems ###
##########################

#####################
### Base: Class 1 ###
#####################

# Binomial expansions of Class 1 polynomials;

###############
### Order 3 ###
###############

### Base:
# x = k^2 + b11*k
# y = k^2 + b21*k
# where k^3 = K;

### System:
# x^3 + y^3 - 3*b11*K*x - 3*b21*K*y = 2*K^2 + (b11^3 + b21^3)*K
# x*y*(x+y) - (b11 + 2*b21)*K*x - (2*b11 + b21)*K*y = 2*K^2 + b11*b21*(b11+b21)*K


### Derivation:

### P(x^3) + P(y^3)
x^3 + y^3 - 3*b11*K*x - 3*b21*K*y - 2*K^2 - (b11^3 + b21^3)*K # = 0

### (x+y)^3 - (x^3+y^3)
3*x*y*(x+y) - 6*(b11+b21)*K*(x+y) + 3*b11*K*x + 3*b21*K*y - 6*K^2 - 3*b11*b21*(b11+b21)*K # = 0
x*y*(x+y) - 2*(b11+b21)*K*(x+y) + b11*K*x + b21*K*y - 2*K^2 - b11*b21*(b11+b21)*K # = 0
x*y*(x+y) - (b11 + 2*b21)*K*x - (2*b11 + b21)*K*y - 2*K^2 - b11*b21*(b11+b21)*K # = 0

### Note:
### Eq 1 + 3*Eq 2 => S = (x+y) =>
#   S^3 - 6*(b11+b21)*K*S + ... = 0;
# - enables solving for S = (x+y);
### Step 2:
# x + y = S;
# x*y*S - b11*K*x - b21*K*y = ... - S^3 / 3;


### Solver:
# simple variant of solver (non-robust);
solve.DP3 = function(K, b, all=FALSE) {
	c.f = function(b) b*K;
	d.f = function(b) (K^2 + b^3*K) / 2;
	# can also use the direct formulas;
	x = solve.Cardano(c.f(b[1]), d.f(b[1]), n=3);
	y = solve.Cardano(c.f(b[2]), d.f(b[2]), n=3);
	if(all) {
		# is actually NOT correct
		sol = expand.grid(x, y);
	} else sol = cbind(x=x, y=y);
	return(sol);
}

### Examples:

### Ex 1:
b = c(1,-2)
K = 3
#
sol = solve.DP3(K, b);
x = sol[,1]; y = sol[,2];
b11 = b[1]; b21 = b[2];
# concrete Example
round0(x^3 + y^3 - 3*K*x + 6*K*y - 2*K^2 + 7*K) # = 0
round0(x*y*(x+y) + 3*K*x - 2*K^2 - 2*K) # = 0


### Ex 2:
b = c(-2,3)
K = 3
#
sol = solve.DP3(K, b, all=T);
x = sol[,1]; y = sol[,2];
b11 = b[1]; b21 = b[2];
### Test
# only 3 solutions are correct, but not necessarily the base-set;
x^3 + y^3 - 3*b11*K*x - 3*b21*K*y # = 2*K^2 + (b11^3 + b21^3)*K
err = x*y*(x+y) - (b11 + 2*b21)*K*x - (2*b11 + b21)*K*y - 2*K^2 - b11*b21*(b11+b21)*K
round0(err)

