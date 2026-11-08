#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## Log( F(Y) )-Types


####################

### Helper Functions

source("Polynomials.Helper.R")
source("DE.ODE.Helper.R")


###################

###################
### Log-Product ###
###################

### log(y + x) * log(y - x) = P(x)

### Check:
ye = expression((x+b0)^(4/3))[[1]]
pe = substitute(expression(log(y + x) * log(y - x)), list(y=ye))[[2]]
x = sqrt(3); b0 = sqrt(2);
params = list(x=x, b0=b0);
#
y = eval(ye, params); dy = eval(D(ye, "x"), params); 
p = eval(pe, params); dp = eval(D(pe, "x"), params);
d2y = eval(D(D(ye, "x"), "x"), params);
d2p = eval(D(D(pe, "x"), "x"), params);

### D =>
(dy + 1)*(y - x)*log(y - x) + (dy - 1)*(y + x)*log(y + x) - (y^2 - x^2)*dp # = 0

### D2 =>
(d2y*(y - x) + dy^2 - 1)*log(y - x) + (d2y*(y + x) + dy^2 - 1)*log(y + x) +
	+ 2*dy^2 - 2 - 2*(y*dy - x)*dp - (y^2 - x^2)*d2p # = 0

### Solve linear system:
div = 2*y^2*d2y - 2*x^2*d2y - 2*x*dy^3 + 2*y*dy^2 + 2*x*dy - 2*y;
# log(y + x)
(- (dp*y^3 - x*dp*y^2 - x^2*dp*y + x^3*dp)*d2y +
	- 2*(y-x)*dy^3 + (dp*y^2 - 2*x*dp*y - 2*y + x^2*dp + 2*x)*dy^2 +
	+ (d2p*y^3 - d2p*x*y^2 + 2*dp*y^2 - d2p*x^2*y - 4*dp*x*y + 2*y +
		+ d2p*x^3 + 2*dp*x^2 - 2*x)*dy +
	+ d2p*y^3 - d2p*x*y^2 + dp*y^2 - d2p*x^2*y - 2*dp*x*y + 2*y +
	+ x^3*d2p + x^2*dp - 2*x) / div;

# log(y - x)
((dp*y^3 + x*dp*y^2 - x^2*dp*y - x^3*dp)*d2y +
	+ 2*(y+x)*dy^3 - (dp*y^2 + 2*x*dp*y + 2*y + x^2*dp + 2*x)*dy^2 +
	- (d2p*y^3 + d2p*x*y^2 - 2*dp*y^2 - d2p*x^2*y +
		- 4*x*dp*y + 2*y - x^3*d2p - 2*x^2*dp + 2*x)*dy +
	+ d2p*y^3 + d2p*x*y^2 - dp*y^2 - d2p*x^2*y - 2*dp*x*y + 2*y +
	- x^3*d2p - x^2*dp + 2*x) / div;

# TODO: simplify / group;


### ODE:
# - Substitute log(y+x) & log(y-x) in eq:
#   log(y+x) * log(y-x) - p = 0;
# TODO: huge;


# Derivation:
pDiv = as.pm("2*(y^2*d2y - x^2*d2y - x*dy^3 + y*dy^2 + x*dy - y)");
# pLp, pLn = expressions for the 2 logs (without the division by div);
pR = as.pm("pLp() * pLn() - p*pDiv()^2");
pR = sort.dpm(pR, "y", "x");
print(pR); # 101 monomials;


### Examples:

### P(x) = x
(dy + 1)*(y - x)*log(y - x) + (dy - 1)*(y + x)*log(y + x) - (y^2 - x^2) # = 0
# =>
(d2y*(y - x) + dy^2 - 1)*log(y - x) + (d2y*(y + x) + dy^2 - 1)*log(y + x) +
	+ 2*(dy^2 - y*dy + x - 1) # = 0
# =>
div = (y^2*d2y - x^2*d2y - x*dy^3 + y*dy^2 + x*dy - y);
2*log(y + x) # =
- (y - x) * (2*(dy + 1)*(dy^2 - y*dy + x - 1) + (d2y*(y - x) + dy^2 - 1)*(y + x)) / div;
2*log(y - x) # =
+ (y + x) * (2*(dy - 1)*(dy^2 - y*dy + x - 1) + (d2y*(y + x) + dy^2 - 1)*(y - x)) / div;
# =>
2*log(y + x) # =
- (y - x)*(y^2*d2y - x^2*d2y + 2*dy^3 - y*dy^2 + x*dy^2 + 2*dy^2 +
	- 2*y*dy + 2*x*dy - 2*dy - y + x - 2) / div;
2*log(y - x) # =
+ (y + x)*(y^2*d2y - x^2*d2y + 2*dy^3 - y*dy^2 - x*dy^2 - 2*dy^2 +
	+ 2*y*dy + 2*x*dy - 2*dy - y - x + 2) / div;


### ODE:
(y^2 - x^2)*(...)*(...) +
	- 4*x*(y^2*d2y - x^2*d2y - x*dy^3 + y*dy^2 + x*dy - y)^2 # = 0

### TODO: check!

######################

### (y^2 - x^2) * log(y + x) * log(y - x) = P(x)

### Check:
ye = expression((x+b0)^(4/3))[[1]]
pe = substitute(expression((y^2 - x^2) * log(y + x) * log(y - x)), list(y=ye))[[2]]
x = sqrt(3); b0 = sqrt(2);
params = list(x=x, b0=b0);
#
y = eval(ye, params); dy = eval(D(ye, "x"), params); 
p = eval(pe, params); dp = eval(D(pe, "x"), params);
d2y = eval(D(D(ye, "x"), "x"), params);
d2p = eval(D(D(pe, "x"), "x"), params);

### D =>
(dy + 1)*(y - x)*log(y - x) + (dy - 1)*(y + x)*log(y + x) +
	+ 2*(y*dy - x) * log(y + x) * log(y - x) - dp # = 0
# =>
(dy + 1)*(y - x)*log(y - x) + (dy - 1)*(y + x)*log(y + x) +
	- dp + 2*p*(y*dy - x) / (y^2-x^2) # = 0


### D2 =>
((y - x)*log(y - x) + (y + x)*log(y + x))*(y^2-x^2)*d2y +
	+ 2*p*(y*d2y + dy^2 - 1) +
	+ (y^2-x^2)*(dy^2 - 1)*(log(y - x) + log(y + x)) +
	+ 2*(y^2-x^2)*(dy^2 - 1) +
	+ 2*(y*dy - x) * (dy+1)  * (y-x) * log(y - x) +
	+ 2*(y*dy - x) * (dy-1)  * (y+x) * log(y + x) +
	- d2p*(y^2-x^2) # = 0


# TODO


#############################
#############################

### log(y) * log(1-y) = k*x^n


###
n = sqrt(3); # n = 3-sqrt(3); # n = 2;
y = exp(-1); k = sqrt(2);
x = (log(y)*log(1-y) / k)^(1/n);
params = list(x=x, k=k, n=n, y=y);
e = expression(log(y) * log(1-y))[[1]];
#
dye = D(e, "y");
dy  = k*n*x^(n-1) / eval(dye, params);
d2y = - (dy^2 * eval(D(dye, "y"), params) - k*n*(n-1)*x^(n-2)) / eval(dye, params);

# D =>
((1-y)*log(1-y) - y*log(y)) * dy - k*n*x^(n-1) * y*(1-y) # = 0
# =>
(y*log(y)^2 - k*x^n * (1-y)) * dy + k*n*x^(n-1) * y*(1-y) * log(y) # = 0

# D2 =>
(y*log(y)^2 - k*x^n * (1-y)) * d2y +
	+ (log(y)^2 * dy + 2*log(y)*dy + k*x^n * dy - k*n*x^(n-1) * (1-y)) * dy +
	+ k*n*x^(n-1) * (1-y) * dy +
	+ k*n*x^(n-1) * (1-2*y) * log(y) * dy +
	+ k*n*(n-1)*x^(n-2) * y*(1-y) * log(y) # = 0

### ODE:
n^2*k^2*x^2*(y-1)^3*y^3*d2y^2 +
 + 4*n*k*x^(3-n)*y^4*dy^3*d2y - 8*n*k*x^(3-n)*y^3*dy^3*d2y +
 + 4*n*k*x^(3-n)*y^2*dy^3*d2y - 2*n^2*k^2*x^2*y^5*dy^2*d2y + 5*n^2*k^2*x^2*y^4*dy^2*d2y +
 - 4*n^2*k^2*x^2*y^3*dy^2*d2y + n^2*k^2*x^2*y^2*dy^2*d2y + 2*n^2*k^2*x*y^6*dy*d2y - 2*n^3*k^2*x*y^6*dy*d2y +
 - 6*n^2*k^2*x*y^5*dy*d2y + 6*n^3*k^2*x*y^5*dy*d2y + 6*n^2*k^2*x*y^4*dy*d2y - 6*n^3*k^2*x*y^4*dy*d2y +
 - 2*n^2*k^2*x*y^3*dy*d2y + 2*n^3*k^2*x*y^3*dy*d2y + 4*x^(4-2*n)*y^2*dy^6 - 4*x^(4-2*n)*y*dy^6 +
 + k*x^(4-n)*dy^6 - 4*n*k*x^(3-n)*y^3*dy^5 + 6*n*k*x^(3-n)*y^2*dy^5 +
 - 2*n*k*x^(3-n)*y*dy^5 + n^2*k^2*x^2*y^4*dy^4 + 4*n*k*x^(2-n)*y^4*dy^4 +
 - 4*n^2*k*x^(2-n)*y^4*dy^4 - 2*n^2*k^2*x^2*y^3*dy^4 - 8*n*k*x^(2-n)*y^3*dy^4 +
 + 8*n^2*k*x^(2-n)*y^3*dy^4 + n^2*k^2*x^2*y^2*dy^4 + 4*n*k*x^(2-n)*y^2*dy^4 +
 - 4*n^2*k*x^(2-n)*y^2*dy^4 - 2*n^2*k^2*x*y^5*dy^3 + 2*n^3*k^2*x*y^5*dy^3 + 5*n^2*k^2*x*y^4*dy^3 +
 - 5*n^3*k^2*x*y^4*dy^3 - 4*n^2*k^2*x*y^3*dy^3 + 4*n^3*k^2*x*y^3*dy^3 + n^2*k^2*x*y^2*dy^3 +
 - n^3*k^2*x*y^2*dy^3 + n^2*k^2*y^6*dy^2 - 2*n^3*k^2*y^6*dy^2 + n^4*k^2*y^6*dy^2 - 3*n^2*k^2*y^5*dy^2 +
 + 6*n^3*k^2*y^5*dy^2 - 3*n^4*k^2*y^5*dy^2 + 3*n^2*k^2*y^4*dy^2 - 6*n^3*k^2*y^4*dy^2 + 3*n^4*k^2*y^4*dy^2 +
 - n^2*k^2*y^3*dy^2 + 2*n^3*k^2*y^3*dy^2 - n^4*k^2*y^3*dy^2 # = 0

# TODO: format/simplify;

