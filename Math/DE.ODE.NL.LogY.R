#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## Log( F(Y) )-Types
##
## v.0.2c


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

### Log(y) * Log(1-y) = k*x^n

### Check
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


##################

##################
### Mixed Sums ###
##################

### x*Log(y) + y*Log(x) = k*x^n

### Check:
n = sqrt(3); # n = 3-sqrt(3); # n = 2; # n = 1;
y = exp(-1/3); k = sqrt(2);
x = 1/sqrt(5);
k = (x*log(y) + y*log(x)) / x^n;
params = list(x=x, k=k, n=n, y=y);
e = expression(x*log(y) + y*log(x))[[1]];
#
dye = D(e, "y"); dxe = D(e, "x");
dy  = (k*n*x^(n-1) - eval(dxe, params)) / eval(dye, params);
d2y = (k*n*(n-1)*x^(n-2) - dy/y + y/x^2 - dy/x    - dy/y + x/y^2*dy^2 - dy/x) / (x/y + log(x));

# D =>
x*y*log(y) + x*log(x) * y*dy + x^2 * dy + y^2 - k*n*x^n * y # = 0

# Quasi-Linear System:
log(y) # ==
(k*x^(n+1)*dy + (x^2 * dy + y^2 - k*n*x^n * y)) / (x*(x*dy - y));
#
log(x) # ==
(k*x^n*y + (x^2 * dy + y^2 - k*n*x^n * y)) / (y*(y - x*dy));

# D2 =>
y*log(y) + x*dy*log(y) +
	+ log(x) * y*dy + x*log(x) * dy^2 + x*log(x) * y*d2y +
	+ x^2 * d2y + 3*x*dy + 3*y*dy - k*n*x^n * dy - k*n^2*x^(n-1) * y # = 0
(x*dy + y) * (k*x^(n+1)*dy + x^2 * dy + y^2 - k*n*x^n * y) * y  / (x*dy - y) +
	- x*(x*y*d2y + x*dy^2 + y*dy) * (k*x^n*y + x^2 * dy + y^2 - k*n*x^n * y) / (x*dy - y) +
	+ x^3 * y*d2y + 3*x^2*y*dy + 3*x*y^2*dy - k*n*x^(n+1) * y*dy - k*n^2*x^(n) * y^2 # = 0
(- x^4*y*dy*d2y - x^2*y^3*d2y - k*x^(n+2)*y^2*d2y + n*k*x^(n+2)*y^2*d2y +
	- x^4*dy^3 - x^2*y^2*dy^2 + n*k*x^(n+2)*y*dy^2 + x^2*y^2*dy + y^4 - n*k*x^n*y^3) / (x*dy - y) +
	+ x^3 * y*d2y + 3*x^2*y*dy + 3*x*y^2*dy - k*n*x^(n+1) * y*dy - k*n^2*x^(n) * y^2 # = 0

### ODE:
x^2*(y - k*(n-1)*x^n + x) * y^2*d2y +
	+ x^4*dy^3 - 2*x^2*y^2*dy^2 - 3*x^3*y*dy^2 +
	+ x*(3*y + k*n*(n-1)*x^n + 2*x) * y^2*dy - y^4 - k*n*(n-1)*x^n*y^3 # = 0


### Special Cases:

### Case: n = 1;
x^2*(y + x) * y^2*d2y + x^4*dy^3 - x^2*(2*y + 3*x) * y*dy^2 +
	+ x*(3*y + 2*x) * y^2*dy - y^4 # = 0

### Case: n = 2;
x^2*(y - k*x^2 + x) * y^2*d2y + x^4*dy^3 - (2*x^2*y + 3*x^3) * y*dy^2 +
	+ x*(3*y + 2*k*x^2 + 2*x) * y^2*dy - y^4 - 2*k*x^2*y^3 # = 0

