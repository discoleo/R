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
print.pm(pR, lead = NA); # 101 monomials;


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

