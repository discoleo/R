#########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Log * Radical
##
## draft v.0.1b


### Log to Power & Radical:
# y = Log(P1(x))^n1 * P2(x)^(1/n2)



### Examples:

# 2*x*(x+b)^2 * y*d2y - x*(x+b)^2 * dy^2 - 2*(x+b)*((k-1)*x - b) * y*dy +
#   + k*(k*x - 2*b) * y^2 = 0;


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### y = Log(x)^2 * (x+b)^k

# Check:
k = 1/sqrt(5); # k = 2; # k = 1/2;
b = exp(-1);
x = sqrt(3); params = list(x=x, b=b, k=k);
e = expression(log(x)^2 * (x+b)^k)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*(x+b)*dy - k*x*y - 2*(x+b) * log(x)*(x+b)^k # = 0

# D2 =>
x^2*(x+b)*d2y - x*((k-2)*x - b)*dy - k*x*y +
	- 2*(k+1)*x * log(x)*(x+b)^k +
	- 2*(x+b) *(x+b)^k # = 0
x^2*(x+b)^2 * d2y - x*(x+b)*((2*k-1)*x - b) * dy +
	+ k*x*(k*x - b) * y - 2*(x+b)^2 *(x+b)^k # = 0
2*x^2*(x+b)^2 * d2y - 2*x*(x+b)*((2*k-1)*x - b) * dy +
	+ 2*k*x*(k*x - b) * y - (x*(x+b)*dy - k*x*y)^2 / y # = 0

### ODE:
2*x*(x+b)^2 * y*d2y - x*(x+b)^2 * dy^2 +
	- 2*(x+b)*((k-1)*x - b) * y*dy + k*(k*x - 2*b) * y^2 # = 0

### Special Cases:

### Case: k = 2;
2*x*(x+b)^2 * y*d2y - x*(x+b)^2 * dy^2 +
	- 2*(x+b)*(x-b) * y*dy + 4*(x-b) * y^2 # = 0

### Case: k = 1/2;
8*x*(x+b)^2 * y*d2y - 4*x*(x+b)^2 * dy^2 +
	+ 4*(x+b)*(x + 2*b) * y*dy + (x - 4*b) * y^2 # = 0


##########################

### y = Log(x)^n * (x+b)^k

# Check:
k = 1/sqrt(5); # k = 2; # k = 1/2;
n = exp(2/3);  # n = 2*k; # n = 1/2; k = 1/4;
b = exp(-1);
x = sqrt(3); params = list(x=x, b=b, k=k);
e = expression(log(x)^n * (x+b)^k)[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
x*(x+b)*dy - k*x*y - n*(x+b) * log(x)^(n-1) * (x+b)^k # = 0

# D2 =>
x^2*(x+b)*d2y - x*((k-2)*x - b)*dy - k*x*y +
	- n*(k+1)*x * log(x)^(n-1) * (x+b)^k +
	- n*(n-1)*(x+b) * log(x)^(n-2) *(x+b)^k # = 0
x^2*(x+b)^2 * d2y - x*(x+b)*((2*k-1)*x - b) * dy +
	+ k*x*(k*x - b) * y +
	- n*(n-1)*(x+b)^2 * log(x)^(n-2) *(x+b)^k # = 0

### ODE:
n*x*(x+b)^2 * y*d2y - (n-1)*x*(x+b)^2 * dy^2 +
	- (x+b)*((2*k-n)*x - b*n) * y*dy +
	+ k*(k*x - b*n) * y^2 # = 0

### Special Cases:

### Case: n = 2*k;
2*k*x*(x+b)^2 * y*d2y - (2*k-1)*x*(x+b)^2 * dy^2 +
	+ 2*k*b*(x+b) * y*dy + k^2*(x-2*b) * y^2 # = 0

# Example: n = 1/2; k = 1/4;
8*x*(x+b)^2 * y*d2y + 8*x*(x+b)^2 * dy^2 +
	+ 8*b*(x+b) * y*dy + (x-2*b) * y^2 # = 0

