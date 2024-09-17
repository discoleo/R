########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs: Polynomial types
## with Radicals
##
## draft v.0.1a



### Examples:

# (x^6+1)*dy^2 - 4*x^6 * y = 8*x^8;


######################

### Cardano-type Roots

### Pow = 6
f0 = 0; df0 = 0; # for simplicity
x = sqrt(3); # Test

y = (sqrt(x^6 + 1) + 1)^(2/3) + (sqrt(x^6 + 1) - 1)^(2/3) + f0;
dy = 2 * x^5 / sqrt(x^6 + 1) * 
	((sqrt(x^6 + 1) + 1)^(-1/3) + (sqrt(x^6 + 1) - 1)^(-1/3)) + df0;

### ODE:
dy - 2*x^3 / sqrt(x^6 + 1) * sqrt(y + 2*x^2 - f0) - df0 # = 0

# Alternative Eq:
(x^6+1)*dy^2 - 2*(x^6+1)*df0*dy - 4*x^6 * y +
	+ (x^6+1)*df0^2 + 4*x^6*f0 - 8*x^8 # = 0


### Special Cases:

###
f0 = x^2; df0 = 2*x;
(x^6+1)*dy^2 - 4*x*(x^6+1)*dy - 4*x^6 * y + 4*x^2 # = 0

###
f0 = 0; df0 = 0;
(x^6+1)*dy^2 - 4*x^6 * y - 8*x^8 # = 0


### Test:
f0  = \(x) x^2; df0 = \(x) 2*x;
yf  = \(x) (sqrt(x^6 + 1) + 1)^(2/3) + (sqrt(x^6 + 1) - 1)^(2/3) + f0(x);
dyf = \(x) 2 * x^5 / sqrt(x^6 + 1) * 
	((sqrt(x^6 + 1) + 1)^(-1/3) + (sqrt(x^6 + 1) - 1)^(-1/3)) + df0(x);
dya = \(x, eps = 1E-4) (yf(x + eps) - yf(x)) / eps;

x = c(1/3, 1/2, 1, 3/2);
dyf(x); dya(x);


