

### Trig: Mixed Fractions


####################

### Helper Functions

int = function(FUN, upper, lower, rel.tol=1E-8, subdivisions=4097, ...) {
	r = integrate(FUN, upper, lower, subdivisions=subdivisions, rel.tol=rel.tol, ...);
	cat(paste0("Err = ", r$abs.error, "; subd = ", r$subdivisions, "\n"));
	return(r$value);
}

####################

### I( 1 / (x - sin(x)) )
# Fails:
integrate(\(x) 1 / (x - sin(x)) - 6/x^3 - 3/10/x, 0, 1)
# pracma: did not finish for > 5 minutes
# pracma::integral(\(x) 1 / (x - sin(x)) - 6/x^3 - 3/10/x, 0, 1)

# Interval split:
mid = 0.0075; up = 1;
integrate(\(x) 1 / (x - sin(x)) - 6/x^3 - 3/10/x, mid, 10*mid, rel.tol=1E-7, subdivisions=65)$value +
	+ integrate(\(x) 1 / (x - sin(x)) - 6/x^3 - 3/10/x, 10*mid, up, rel.tol=1E-8)$value +
	+ 11/2800 * mid^2 + 17/(4*126000) * mid^4;

# Note:
# Base Int correction: - 3/x^2 + 3/10*log(x);
# - Case: up = 1 => - 3;

# TODO: ???


###############

### on [0, Inf]
FUN = \(x) 1 / (x - sin(x)) - 6/x^3 - 3/10/x - 7/10/(x+1);
mid = 0.01;
int(FUN, mid, 1, rel.tol=1E-7, subdiv=8025) +
int(FUN, 1, 1E+5, rel.tol=1E-7, subdiv=8025) +
int(FUN, 1E+5, 1E+6, rel.tol=1E-6, subd=1025) +
int(FUN, 1E+6, 5E+6, rel.tol=1E-6, subd=1025) +
+ 11/2800 * mid^2 + 17/(4*126000) * mid^4;

# TODO: ???


###
integrate(\(x) 1 / (x - sin(pi/2*x)) - 1/(1-pi/2)/x - 1/(x-1), 0, 1)

# TODO: ???


##################

###
# 1. Maths 505: Feynman's technique is unreasonably OP!
#    https://www.youtube.com/watch?v=xOtQ8Mh0cvg
#  => the sub-integral (NOT Feynman);
# 2. Maths 505: 2 ridiculously awesome integrals!
#    https://www.youtube.com/watch?v=SwizwPy-GmE

### I( sin(x) / (x * (1 + a*cos(x))) ) on [0, Inf]
a = sqrt(3)/5;
# upper = Inf: numerical issues;
integrate(\(x) sin(x) / (x * (1 + a*cos(x))), 0, 2000, subdivisions = 1025)
pi/(2*a) + (a-1)*pi/(2*a*sqrt(1-a^2))

# Gen: I( sin(x) / (x * (b + a*cos(x))) )
a = sqrt(3)/5; b = sqrt(5);
# upper = Inf;
integrate(\(x) sin(x) / (x * (b + a*cos(x))), 0, 2000, subdivisions = 1025)
pi/(2*a) + (a - b)*pi/(2*a*sqrt(b^2 - a^2))


### I( sin(x) / (x * (b + cos(x)^2)) )
b = 1
# upper = Inf: numerical issues;
integrate(\(x) sin(x) / (x * (b + cos(x)^2)), 0, 16000, subdivisions=10029)
pi/2 / sqrt(b*(b+1))


###
b = sqrt(3)
# upper = Inf: numerical issues;
integrate(\(x) sin(x) / (x * (b + cos(x)^2)), 0, 20000, subdivisions=10029)
pi/2 / sqrt(b*(b+1))


### I( sin(x^n) / (x * (b + cos(x^n)^2)) )
b = sqrt(3); n = sqrt(5);
# upper = Inf: extreme numerical issues;
# takes 1-2 minutes with pracma;
pracma::integral(\(x) sin(x^n) / (x * (b + cos(x^n)^2)), 0, 1000)
pi/(2*n) / sqrt(b*(b+1))


### I( sin(x) / (x * (b + cos(x)^2)^2) )
b = sqrt(3)
# upper = Inf: numerical issues;
integrate(\(x) sin(x) / (x * (b + cos(x)^2)^2), 0, 16000, subdivisions=10029)
pi/4 * (2*b+1) / (b*(b+1))^(3/2)

