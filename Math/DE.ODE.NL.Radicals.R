########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Exponentials w. Radicals
##
## draft v.0.1a




#########################

### y = exp(sqrt(log(x)))

x = sqrt(5)
params = list(x=x);
e = expression(exp(sqrt(log(x))))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
2*x*dy - 1/sqrt(log(x)) * y # = 0

# D2 =>
4*x^2*d2y + 4*x*dy - 1/log(x) * y + 1 / sqrt(log(x))^3 * y # = 0
x^2*y^2*d2y + x*y^2*dy - x^2*y*dy^2 + 2*x^3*dy^3 # = 0

### ODE:
x*y^2*d2y + 2*x^2*dy^3 - x*y*dy^2 + y^2*dy # = 0


#########################

### y = exp(k * sqrt(log(x))^3)

x = sqrt(5); k = sqrt(2);
params = list(x=x, k=k);
e = expression(exp(k * sqrt(log(x))^3))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
2*x*dy - 3*k*sqrt(log(x)) * y # = 0

# D2 =>
4*x^2*d2y + 4*x*dy - 9*k^2*log(x) * y +
	- 3*k / sqrt(log(x)) * y # = 0

### ODE:
8*x^3*y*dy*d2y - 8*x^3*dy^3 + 8*x^2*y*dy^2 - 9*k^2*y^3 # = 0


#############################

### y = exp(x * sqrt(log(x)))

x = sqrt(5); k = sqrt(2);
params = list(x=x, k=k);
e = expression(exp(k*x * sqrt(log(x))))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
2*dy - k*(2*sqrt(log(x)) + 1/sqrt(log(x))) * y # = 0
#
sqrt(4*dy^2 - 8*k^2*y^2) - k*(2*sqrt(log(x)) - 1/sqrt(log(x))) * y # = 0


# D2 =>
4*x*d2y - k^2*x * (2*sqrt(log(x)) + 1/sqrt(log(x)))^2 * y +
	- k * (2/sqrt(log(x)) - 1/sqrt(log(x))^3) * y # = 0
4*x*y*d2y - 4*x*dy^2 +
	- k * (2/sqrt(log(x)) - 1/sqrt(log(x))^3) * y^2 # = 0
4*x*y*d2y - 4*x*dy^2 - 2*y*dy / log(x) +
	+ 2*k / sqrt(log(x))^3 * y^2 # = 0
# Alternative:
2*x*y*d2y - 2*x*dy^2 +
	- sqrt(dy^2 - 2*k^2*y^2) * y / log(x) # = 0
2*x*log(x) * (y*d2y - dy^2) - sqrt(dy^2 - 2*k^2*y^2) * y # = 0

# TODO: find nice form;


############################

### y = x * exp(sqrt(log(x)))

# TODO

