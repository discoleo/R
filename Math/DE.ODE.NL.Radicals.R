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


