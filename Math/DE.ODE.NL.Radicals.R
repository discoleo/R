########################
##
## Leonard Mada
## [the one and only]
##
## Differential Equations
## NL ODEs - Exponentials w. Radicals
##
## draft v.0.1b




#############################

### y = exp(sqrt(log(x)) / k)

x = sqrt(5); k = sqrt(2);
params = list(x=x);
e = expression(exp(sqrt(log(x)) / k))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
2*x*dy - 1/k / sqrt(log(x)) * y # = 0

# D2 =>
4*x^2*d2y + 4*x*dy - 1/k^2 / log(x) * y + 1/k / sqrt(log(x))^3 * y # = 0
x^2*y^2*d2y + x*y^2*dy - x^2*y*dy^2 + 2*k^2*x^3*dy^3 # = 0

### ODE:
x*y^2*d2y + 2*k^2*x^2*dy^3 - x*y*dy^2 + y^2*dy # = 0


### Note:
# NO Simple solution: y = x^p;
# =>
# p*(p-1) + 2*k^2*p^3 - p^2 + p = 0
# 2*k^2*p^3 # = 0


###############################

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


### Scaled Variant:
### y = exp(sqrt(8)/3 * k * sqrt(log(x))^3)

x = sqrt(5); k = sqrt(2);
params = list(x=x, k=k);
e = expression(exp(sqrt(8)/3 * k * sqrt(log(x))^3))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);


### ODE:
x^3*y*dy*d2y - x^3*dy^3 + x^2*y*dy^2 - k^2*y^3 # = 0


########
### Gen:

### y = exp(sqrt(log(x + b)) / k)

x = sqrt(5); k = sqrt(2); b = sqrt(2);
params = list(x=x);
e = expression(exp(sqrt(log(x + b)) / k))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
2*(x + b)*dy - 1/k / sqrt(log(x + b)) * y # = 0

# D2 =>
4*(x+b)^2*d2y + 4*(x+b)*dy - 1/k^2 / log(x+b) * y + 1/k / sqrt(log(x+b))^3 * y # = 0


### ODE:
(x+b)*y^2*d2y + 2*k^2*(x+b)^2*dy^3 - (x+b)*y*dy^2 + y^2*dy # = 0


################
### Gen Order 2:

### y = exp(sqrt(log(x^2 + b)) / k)

x = sqrt(5); k = sqrt(2); b = sqrt(2);
params = list(x=x, k=k, b=b);
e = expression(exp(sqrt(log(x^2 + b)) / k))[[1]];
#
y   = eval(e, params);
dy  = eval(D(e, "x"), params);
d2y = eval(D(D(e, "x"), "x"), params);

# D =>
(x^2 + b)*dy - x/k / sqrt(log(x^2 + b)) * y # = 0

# D2 =>
(x^2+b)^2*y*d2y + 2*x*(x^2+b)*y*dy - (x^2+b)^2 * dy^2 +
	- (x^2+b)/k / sqrt(log(x^2 + b)) * y^2 + x^2/k / sqrt(log(x^2+b))^3 * y^2 # = 0
x*(x^2+b)^2*y^2*d2y + k^2*(x^2 + b)^3 * dy^3 +
	- x*(x^2+b)^2 * y*dy^2 + (x^2-b)*(x^2+b) * y^2*dy # = 0

### ODE:
x*(x^2+b) * y^2*d2y + k^2*(x^2 + b)^2 * dy^3 +
	- x*(x^2+b) * y*dy^2 + (x^2-b) * y^2*dy # = 0


################
################

################
### Variants ###

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

