########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODE From another Equation
###
### draft v.0.2b

# - usually non-linear ODEs;


####################

### Helper Functions

source("Polynomials.Helper.ODE.R")
source("DE.ODE.Helper.R")


#########################
#########################

### Basic Examples
# see file: DE.ODE.Polynomial.R


### Other Examples:
# - spread across various files (DE.ODE. ...);


####################

### 2*atan((y + 2) / (1 - x)) - log(y + 2) = Const
# Example from:
# Maths 505: A very interesting differential equation
# https://www.youtube.com/watch?v=zDw4FjmLv0Y

### ODE:
dy - 2*(y + 2)^2 / (y + x + 1)^2 # = 0


##################

### y * exp(x*y) + x * exp(y^2) = F(x)

# for testing:
x = sqrt(5); y = x; f = 2*x*exp(x^2); dy = 1; d2y = 0;
df = 2*(1 + 2*x^2)*exp(x^2); d2f = 2*(4*x + 2*x*(1 + 2*x^2))*exp(x^2);

# D =>
(dy + x*y*dy + y^2) * exp(x*y) + (1 + 2*x*y*dy) * exp(y^2) - df # = 0
# x*exp(y^2) = - y * exp(x*y) + f;
x*(dy + x*y*dy + y^2) * exp(x*y) - (2*x*y*dy + 1) * (y * exp(x*y) - f) - x*df # = 0
(x*dy - 2*x*y^2*dy + x^2*y*dy + x*y^2 - y) * exp(x*y) + 2*x*f*y*dy + f - x*df # = 0

# D2 =>
(x^2*y*d2y - 2*x*y^2*d2y + x*d2y - 2*x^2*y^2*dy^2 + x^3*y*dy^2 - 4*x*y*dy^2 + 2*x^2*dy^2 +
	- 2*x*y^3*dy + 2*x^2*y^2*dy - 2*y^2*dy + 2*x*y*dy + 2*x*y*dy + x*y^3) * exp(x*y) +
	+ 2*x*f*y*d2y + 2*x*f*dy^2 + 2*f*y*dy + 2*x*df*y*dy - x*d2f # = 0


# TODO:
# - substitute exp(x*y) and check;


######################
######################

### P(x, y) = B(x) * Trig( Log ) + F(x)

### y^2 + x*y = sin(log(x)) + x^3
genODE.TrigLog.pm(as.pm(1), 0, pT=as.pm("x"), f0=as.pm("x^3"), pMxy = as.pm("y^2 + x*y"))
# ODE:
2*x^2*y*d2y + x^3*d2y + 2*x^2*dy^2 + 2*x*y*dy + 3*x^2*dy + y^2 + 2*x*y - 10*x^3 # = 0

# TODO: check;


### y^2 + x*y = x * sin(log(x)) + x^3
genODE.TrigLog.pm(as.pm("x"), 0, pT=as.pm("x"), f0=as.pm("x^3"), pMxy = as.pm("y^2 + x*y"))
# ODE:
2*x^2*y*d2y + x^3*d2y + 2*x^2*dy^2 - 2*x*y*dy + x^2*dy + 2*y^2 + x*y - 5*x^3 # = 0

# TODO: check;


### y^2 + x*y = (x+1) * sin(log(x+1)) + 1
tmp = genODE.TrigLog.pm(as.pm("x+1"), 0, pT=as.pm("x+1"), f0=as.pm(1), pMxy = as.pm("y^2 + x*y"))
print.dpm(div.pm(tmp, "(x+1)^2")$Rez)
# ODE:
2*(x + 1)^2*y*d2y + x*(x + 1)^2*d2y + 2*(x + 1)^2*dy^2 +
	- 2*(x + 1)*y*dy + (x^2 + 3*x + 2)*dy + 2*y^2 + x*y - y - 2 # = 0

# TODO: check;


########################

### y + x*tan(y) = F0(x)

# D =>
dy + tan(y) + x*(1 + tan(y)^2)*dy - df0 # = 0
x*dy + x*tan(y) + x^2*dy + x^2*tan(y)^2*dy - x*df0 # = 0
x*dy + f0 - y + x^2*dy + (f0 - y)^2*dy - x*df0 # = 0
### ODE:
y^2*dy - 2*f0*y*dy + (f0^2 + x^2)*dy - y - x*df0 + f0 # = 0

# TODO: check;

# Case: f0 = x
y^2*dy - 2*x*y*dy + 2*x^2*dy - y # = 0


##########################

### y + x*tan(y^2) = F0(x)

# D =>
dy + tan(y^2) + 2*x*(1 + tan(y^2)^2)*y*dy - df0 # = 0
x*dy + x*tan(y^2) + 2*x^2*y*dy + 2*x^2*tan(y^2)^2*y*dy - x*df0 # = 0
x*dy + f0 - y + 2*x^2*y*dy + 2*(f0 - y)^2*y*dy - x*df0 # = 0
### ODE:
2*y^3*dy - 4*f0*y^2*dy + 2*(f0^2 + x^2)*y*dy - y - x*df0 + f0 # = 0

# TODO: check;

# Case: f0 = x
2*y^3*dy - 4*x*y^2*dy + 4*x^2*y*dy - y # = 0
# =>
2*y^2*dy - 4*x*y*dy + 4*x^2*dy - 1 # = 0


########################

### y + x*tan(log(y)) = F0(x)

# D =>
dy + tan(log(y)) + x*(1 + tan(log(y))^2)*dy / y - df0 # = 0
x*y*dy + x*y*tan(log(y)) + x^2*dy + x^2*tan((log(y)))^2*dy - x*df0*y # = 0
x*y*dy + (f0 - y)*y + x^2*dy + (f0 - y)^2*dy - x*df0*y # = 0
### ODE:
y^2*dy - (2*f0 - x)*y*dy + (f0^2 + x^2)*dy - y^2 - (x*df0 - f0)*y # = 0

# TODO: check;

# Case: f0 = x
y^2*dy - x*y*dy + 2*x^2*dy - y^2 # = 0

