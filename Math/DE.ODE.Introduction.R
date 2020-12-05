########################
###
### Leonard Mada
### [the one and only]
###
### Differential Equations
### ODEs - Introduction
###
### draft v.0.1a


### Linear & Non-Linear ODEs

# - a systematic approach to ODEs;

### TODO:
# - aggregate & systematize ALL existing results;


###############

###############
### History ###
###############

### draft v.0.1a:
# - started to aggregate DE.ODE.Trigonometric.R;


###################

###################
### Linear ODEs ###
###################

############
### Case ###

### Simple Trigonometric
### y = a1*sin(P(x)) + a2*cos(P(x))
dP*d2y - d2P*dy + dP^3 * y # = 0;

### Source
# see DE.ODE.Trigonometric.R;

### Examples
### P(x) = x^2
x*d2y - dy + 2*x^2 * y # = 0
### P(x) = x^n => * x^(2-n)
x*d2y - (n-1)*dy + n*x^n * y # = 0
### P(x) = x^n + b*x
(n*x^(n-1)+b)*d2y - n*(n-1)*x^(n-2)*dy + (n*x^(n-1)+b)^3 * y # = 0
# see DE.ODE.Trigonometric.R;


#################
### Special Case:
# based on substitution: P(x) => log(P(x));

### y = a1*sin(log(P(x))) + a2*cos(log(P(x)))
### Eq:
dP*P^2 * d2y + P*(dP^2 - P*d2P)*dy + dP^3*y # = 0

### Examples:
### P(x) = x + b
(x+b)^2 * d2y + (x+b)*dy + y # = 0
### P(x) = x^n + b0
x*(x^n+b0)^2 * d2y + (x^(2*n) - (n-2)*b0*x^n - (n-1)*b0^2)*dy + n^2*x^(2*n-1)*y # = 0
# see DE.ODE.Trigonometric.R;


#######################
#######################

#######################
### Non-Linear ODEs ###
#######################


############
### Case ###

### Simple Inverse Trigonometric
### sin(y^n) = f(x)
df*y*d2y - n^2*f*y^(2*n-1)*dy^3 + (n-1)*df*dy^2 - d2f*y*dy = 0
### Alternative Eq:
(1-f^2)*df*y^(2*n-1)*d2y + (f^2*d2f - f*df^2 - d2f)*y^(2*n-1)*dy + (n-1)/n^2*df^3 = 0

### Source
# see DE.ODE.Trigonometric.R;

### Examples:

### n = 2
df*y*d2y - 4*f*y^3*dy^3 + df*dy^2 - d2f*y*dy = 0
# alternative:
df*(1 - f^2)*y^3*d2y + (d2f*f^2 - df^2*f - d2f)*y^3*dy + 1/4 * df^3 # = 0
### n = 2; f(x) = x^2;
x*(1-x^4)*y^3*d2y - (x^4+1)*y^3*dy + x^3 # = 0
# see DE.ODE.Trigonometric.R;

