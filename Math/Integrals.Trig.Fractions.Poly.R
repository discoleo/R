


### Helper

Catalan = 0.915965594177219015054603514;


int.FrU01 = function(n, p=0) {
	(digamma(((p+1)/n + 1)/2) - digamma((p+1)/n/2)) / (2*n);
}
int.FrDU01 = function(n, p=0) {
	# TODO: p != 0
	digamma(1/n) + Euler + log(n);
}

zeta = function(x) {
	pracma::zeta(x);
}


####################
####################

###
# Michael Penn: Check out this crazy integral trick!
# https://www.youtube.com/watch?v=nBGcO3VL8Kw

# Upper = Inf: numeric instability;
pracma::integral(\(x) sin(x) / (x + 1), 0, 1000000)
integrate(\(x) 1 / (log(x)^2 + 1), 0, 1)
# TODO: ???


####################
####################

# qncubed3:  TWO WAYS to destroy this INSANE integral feat. @maths_505
# https://www.youtube.com/watch?v=Ry5q4NsZDx0

integrate(function(x) sin(x)^2 / (x^2*(x^2+1)), 0, Inf)
pi/4*(1 + exp(-2))


### Helper:
k = sqrt(3)
integrate(\(x) cos(k*x) / (x^2 + 1), 0, Inf, subdivisions=4096*2, rel.tol=1E-6)
pi/2 * exp(-k)


### I( x * sin(k*x) / (x^2 + 1) )
k = sqrt(3)
# Upper = Inf: numerical instability;
pracma::integral(\(x) x * sin(k*x) / (x^2 + 1), 0, 200000)
pi/2 * exp(-k)



### Gen 1:
k = sqrt(3)
b = sqrt(5)
integrate(function(x) cos(k*x) / (x^2 + b^2), 0, Inf, subdivisions=4096*2, rel.tol=1E-6)
pi/2 * exp(-k*b)/b


### TODO
k = 1
integrate(function(x) sin(k*x) / (x^2 + 1), 0, Inf, subdivisions=4096*2, rel.tol=1E-5)


### Gen 1:
k = sqrt(3)
integrate(function(x) sin(k*x)^2 / (x^2 + 1), 0, Inf, subdivisions=4096, rel.tol=1E-6)
pi/4*(1 - exp(-2*k))


### Gen 2:
k = sqrt(3)
integrate(function(x) sin(k*x)^2 / (x^2*(x^2 + 1)), 0, Inf, subdivisions=1024, rel.tol=1E-8)
pi/4*(2*k - 1 + exp(-2*k))

# f = (1 - exp(2i*k*z)) / (z^2 * (z^2 + 1));
# Res at 1i = pi*(1 - exp(-2*k));


### I( sin(k*x) / (x*(x^2 + 1)) )
# D(k) =>
k = sqrt(3)
# Upper = Inf: Numeric instability;
pracma::integral(function(x) sin(k*x) / (x*(x^2 + 1)), 0, 100000)
pi/2*(1 - exp(-k))


### I( sin(k*x)^3 / (x*(x^2 + 1)) )
pracma::integral(function(x) sin(k*x)^3 / (x*(x^2 + 1)), 0, 100000)
pi/8*(2 - 3*exp(-k) + exp(-3*k))


### Helper:

### I( cos(k*x)/x - exp(-x)/x )
k = 3
# upper = Inf; numerical issues;
pracma::integral(\(x) cos(k*x)/x - exp(-x)/x, 0, 100000, no_=1024)
- log(k)

###
k = 3
# upper = Inf; numerical issues;
pracma::integral(\(x) cos(k*x)/x^2 - 1/x^2, 0, 200000)
- k*pi/2


########################

### Numerical Stability:

k = sqrt(3)
start = 0; steps = 1000
r = sapply(seq(start, steps), function(id) {
	startI = 2*pi*id/k;
	endI = 2*pi*(id+1)/k;
	integrate(function(x) cos(k*x) / (x^2 + 1), startI, endI, subdivisions=4096*2, rel.tol=1E-8)$value
})
r = sum(r)
# FAILS miserably!
r + integrate(function(x) cos(k*x) / (x^2 + 1), 2*pi*(steps+1)/k, Inf, subdivisions=4096*2, rel.tol=1E-8)$value
# residual: should be < 1E-8;
print(r, 12)
print(pi/2*exp(-k), 12)


######################
######################

### Basic

### sin(x) / x^(s+1)
# 1. Maths 505: A RIDICULOUSLY AWESOME INTEGRAL: Ramanujan vs Maths 505
#    https://www.youtube.com/watch?v=_VkRvuSxF18
# 2. Michael Penn: A nice integral.
#    https://www.youtube.com/watch?v=nkaZEI_e2SU

# see also:
# 3Blue1Brown: Researchers thought this was a bug (Borwein integrals)
# https://www.youtube.com/watch?v=851U557j6HE


# s in (-1, 1)
s = 1 - sqrt(2)
# FAILS:
integrate(\(x) sin(x) / x^(s+1), 0, Inf, subdivisions=128)
integrate(\(x) sin(x) / x^(s+1), 0, 200000, subdivisions=400000)
pracma::integral(\(x) sin(x) / x^(s+1), 0, 2^16*pi, no_intervals=1024)
- gamma(-s) * sin(pi*s/2)

# s in (-1, 1)
s = 2 - sqrt(2)
integrate(\(x) sin(x) / x^(s+1), 0, 200000, subdivisions=40000)
- gamma(-s) * sin(pi*s/2)

# s = 0
# integrate(\(x) sin(x) / x, 0, Inf, subdivisions=128)
integrate(\(x) sin(x) / x, 0, 200000, subdivisions=400000)
pi/2


### I( cos(k*x) / x^s )
# s in (0, 2) (with Gamma-term);
# but numerically until s = 1.52...;
s = 2 - sqrt(2); k = 3;
pracma::integral(\(x) cos(k*x) / x^s - exp(-x)/x^s, 0, 1, no_=1025) +
	pracma::integral(\(x) cos(k*x) / x^s - exp(-x)/x^s, 1, 200000);
gamma(1-s) * sin(pi*s/2) * k^(s-1) - gamma(1-s);

# Lim: s -> 1
s = 1 - 1E-6;
gamma(1-s) * sin(pi*s/2) * k^(s-1) - gamma(1-s);
- log(k);

# Note: Compute:
# lim (gamma(x) * sin(pi*x)) as x -> 0;
# lim (gamma(x) * sin(x)) as x -> 0;


### Fresnel-Type
# see file: Integrals.Trig.Fresnel.R;

### n = 1:
# formula works for n = 1 as well;
# 0 < p < 2;
p = 3/4
# Up = Inf; numerically unstable;
pracma::integral(\(x) sin(x) / x^p, 0, 20000)
sin(pi*(1 - p)/2) * gamma(1 - p)


### Lim: p -> 1
pracma::integral(\(x) sin(x) / x, 0, 20000)
# sin(pi*(1 - p)/2) * gamma(1 - p)
pi/2


### Pow = 2
### I( sin(k*x)^2 / x^(p+1) )

###
# p in (0, 2)
p = sqrt(2)
# Up = Inf; numerical issues;
integrate(\(x) sin(x)^2 / x^(p+1), 0, 200000, subdivisions=40000)
integrate(\(x) 1/p * sin(2*x) / x^p, 0, 200000, subdivisions=40000)
- gamma(-(p-1)) * sin(pi*(p-1)/2) * 2^(p-1) / p

###
p = sqrt(2)
k = 1/5; # numerical issues with irrational numbers;
integrate(\(x) sin(k*x)^2 / x^(p+1), 0, 200000, subdivisions=40000)
- gamma(-(p-1)) * sin(pi*(p-1)/2) * 2^(p-1) * k^p / p


###
k = sqrt(3)
integrate(\(x) sin(k*x)^2 / x^2, 0, Inf, subdivisions=1025)
pi*k/2


### Pow = 3
### I( sin(k*x)^3 / x^p )

###
# p in (0, 4)
p = sqrt(2)
integrate(\(x) sin(x)^3 / x^p, 0, Inf, subdivisions=4000)
gamma(-(p-1)) * sin(pi*(p-1)/2) * (3^(p-1) - 3) / 4;


###
p = sqrt(2); k = 7/16;
# extensive numerical issues;
integrate(\(x) sin(k*x)^3 / x^p, 0, Inf, rel.tol=1E-5, subdivisions=40000)
pracma::integral(\(x) sin(k*x)^3 / x^p, 0, 32*pi*1000) +
	pracma::integral(\(x) sin(k*x)^3 / x^p, 32*pi*1000, 32*pi*10000)
gamma(-(p-1)) * sin(pi*(p-1)/2) * (3^(p-1) - 3) * k^(p-1) / 4;


# Note: FAILS
p = sqrt(2); k = 7/16;
integrate(\(x) sin(k*x)^3 / x^p, 32*pi*1000, 32*pi*10000)


### Borwein Integrals:
# - already numerical issues;
integrate(\(x) sin(x) / x * sin(x/3) / x, 0, Inf, rel.tol=1E-6, subdivisions=4000)
pi/2 / 3


################
### on [0, pi/2]

###
integrate(\(x) sin(x) / x, 0, pi/2)
pi/2 - integrate(\(x) x * exp(-pi/2 * x) / (x^2+1), 0, Inf, rel.tol=1E-8)$value
# TODO: ???
# Note: Walpha reverses the order of the integrals;


######################
######################

#################
### "Inverse" ###
#################

library(pracma)

Catalan = 0.915965594177219015054603514;

### Refs:
# 1) Maths 505: Integral of x^2/sin(x) from zero to pi/2
#    https://www.youtube.com/watch?v=g4aKTQyETZw
# 2) Maths 505: A cool integral for Apery's constant
#   (ζ(3)): int 0 to 1 (x(1-x))/sin(πx)
#    https://www.youtube.com/watch?v=TNAtG8gXWuU


### I( 1 / sin(x) )
integrate(function(x) 1 / sin(x) - 1/x, 0, pi/2)
- log(pi/4)

###
integrate(function(x) x / sin(x), 0, pi/2)
2*Catalan

###
integrate(function(x) x^2 / sin(x), 0, pi/2)
2*pi*Catalan - 7/2*pracma::zeta(3)

###
integrate(function(x) x^3 / sin(x), 0, pi/2)
(pracma::psi(3, 3/4) - pracma::psi(3, 1/4)) / 128 + 3/2*pi^2*Catalan

###
integrate(function(x) x^4 / sin(x), 0, pi/2)
3/2*31*zeta(5) + pi^3*Catalan +
	+ pi*(pracma::psi(3, 3/4) - pracma::psi(3, 1/4)) / 64;


### Composed: on [0, pi]
# - some terms cancel out;
integrate(function(x) x*(pi - x) / sin(x), 0, pi)
7*zeta(3)

###
integral(function(x) x*(pi - x)*(pi + x) / sin(x), -pi, pi)
3*7*pi*zeta(3)

###
integral(function(x) x^2*(pi - x)*(pi + x) / sin(x), 0, pi)
2*7*pi^2*zeta(3) - 3*31*zeta(5)

###
integral(function(x) x^3*(pi - x)*(pi + x) / sin(x), -pi, pi)
7*7*pi^3*zeta(3) - 15*31*pi*zeta(5)

###
integral(function(x) x^4*(pi - x)*(pi + x) / sin(x), 0, pi)
9/2*7*pi^4*zeta(3) - 39/2*31*pi^2*zeta(5) + 45/2*127*zeta(7)

###
integral(function(x) x^5*(pi - x)*(pi + x) / sin(x), -pi, pi)
11*7*pi^5*zeta(3) - 90*31*pi^3*zeta(5) + 315/2*127*pi*zeta(7)


########
### TAN:
integrate(function(x) x / tan(x), 0, pi/2)
pi*log(2)/2
#
integrate(function(x) 2*x * (1/tan(x) - tan(x)), 0, pi/4)

# =>
integrate(function(x) x / tan(x), 0, pi/4)
pi*log(2)/8 + Catalan/2


### on [0, pi/3]
# - see file: Integrals.Log.Trig.R;
integrate(\(x) x / tan(x), 0, pi/3)
pi*log(3)/6 - sqrt(3)/(8*6^2) *
	( pracma::psi(1, 1/12) - pracma::psi(1, 11/12) +
	- pracma::psi(1, 5/12) + pracma::psi(1, 7/12) +
	- 5 * pracma::psi(1, 1/6) + 5 * pracma::psi(1, 5/6) +
	- 3 * pracma::psi(1, 2/6) + 3 * pracma::psi(1, 4/6));


###
integrate(function(x) x*(pi-x)*(pi+x) / tan(x), 0, pi)
- 3/2 * pi * zeta(3)

###
integrate(function(x) x*(pi-x)*(pi+x) / tan(x), 0, pi/2)
9/16 * pi * zeta(3) + 3/8 * pi^3 * log(2)

###
integrate(function(x) x^2 / tan(x), 0, pi/2)
- 7/8 * zeta(3) + 1/4 * pi^2 * log(2)

###
integrate(function(x) x^3 / tan(x), 0, pi/2)
- 9/16 * pi * zeta(3) + 1/8 * pi^3 * log(2)

###
integrate(function(x) x^4 / tan(x), 0, pi/2)
93/32 * zeta(5) - 9/16 * pi^2 * zeta(3) + 1/16 * pi^4 * log(2)


### Derived
# on [0, pi]

###
integrate(function(x) x*(pi - x) / sin(x), 0, pi)
7*zeta(3)

###
integrate(function(x) x^2*(pi - x) / sin(x), 0, pi)
1/2*7*pi*zeta(3)

###
integral(function(x) x^3*(pi - x) / sin(x), 0, pi)
3/2*7*pi^2*zeta(3) - 3*31*zeta(5)

###
integral(function(x) x^4*(pi - x) / sin(x), 0, pi)
4/2*7*pi^3*zeta(3) - 9/2*31*pi*zeta(5)

###
integral(function(x) x^5*(pi - x) / sin(x), 0, pi)
5/2*7*pi^4*zeta(3) - 30/2*31*pi^2*zeta(5) + 45/2*127*zeta(7)


### Derived
# on [0, pi/2]

###
integrate(function(x) x*(pi - x) / sin(x), 0, pi/2)
7/2*zeta(3)

###
integrate(function(x) x^2 / sin(x), 0, pi/2)
2*pi*Catalan - 7/2*zeta(3)

###
integrate(function(x) x^2*(pi - x) / sin(x), 0, pi/2)
1/2*pi^2*Catalan - 7/2*pi*zeta(3) +
	- (pracma::psi(3, 3/4) - pracma::psi(3, 1/4)) / 128;

###
integrate(function(x) x^3*(pi - x) / sin(x), 0, pi/2)
3/2*7*pi^2*zeta(3) - 3*31*zeta(5) - integrate(function(x) x*(pi - x)^3 / sin(x), 0, pi/2)$value
- 3/2*31*zeta(5) + 1/2*pi^3*Catalan +
	- pi*(pracma::psi(3, 3/4) - pracma::psi(3, 1/4)) / 128;

###
integrate(function(x) x^2*(pi - x)^2 / sin(x), 0, pi/2)
- 1/2*7*pi^2*zeta(3) + 3/2*31*zeta(5)


#################

### Higher Powers

### I( x^p / sin(x)^2 )
# Maths 505: A surprisingly fascinating integral
# https://www.youtube.com/watch?v=8DWa0zIM9lY
# - see also file: Integrals.Log.Trig.R;


### I( x^2 / sin(x)^2 )
integrate(function(x) x^2 / sin(x)^2, 0, pi/2)
pi*log(2)


### I( x^3 / sin(x)^2 )
integrate(function(x) x^3 / sin(x)^2, 0, pi/2)
3/4*pi^2*log(2) - 3*7/8 * zeta(3)


### I( x^4 / sin(x)^2 )
integrate(function(x) x^4 / sin(x)^2, 0, pi/2)
1/2*pi^3*log(2) - 9/4 * pi * pracma::zeta(3)


### Varia:

# By Parts:
integrate(function(x) x^2 / tan(x)^2, 0, pi/2)
pi*log(2) - pi^3/24

# By Parts: => 1/2 * I( x^2 * (1 + 1/tan(x)^2) )
integrate(function(x) x / tan(x), 0, pi/2)
pi*log(2)/2

# By Parts:
integrate(function(x) log(sin(x)), 0, pi/2)
- pi*log(2)/2

# - for other intervals, see file:
#   Integrals.Log.Trig.R;


#######################
#######################

### ATAN

### I( atan(x^p) / x )
integrate(\(x) atan(x) / x, 0, 1)
Catalan
#
integrate(\(x) atan(x^2) / x, 0, 1)
Catalan / 2
#
p = sqrt(5)
integrate(\(x) atan(x^p) / x, 0, 1)
Catalan / p

### I( x^p * atan(x^n) )
n = sqrt(5)
integrate(\(x) atan(x^n), 0, 1)
pi/4 - n * int.FrU01(2*n, n)

###
n = sqrt(5); p = sqrt(3);
integrate(\(x) x^p * atan(x^n), 0, 1)
pi/(4*(p+1)) +
	- (digamma(((p+1)/n + 3)/4) - digamma(((p+1)/n + 1)/4)) / (4*(p+1));
pi/(4*(p+1)) - n / (p+1) * int.FrU01(2*n, n + p)


###
# see file: Integrals.Trig.Tan.R;
integrate(\(x) x * atan(x) / (x^2 + 1), 0, 1)
Catalan / 2 - 1/8 * pi*log(2)

###
integrate(\(x) atan(x) / (x^2 + 1), 0, 1)
pi^2 / 32


### I( x^p * atan(x^n) * log(x) )
n = sqrt(3); p = sqrt(5);
# works only for n > 0;
# see also file: Integrals.Trig.Tan.R;
integrate(\(x) x^p * atan(x^n) * log(x), 0, 1)
- pi/(4*(p+1)^2) +
	+ (digamma(((p+1)/n + 3)/4) - digamma(((p+1)/n + 1)/4)) / (4*(p+1)^2) +
	- (pracma::psi(1, ((p+1)/n + 3)/4) - pracma::psi(1, ((p+1)/n + 1)/4)) / (16*(p+1)*n);

