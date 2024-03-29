


### Helper

Catalan = 0.915965594177219015054603514;


int.FrU01 = function(n, p=0) {
	(digamma(((p+1)/n + 1)/2) - digamma((p+1)/(2*n))) / (2*n);
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

### I( sin(x) / (x^2 + 1) )
# Various variants:

# 1) qncubed3:  TWO WAYS to destroy this INSANE integral feat. @maths_505
#    https://www.youtube.com/watch?v=Ry5q4NsZDx0
# 2) Maths 505: INSANE integral solved using 2 different methods
#    (feat. Feynman's technique and Lobachevsky)
#    https://www.youtube.com/watch?v=GMcCiR0y3wg


###
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
integrate(\(x) sin(k*x) / (x^2 + 1), 0, Inf, subdivisions=4096*2, rel.tol=1E-5)
integrate(\(x) sapply(x, \(x) {
	1/(4*pi) * Im(sin(k*x) * (pracma::psi((x + 1i)/(2*pi)) - pracma::psi((x - 1i)/(2*pi))))
	}), 0, 2*pi)


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

### I( sin(tan(x)) / x )
# takes very long: upper = Inf
pracma::integral(\(x) sin(tan(x)) / x, 0, 2000)
pracma::integral(\(x) sin(x) / (x*(x^2 + 1)), 0, 10000)
pi/2*(1 - exp(-1))


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
# 3. Maths 505: The generalised Dirichlet integral:
#    Integral of (sinx)^n/x^n from zero to infinity
#    https://www.youtube.com/watch?v=lLDiIf_BB-I

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

# - some terms cancel out when using: x*(pi - x);
# - for intuition of underlying principle, see e.g:
#   blackpenredpen: I saw this NTU entrance exam integral on Dcard
#   and I just had to integrate it
#   https://www.youtube.com/watch?v=PWGJGnji5Nk

### I( 1 / sin(x) )
integrate(\(x) 1 / sin(x) - 1/x - 1/(pi-x), 0, pi)
- 2*log(pi/2)

### I( x / sin(x) )
integrate(\(x) x / sin(x) - pi / (pi-x), 0, pi)
- pi*log(pi/2)

### I( x^2 / sin(x) )
integrate(\(x) x^2 / sin(x) - pi^2 / (pi-x), 0, pi)
- 7*pracma::zeta(3) - pi^2*log(pi/2)

###
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


### on Other Intervals

###
integrate(function(x) 1 / sin(x) - 1/x, 0, pi/3)
log(3)/2 - log(pi/2)


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
integrate(\(x) cos(x) / sin(x)^2 - 1/x^2, 0, pi/3)
-2/sqrt(3) + 3/pi

### [solvable on any interval]
integrate(\(x) x * cos(x) / sin(x)^2 - 1/x, 0, pi/3)
log(3)/2 - log(pi/2) - pi/3 / sin(pi/3) + 1


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


### Helper:
k = 1/5;
fk = function(x, k) log((1 - cos(k*x)) / (1 + cos(k*x))) / (2*k);
integrate(\(x) 1 / sin(k*x), pi/5, pi/3)
fk(pi/3, k) - fk(pi/5, k)


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

###
integrate(function(x) 1 / sin(x)^2 - 1/x^2, 0, pi/2)
2/pi

### I( x / sin(x)^2 )
integrate(function(x) x / sin(x)^2 - 1/x, 0, pi/2)
1 - log(pi/2)

### I( x^2 / sin(x)^2 )
integrate(function(x) x^2 / sin(x)^2, 0, pi/2)
pi*log(2)


### I( x^3 / sin(x)^2 )
integrate(function(x) x^3 / sin(x)^2, 0, pi/2)
3/4*pi^2*log(2) - 3*7/8 * zeta(3)


### I( x^4 / sin(x)^2 )
integrate(function(x) x^4 / sin(x)^2, 0, pi/2)
1/2*pi^3*log(2) - 9/4 * pi * pracma::zeta(3)


### I( 1 / sin(x)^3 )
integrate(function(x) 1 / sin(x)^3 - 1/x^3 - 1/(2*x), 0, pi/2)
log(2)/2 - 1/12 + (2/pi^2 - log(pi/2)/2);

# Note:
# I( 1 / sin(x)^(2*n+1) ) can be computed on arbitrary intervals;
### I( 1 / sin(x)^3 )
up = pi/5;
integrate(function(x) 1 / sin(x)^3 - 1/x^3 - 1/(2*x), 0, up)
(log( (1 - cos(up)) / (1 + cos(up)) ) +
	+ 1/(cos(up) - 1) + 1/(cos(up) + 1) + 2*log(2) - 1/3) / 4 +
	+ (1/(2*up^2) - log(up)/2);

# Lim: x -> 0
x = 1E-3; # but numerical issues!
(-log(1 - cos(x)) - 1/(cos(x)-1)) / 2 - 1/(x^2) + log(x)
log(2)/2 + 1/12


### on [0, pi/4]

### I( 1 / sin(x)^2 )
integrate(function(x) 1 / sin(x)^2 - 1/x^2, 0, pi/4)
4/pi - 1
#
integrate(function(x) 1 / sin(x)^2, pi/4, pi/2)
1;

### I( x / sin(x)^2 )
integrate(function(x) x / sin(x)^2 - 1/x, 0, pi/4)
1 - pi/4 - 1/2 * log(2) - log(pi/4)
# Note: log(pi/4) is from 1/x;
integrate(function(x) x / sin(x)^2, pi/4, pi/2)
pi/4 + log(2)/2

### I( x^2 / sin(x)^2 ) on [0, pi/4]
integrate(function(x) x^2 / sin(x)^2, 0, pi/4)
- pi^2 / 16 + pi/4 * log(2) + Catalan


### Arbitrary Interval:

### Base: I( 1 / tan(k*x) )
k = pi/7
integrate(\(x) 1 / tan(k*x) - 1/(k*x), 0, 1)
log(sin(k)) / k - log(k) / k;


### Base: I( 1 / sin(k*x)^3 )
k = pi/7
integrate(\(x) 1 / sin(k*x)^3 - 1/(k^3*x^3) - 1/(2*k*x), 0, 1)
(log( (1 - cos(k)) / (1 + cos(k)) ) +
	+ 1/(cos(k) - 1) + 1/(cos(k) + 1) + 2*log(2) - 1/3) / (4*k) +
	+ 1/(2*k^3) - log(k)/(2*k);

### Base: I( cos(k*x) / sin(k*x)^3 )
k = pi/7
integrate(\(x) cos(k*x) / sin(k*x)^3 - 1/(k^3*x^3), 0, 1)
- sin(k)^(-2) / (2*k) + 1/(2*k^3) + 1/(6*k);


# =>

### I( x / sin(k*x)^2 )
k = pi/7
integrate(\(x) x / sin(k*x)^2 - 1/(k^2*x), 0, 1)
- cos(k) / (k * sin(k)) + log(sin(k)) / k^2 - log(k) / k^2 + 1/k^2;


### I( x * cos(k*x) / sin(k*x)^4 )
k = pi/2 + 1/11;
integrate(\(x) x * cos(k*x) / sin(k*x)^4 - 1/(k^4*x^3) - 1/(6*k^2*x), 0, 1)
- (2/sin(k) + 2*(cos(k)^2 + 1)/sin(k)^3) / (12*k) +
	+ (log( (1 - cos(k)) / (1 + cos(k)) ) +
		- 2*cos(k)/sin(k)^2 + 2*log(2) - 1/3) / (12*k^2) +
	+ 1/(2*k^4) + 1/(6*k^2) - log(k)/(6*k^2);


### I( x / sin(k*x)^4 )
k = pi/7
integrate(\(x) x / sin(k*x)^4 - 1/(k^4*x^3) - 2/(3*k^2*x), 0, 1)
- cos(k) / (3*k * sin(k)^3) - 1 / (6*k^2 * sin(k)^2) - 2*cos(k) / (3*k * sin(k)) +
	+ (12*log(sin(k)) - 12*log(k) + 9/k^2 + 13) / (18*k^2);


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
pi/(4*(p+1)) - n / (p+1) * int.FrU01(2*n, n + p);


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

