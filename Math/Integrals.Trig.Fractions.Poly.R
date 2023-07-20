


### Helper

Catalan = 0.915965594177219015054603514;


int.FrU01 = function(n, p=0) {
	(digamma(((p+1)/n + 1)/2) - digamma((p+1)/n/2)) / (2*n);
}
int.FrDU01 = function(n, p=0) {
	# TODO: p != 0
	digamma(1/n) + Euler + log(n);
}


####################
####################

# qncubed3:  TWO WAYS to destroy this INSANE integral feat. @maths_505
# https://www.youtube.com/watch?v=Ry5q4NsZDx0

integrate(function(x) sin(x)^2 / (x^2*(x^2+1)), 0, Inf)
pi/4*(1 + exp(-2))


### Helper:
k = sqrt(3)
integrate(function(x) cos(k*x) / (x^2 + 1), 0, Inf, subdivisions=4096*2, rel.tol=1E-6)
pi/2*exp(-k)


### Gen 1:
k = sqrt(3)
b = sqrt(5)
integrate(function(x) cos(k*x) / (x^2 + b^2), 0, Inf, subdivisions=4096*2, rel.tol=1E-6)
pi/2*exp(-k*b)/b


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

### sin(x) / x^(s+1)
# A RIDICULOUSLY AWESOME INTEGRAL: Ramanujan vs Maths 505
# https://www.youtube.com/watch?v=_VkRvuSxF18

# s in (-1, 1)
s = 1 - sqrt(2)
# FAILS:
integrate(\(x) sin(x) / x^(s+1), 0, Inf, subdivisions=128)
integrate(\(x) sin(x) / x^(s+1), 0, 200000, subdivisions=400000)
pracma::integral(\(x) sin(x) / x^(s+1), 0, 2^16*pi, no_intervals=1024)
- gamma(-s)*sin(pi*s/2)

# s in (-1, 1)
s = 2 - sqrt(2)
integrate(\(x) sin(x) / x^(s+1), 0, 200000, subdivisions=40000)
- gamma(-s)*sin(pi*s/2)

# s = 0
# integrate(\(x) sin(x) / x, 0, Inf, subdivisions=128)
integrate(\(x) sin(x) / x, 0, 200000, subdivisions=400000)
pi/2

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


### TAN:
integrate(function(x) x / tan(x), 0, pi/2)
pi*log(2)/2
#
integrate(function(x) 2*x * (1/tan(x) - tan(x)), 0, pi/4)

# =>
integrate(function(x) x / tan(x), 0, pi/4)
pi*log(2)/8 + Catalan/2

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
# TODO: seems to involve polygamma;
integrate(function(x) x^2*(pi - x) / sin(x), 0, pi/2)

###
integrate(function(x) x^3*(pi - x) / sin(x), 0, pi/2)
3/2*7*pi^2*zeta(3) - 3*31*zeta(5) - integrate(function(x) x*(pi - x)^3 / sin(x), 0, pi/2)$value
1/2*7*pi^2*zeta(3) - 3/2*31*zeta(5) + pi * integrate(function(x) x^2*(pi - x) / sin(x), 0, pi/2)$value
###
integrate(function(x) x^2*(pi - x)^2 / sin(x), 0, pi/2)
- 1/2*7*pi^2*zeta(3) + 3/2*31*zeta(5)


### Square:
integrate(function(x) x^2 / sin(x)^2, 0, pi/2)
pi*log(2)

# By Parts:
integrate(function(x) x^2 / tan(x)^2, 0, pi/2)
pi*log(2) - pi^3/24

# By Parts:
integrate(function(x) x / tan(x), 0, pi/2)
pi*log(2)/2

# By Parts:
integrate(function(x) log(sin(x)), 0, pi/2)
- pi*log(2)/2


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
pi/(4*(p+1)) - n / (p+1) * int.FrU01(2*n, n + p)


###
# see file: Integrals.Trig.Tan.R;
integrate(\(x) atan(x) * x / (x^2 + 1), 0, 1)
Catalan / 2 - 1/8 * pi*log(2)

###
integrate(\(x) atan(x) / (x^2 + 1), 0, 1)
pi^2 / 32

