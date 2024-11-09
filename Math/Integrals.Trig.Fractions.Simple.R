##################
##
## Integrals: Trig
## Polynomial Fractions
## Basic Type
##
## Leonard Mada


##################

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
################

### on [0, pi/2]

###
integrate(\(x) sin(x) / x, 0, pi/2)
pi/2 - integrate(\(x) x * exp(-pi/2 * x) / (x^2+1), 0, Inf, rel.tol=1E-8)$value
# TODO: ???
# Note: Walpha reverses the order of the integrals;

