

# Note:
# erf(1) = 2 * (pnorm(1, sd = 1/sqrt(2)) - 0.5);
# erf(1) = 2 * (pnorm(sqrt(2)) - 0.5);
# I(exp(-x^2), 0, 1) = erf(1) * sqrt(pi)/2;


################

### I( exp(-x^n) / (x^n + 1) ) on [0, Inf]
# Maths 505: A crazy yet perfect integral
# https://www.youtube.com/watch?v=regnh-HOR0c
# - the sub-integral;
# - using Feynman's technique;

### I( exp(-x^2) / (x^2 + 1) ) on [0, Inf]
integrate(\(x) exp(-x^2) / (x^2 + 1), 0, Inf)
exp(1) * (pi/2 - gamma(1/2) * integrate(\(x) exp(-x^2), 0, 1)$value)
exp(1) * (pi/2 - gamma(1/2) * (pnorm(sqrt(2)) - 0.5)*sqrt(pi))


###
integrate(\(x) exp(-x^3) / (x^3 + 1), 0, Inf)
exp(1) * (pi/3/sin(pi/3) - gamma(1/3) * integrate(\(x) x * exp(-x^3), 0, 1)$value)
exp(1) * (pi/3/sin(pi/3) - gamma(1/3)/2 * integrate(\(x) exp(-x^1.5), 0, 1)$value)


###
integrate(\(x) exp(-x^4) / (x^4 + 1), 0, Inf)
exp(1) * (pi/4/sin(pi/4) - gamma(1/4) * integrate(\(x) x^2 * exp(-x^4), 0, 1)$value)
exp(1) * (pi/4/sin(pi/4) - gamma(1/4)/3 * integrate(\(x) exp(-x^(4/3)), 0, 1)$value)

### Generalization:

###
n = sqrt(3)
integrate(\(x) exp(-x^n) / (x^n + 1), 0, Inf)
int = integrate(\(x) x^(-1/n) * exp(-x), 0, 1)$value
exp(1) * (pi/sin(pi/n) - int*gamma(1/n))/n

### [stable]
n = sqrt(3) - sqrt(2)
integrate(\(x) exp(-x^n) / (x^n + 1), 0, Inf)
int = integrate(\(x) n*x^(n-2) * exp(-x^n), 1, Inf)$value
exp(1) * int * gamma(1/n) / n

# [unstable/divergent] only: n > 1
# int = integrate(\(x) x^(-1/n) * exp(-x), 0, 1)$value
int = integrate(\(x) n*x^(n-2) * exp(-x^n), 0, 1)$value
exp(1) * (pi/sin(pi/n) - int*gamma(1/n)) / n


### I( x^p * exp(-x^n) / (x^n + 1) )
p = sqrt(3); n = sqrt(11);
integrate(\(x) x^p * exp(-x^n) / (x^n + 1), 0, Inf)
# numerical issues when p+1 >= n;
int = integrate(\(x) x^(-(p+1)/n) * exp(-x), 0, 1)$value
exp(1) * (pi/sin(pi*(p+1)/n) - int*gamma((p+1)/n)) / n


### Special Cases:
### p == n
n = sqrt(11);
integrate(\(x) x^n * exp(-x^n) / (x^n + 1), 0, Inf)
int = integrate(\(x) x^(-1/n) * exp(-x), 0, 1)$value
gamma(1/n) / n - exp(1) * (pi/sin(pi/n) - int*gamma(1/n)) / n

# TODO:
n = sqrt(11); p = n-1;
integrate(\(x) x^p * exp(-x^n) / (x^n + 1), 0, Inf)
# Exponential Integral:
# - exp(1) * Ei(-1) / n;


##################

### I( exp(-x^n) / (x^n - 1) )
n = sqrt(3)
# Note: numerical issues;
integrate(\(x) exp(-x^n) / (x^n - 1), 0, 1 - 1E-7, rel.tol=1E-8)$value +
	integrate(\(x) exp(-x^n) / (x^n - 1), 1 + 1E-7, Inf, rel.tol=1E-8)$value;
# Partial Gamma:
int = integrate(\(x) abs(x)^(-1/n) * exp(-x), 0, -1)$value;
int = - n/(n-1) * integrate(\(x) exp(x^(n/(n-1))), 0, 1, rel.tol = 1E-8)$value;
- (pi/tan(pi/n) - int*gamma(1/n)) / (n * exp(1))


#########################

### Exponential Integrals

### exp(1) * Ei(-1)
integrate(\(x) exp(-x) / (x + 1), 0, Inf)
integrate(\(x) exp(1 - exp(x)), 0, Inf)
exp(1) * pracma::expint(1)

integrate(\(x) exp(1-x)/x - 1 / (x*(x+1)), 1, Inf)
exp(1) * pracma::expint(1) - log(2)

integrate(\(x) exp(1-x)/x - 1 / (x*(x+1)), 1/2, Inf)
exp(1) * pracma::expint(1/2) - log(3)

integrate(\(x) exp(x)/x - 1/x/(x+1), 0, 1/2)
- Re(pracma::expint(-1/2)) - Euler + log(3)

integrate(\(x) exp(1/2-abs(x))/(1/2+x) - 1 / ((1/2+x)*(x+3/2)), -1/2, Inf)
exp(1) * pracma::expint(1/2) - Re(pracma::expint(-1/2)) - Euler

integrate(\(x) exp(1/2-abs(x))/(1/2+x) + 1 / ((x+1/2)*(x-1/2)), -Inf, -1/2)
Euler

