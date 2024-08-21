

### Integrals: Log(Poly)


################
################

### I( log(x^2 + 2*sin(a)*x + 1) / x ) on [0,1]
# Maths 505: A beautifully unusual integral calculation
# https://www.youtube.com/watch?v=Na5kpJY0s9g

a = asin(sqrt(3)/4)
integrate(\(x) log(x^2 - 2*sin(a)*x + 1) / x, 0, 1)
- a^2/2 - pi/2*a + pi^2 / 24

###
a = asin(1/pi)
integrate(\(x) log(x^2 - 2*sin(a)*x + 1) / x, 0, 1)
- a^2/2 - pi/2*a + pi^2 / 24

###
a = asin(1/pi)
integrate(\(x) log(x^2 + 2*sin(a)*x + 1) / x, 0, 1)
- a^2/2 + pi/2*a + pi^2 / 24

###
a = asin(pi + 0i)
integrate(\(x) log(x^2 + 2*Re(sin(a))*x + 1) / x, 0, 1)
- a^2/2 + pi/2*a + pi^2 / 24


### Special: x^n + 1
# see file: Integrals.Log.Fractions.R
n = sqrt(5); p = sqrt(3);
integrate(function(x) log(x^n + 1) / x^(p+1), 0, 1)
(digamma(((n-p)/n + 1)/2) - digamma((n-p)/(2*n))) / (2*p) - log(2)/p;

# p -> 0
p = 1E-4
(digamma(1 - p/(2*n)) - digamma(1/2 - p/(2*n))) / (2*p) - log(2)/p;
pi^2/(12*n)

###
integrate(\(x) log(x^3 + 1) / x, 0, 1)
a = 2*pi/3
- a^2 / 2 + pi^2 / 4

# Deriv:
- (pi/2 - a)^2/2 + pi/2*(pi/2 - a) + pi^2 / 24 + pi^2/12

###
n = sqrt(11)
integrate(\(x) log(x^n + 1) / x, 0, 1)
pi^2 / (12*n)

