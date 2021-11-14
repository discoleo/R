########################
###
### Leonard Mada
### [the one and only]
###
### Infinite Sums: Fractions
### Radicals
###
### draft v.0.1a


######################

### I( sqrt(x^3 + 1) )

r1 = line_integral(function(x)  sqrt(exp(3i*x) + 1) * exp(1i*x), c(0,  2*pi/3))
r2 = line_integral(function(x)  sqrt(exp(3i*x) + 1) * exp(1i*x), c(0, -2*pi/3))
r2 = line_integral(function(x) - sqrt(exp(-3i*x) + 1) * exp(-1i*x), c(0, 2*pi/3))
- (r1 + r2) * 1i/3

integrate(function(x) sqrt(x^3 + 1), lower=0, upper=1)

# TODO:
# - understand + 0i discrepancy!
rh1 = line_integral(function(x)   sqrt(2) * sqrt(cos(3*x/2 + 0i)) * exp(7i*x/4), c(0, 2*pi/3))
rh2 = line_integral(function(x) - sqrt(2) * sqrt(cos(3*x/2) + 0i) * exp(-7i*x/4), c(0, 2*pi/3))
- (rh1 + rh2) * 1i/3
# TODO: debug!
line_integral(function(x)   2i* sqrt(2) * sqrt(cos(3*x/2)) * sin(7*x/4), c(0, pi/3)) +
line_integral(function(x)   2i*sqrt(2) * sqrt(-cos(3*x/2)) * sin(7*x/4), c(pi/3, 4*pi/7)) +
line_integral(function(x) - 2i*sqrt(2) * sqrt(-cos(3*x/2)) * sin(7*x/4), c(4*pi/7, 2*pi/3))

