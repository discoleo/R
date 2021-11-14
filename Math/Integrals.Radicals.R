########################
###
### Leonard Mada
### [the one and only]
###
### Infinite Sums: Fractions
### Radicals
###
### draft v.0.1c


###########################

###########################
### I ( sqrt(x^n + 1) ) ###
###########################

### n = 3
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
# use cos() instead of sin():
r = line_integral(function(x)  sqrt(cos(3*x/2)) * sin(7*x/4), c(0, pi/3)) +
line_integral(function(x) - sqrt(-cos(3*x/2)) * cos(7*x/4), c(pi/3, 2*pi/3));
r * 2*sqrt(2) / 3
#
r = line_integral(function(x)  sqrt(cos(6*x)) * sin(7*x), c(0, pi/12)) +
line_integral(function(x) - sqrt(-cos(6*x)) * cos(7*x), c(pi/12, pi/6));
r * sqrt(2) * 8/3
- (rh1 + rh2) * 1i/3


#######################
#######################

### n = 5

r1 = line_integral(function(x)  sqrt(exp(5i*x) + 1) * exp(1i*x), c(0,  2*pi/5))
r2 = line_integral(function(x)  sqrt(exp(5i*x) + 1) * exp(1i*x), c(0, -2*pi/5))
r3 = line_integral(function(x)  sqrt(exp(5i*x) + 1) * exp(1i*x), c(0,  4*pi/5))
r4 = line_integral(function(x)  sqrt(exp(5i*x) + 1) * exp(1i*x), c(0, -4*pi/5))
- (r1 + r2 + r3 + r4) * 1i/5

integrate(function(x) sqrt(x^5 + 1), lower=0, upper=1)


#######################
#######################

### n = 7
n = 7

id = seq(n %/% 2); id = c(-id, id);
r = sapply(id, function(id) {
	line_integral(function(x)  sqrt(exp(n * 1i*x) + 1) * exp(1i*x), c(0,  2*id*pi/n))
})
sum(r) * -1i/n

integrate(function(x) sqrt(x^n + 1), lower=0, upper=1)

