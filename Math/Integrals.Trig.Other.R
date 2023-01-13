########################
###
### Leonard Mada
### [the one and only]
###
### Integrals:
### Trigonometric: Other
###
### draft v.0.1d



#####################
#####################

### I( atan(1/x0^n) )
# on [0, Inf]


### Ref:
# Michael Penn: An inverse tangent integral.
# https://www.youtube.com/watch?v=oh0qCeV7o9k
# https://www.youtube.com/watch?v=7wiybMkEfbc


# [work in progress]

n = 3
integrate(function(x) atan(1/x)^n, 0, Inf, subdivisions=1024)
integrate(function(x) n*x^(n-1) / tan(x), 0, pi/2, subdivisions=1024)


n = 5
integrate(function(x) atan(1/x)^n, 0, Inf, subdivisions=1024)
integrate(function(x) n*x^(n-1) / tan(x), 0, pi/2, subdivisions=1024)


#####################
#####################

### n = 2
integrate(function(x) log(cos(x)), 0, pi/2, subdivisions=1024)
integrate(function(x) log(sin(x)), 0, pi/2, subdivisions=1024)
- (pi/2)*log(2)


### n = 3
integrate(function(x) - x*log(cos(x)), 0, pi/2, subdivisions=1024)
integrate(function(x) x*log(sin(x)), 0, pi/2, subdivisions=1024)$value + (pi/2)^2*log(2)
# unfortunately I[cos] = - I[sin] + ...;

integrate(function(x) 4 * x*log(sin(2*x)), 0, pi/2, subdivisions=1024)
integrate(function(x) x*log(sin(x)), 0, pi, subdivisions=1024)
- pi^2*log(2) / 2

integrate(function(x) x*log(sin(x)), 0, pi/2, subdivisions=1024)$value + 
	integrate(function(x) (x+pi/2)*log(cos(x)), 0, pi/2, subdivisions=1024)$value
- pi^2*log(2) / 2
# but: I[sin] + I[cos] = constant!

### Contour: Key-hole
line_integral(\(x) x*log(sin(1i * x)), c(0, pi/2)) +
	(pi/2)^2*1i * line_integral(\(x) exp(2i*x)*log(sin(pi/2*exp(1i*x))), c(0, pi/2))
integrate(function(x) - x*log(sin(x)), 0, pi/2, subdivisions=1024)
# Note: needs to be multiplied by: pi/2 * R * 1i = (pi/2)^2 * 1i;
Re(line_integral(\(x) exp(2i*x)*log(sin(pi/2*exp(1i*x))), c(0, pi/2)))
- pi / 4
#
line_integral(\(x) exp(2i*x)*log(sin(pi/2*exp(1i*x))), c(0, 2*pi))
pi # * (pi/2)^2 * 1i => pi^3/4 * 1i

###
integrate(function(x) x*log(sin(x)/cos(x)), 0, pi/2, subdivisions=1024)
integrate(function(x) (x-pi/2)*log(sin(x)/cos(x)), 0, pi/2, subdivisions=1024)

# Note: I == 0
integrate(function(x) log(sin(x)/cos(x)), 0, pi/2);

integrate(function(x) log(cos(x)/sin(x)), 0, pi/4)
integrate(function(x) - log(x)/(1+x^2), 0, tan(pi/4))
# TODO: How?


########################
########################

### x*sin(k*x) / (a^2 + 1 - 2*a*cos(k*x))

### Ref:
# qncubed3: Complex Analysis: Viewer Suggests INSANE Integral
# https://www.youtube.com/watch?v=Q4il1GoCJVE

# Generalization:

###
# Case: a >= 1;
a = sqrt(3)
k = 3; # odd
integrate(function(x) x*sin(k*x) / (a^2 + 1 - 2*a*cos(k*x)), 0, pi)
pi*log((a+1)/a) / (k*a)
pi*log((a - exp(-1i*pi*k))/a) / (k*a)

###
a = sqrt(3)
k = 2; # even
integrate(function(x) x*sin(k*x) / (a^2 + 1 - 2*a*cos(k*x)), 0, pi)
pi*log((a-1)/a) / (k*a)
pi*log((a - exp(-1i*pi*k))/a) / (k*a)

###
a = sqrt(3)
k = 4; # even
integrate(function(x) x*sin(k*x) / (a^2 + 1 - 2*a*cos(k*x)), 0, pi)
pi*log((a-1)/a) / (k*a)
pi*log((a - exp(-1i*pi*k))/a) / (k*a)
#
integrate(function(x) x*sin(x) / (a^2 + 1 - 2*a*cos(x)), 0, pi*k)
k*pi*log((a-1)/a) / a

###
# TODO: k = non-integer;


### Gen 2:
b = c(5, sqrt(3))
k = 3; # odd
integrate(function(x) x*sin(k*x) / (b[1] - b[2]*cos(k*x)), 0, pi)
bb = b[1]/b[2]; a = bb + sqrt(bb^2 - 1)
2*pi*log((a+1)/a) / (k*b[2])

###
b = c(5, sqrt(3))
k = 4; # even
integrate(function(x) x*sin(k*x) / (b[1] - b[2]*cos(k*x)), 0, pi)
bb = b[1]/b[2]; a = bb + sqrt(bb^2 - 1)
2*pi*log((a-1)/a) / (k*b[2])


### Analysis

k = 1/2; # FIXED value!
a = sqrt(3)
# Factor 1i omitted:
line_integral(function(z) z/(a - exp(-1i*k*z)), pi+c(0, 30i)) +
	+ line_integral(function(z) z/(a - exp(-1i*k*z)), -pi+c(30i, 0))
pi/(k*a)*log(a^2 + 1) + 2/a/k^2*log(a)*(pi/2 - atan(1/a)) +
	+ 2/a/k^2*integrate(function(z) log(z) / (z^2 + 1), 1/a, Inf)$value;

# Derivation:
line_integral(function(z) z/(a - exp(-1i*k*z)), pi+c(0, 30i)) +
	+ line_integral(function(z) z/(a - exp(-1i*k*z)), -pi+c(30i, 0))
pi/(k*a)*log(a^2 + 1) +
	+ 2*integrate(function(z) z / (a^2*exp(-k*z) + exp(k*z)), 0, Inf)$value;
pi/(k*a)*log(a^2 + 1) +
	+ 2/k^2*integrate(function(z) log(z) / (z^2 + a^2), 1, Inf)$value;

