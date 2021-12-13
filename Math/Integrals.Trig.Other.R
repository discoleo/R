########################
###
### Leonard Mada
### [the one and only]
###
### Integrals:
### Trigonometric: Other
###
### draft v.0.1b



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


#########

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
# Note: needs to be multiplied by: pi/2 * R = (pi/2)^2
Re(line_integral(\(x) exp(2i*x)*log(sin(pi/2*exp(1i*x))), c(0, pi/2)))
- pi / 4
#
line_integral(\(x) exp(2i*x)*log(sin(pi/2*exp(1i*x))), c(0, 2*pi))
pi # * (pi/2)^2 => pi^3/4

###
integrate(function(x) x*log(sin(x)/cos(x)), 0, pi/2, subdivisions=1024)
integrate(function(x) (x-pi/2)*log(sin(x)/cos(x)), 0, pi/2, subdivisions=1024)

# Note: I == 0
integrate(function(x) log(sin(x)/cos(x)), 0, pi/2);

integrate(function(x) log(cos(x)/sin(x)), 0, pi/4)
integrate(function(x) - log(x)/(1+x^2), 0, tan(pi/4))
# TODO: How?

