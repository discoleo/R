########################
###
### Leonard Mada
### [the one and only]
###
### Integrals:
### Trigonometric: Other
###
### draft v.0.1e



#####################

### TRIG( TRIG )

### I( sin(tan(x)) * cos(x) / x )
# Maths 505: A beautiful nested trig integral
# https://www.youtube.com/watch?v=IHjGVTeFdFA
# Lobachevsky Transform => Subst => Feynman;

# Note: fails;
integrate(\(x) sin(tan(x)) * cos(x) / x, 0, Inf)
integrate(\(x) sin(tan(x)) / tan(x), 0, pi/2)
integrate(\(x) sin(x) / x / (x^2+1), 0, Inf, rel.tol=1E-9, subdivisions = 2000)
pi/2 * (1 - exp(-1))

### Gen: I( sin(k*tan(x)) * cos(x) / x )
k = 1/3;
integrate(\(x) sin(k*tan(x)) * cos(x) / x, 0, Inf)
integrate(\(x) sin(k*tan(x)) / tan(x), 0, pi/2)
integrate(\(x) sin(k*x) / x / (x^2+1), 0, Inf, rel.tol=1E-9, subdivisions = 2000)
pi/2 * (1 - exp(-k));


#
tmp = lapply(seq(2000), \(id, subdiv = 2025) {
	tmp = integrate(\(x) sin(tan(x)) * cos(x) / x, (id-1) * pi/2, id * pi/2,
		rel.tol=1E-9, subdivisions = subdiv);
	class(tmp) = "list"; tmp$call = NULL;
	as.data.frame(tmp);
})
tmp = do.call(rbind, tmp);
sum(tmp$value); table(cut(tmp$subdivisions, c(0, 100, 200, 500, 10000)))


### I( Trig(Trig) )
# Maths 505: The coolest trig integral!
# https://www.youtube.com/watch?v=4pvK0DTPPg8
# Note: series expansion of cos(t);
# => Beta function & Duplication formula for Gamma & Bessel J;

integrate(\(x) cos(cos(x)) + sin(sin(x)), 0, pi)
# TODO

### I( cos(cos(x)) )
integrate(\(x) cos(cos(x)), 0, pi)
integrate(\(x) 2 * cos(cos(x)), 0, pi/2)
pi * besselJ(1, 0)

### I( sin(sin(x)) )
integrate(\(x) sin(sin(x)), 0, pi)
integrate(\(x) 2 * sin(sin(x)), 0, pi/2)
# Struve function:
id = seq(0, 100)
pi * sum((-1)^id * (1/2)^(2*id+1) / factorial(id+1/2)^2)
# TODO: closed formula;


### Struve Function
# https://en.wikipedia.org/wiki/Struve_function
x = 1/3; n = 1;
id = seq(0, 100)
#
sum((-1)^id * (x/2)^(2*id+n+1) / gamma(id+3/2) / gamma(id+3/2+n))


### Bessel J
x = 1/3; n = 1;
id = seq(0, 100)
#
besselJ(x, n)
sum((-1)^id * (x/2)^(2*id+n) / factorial(id) / factorial(id+n))

# Note:
besselJ(1, 1/2) / besselJ(1, -1/2) # ==
tan(1);


#####################
#####################

### I( atan(1/x0^n) )
# on [0, Inf]


### Ref:
# Michael Penn: An inverse tangent integral.
# https://www.youtube.com/watch?v=oh0qCeV7o9k
# https://www.youtube.com/watch?v=7wiybMkEfbc


# TODO:
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
integrate(\(x) log(cos(x)), 0, pi/2, subdivisions=1024)
integrate(\(x) log(sin(x)), 0, pi/2, subdivisions=1024)
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

### x * sin(k*x) / (a^2 + 1 - 2*a*cos(k*x))

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
pi/(k*a)*log(a^2 + 1) + 2/a/k^2*log(a)*atan(a) +
	+ 2/a/k^2*integrate(function(z) log(z) / (z^2 + 1), 1/a, Inf)$value;

# Derivation:
line_integral(function(z) z/(a - exp(-1i*k*z)), pi+c(0, 30i)) +
	+ line_integral(function(z) z/(a - exp(-1i*k*z)), -pi+c(30i, 0))
pi/(k*a)*log(a^2 + 1) +
	+ 2*integrate(function(z) z / (a^2*exp(-k*z) + exp(k*z)), 0, Inf)$value;
pi/(k*a)*log(a^2 + 1) +
	+ 2/k^2*integrate(function(z) log(z) / (z^2 + a^2), 1, Inf)$value;


###
a = sqrt(3) + 5^(1/3)
k = 1/2; # FIXED value!
integrate(function(x) x*sin(k*x) / (a^2 + 1 - 2*a*cos(k*x)), 0, pi)
- pi*log(a)/(k^2*a) + pi*log(a^2 + 1)/(2*k*a) + log(a)*atan(a)/(k^2*a) +
	+ 1/a/k^2*integrate(function(z) log(z) / (z^2 + 1), 1/a, Inf)$value;


#####################
#####################

### I( sin(x) * sin(1/x) )
# Upper = Inf: numerical instability;
pracma::integral(\(x) sin(x) * sin(1/x), 0, 320000)
pi/2 * Rmpfr::jn(1,2)

