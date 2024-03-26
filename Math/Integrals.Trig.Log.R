


### I( sin(k*x) * log(x) / x^(p+1) )

###
p = sqrt(3) - sqrt(2)
# Up = Inf; numerical issues;
integrate(\(x) sin(x) * log(x) / x^(p+1), 0, 200000, subdivisions=40000)
- gamma(-p) * digamma(-p) * sin(pi*p/2) + pi/2 * gamma(-p) * cos(pi*p/2)


###
p = sqrt(3) - sqrt(2)
k = sqrt(3)
# Up = Inf; numerical issues;
integrate(\(x) sin(k*x) * log(x) / x^(p+1), 0, 100000, subdivisions=40000)
- gamma(-p) * digamma(-p) * sin(pi*p/2) * k^p +
	+ pi/2 * gamma(-p) * cos(pi*p/2) * k^p + gamma(-p) * sin(pi*p/2) * log(k) * k^p;


### Base:
# p in (-1, 1)
# see file: Integrals.Trig.Fractions.Poly.R;
p = 2 - sqrt(2)
k = sqrt(3)
integrate(\(x) sin(k*x) / x^(p+1), 0, 200000, subdivisions=40000)
- gamma(-p) * sin(pi*p/2) * k^p;


#################
#################

### I( cos(log(x)) / ... )

# Maths 505: A nice integral result
# https://www.youtube.com/watch?v=b-5f66DSD7g

# Generalization:

p = sqrt(3); n = 5*sqrt(11); k = 1/sqrt(2);
integrate(function(x) cos(p*log(x)) /(x^n + 1)^k, lower=0, upper=Inf)
(pracma::gammaz((1i*p+1)/n) * pracma::gammaz(k - (1i*p+1)/n) +
	pracma::gammaz((-1i*p+1)/n) * pracma::gammaz(k - (-1i*p+1)/n) ) / (2*n*gamma(k));


### Power 2: I( cos(p * log(x))^2) / ... )

p = sqrt(3); n = 5*sqrt(11); k = 1/sqrt(2);
integrate(function(x) cos(p*log(x))^2 /(x^n + 1)^k, lower=0, upper=Inf)
(pracma::gammaz((2i*p+1)/n) * pracma::gammaz(k - (2i*p+1)/n) +
	+ pracma::gammaz((-2i*p+1)/n) * pracma::gammaz(k - (-2i*p+1)/n) +
	+ 2*gamma(1/n) * gamma(k - 1/n) ) / (4*n*gamma(k));
