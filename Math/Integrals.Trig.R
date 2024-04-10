

### Integrals: Basic Trig


### I( x^(p-1) * sin(k*x) )
# Maths 505: An absolutely wonderful integral!
# https://www.youtube.com/watch?v=xBSxDLQHf3o

# Note:
# - see Mellin Transform;
# - see also Integrals.Trig.Fresnel.R;

p = 1/sqrt(3)
# upper = Inf
integrate(\(x) x^(p-1) * sin(x), 0, 240000, subdivisions = 20024)
# Fresnel integral behaves very badly!
integrate(\(x) sin(x^(1/p)) / p, 0, 1000, subdivisions = 20024)
gamma(p) * sin(pi*p/2)

# upper = Inf
integrate(\(x) x^(p-1) * cos(x), 0, 240000, subdivisions = 20024)
gamma(p) * cos(pi*p/2)

###
p = 1/sqrt(3); k = sqrt(5);
# upper = Inf; drastic increase in subdivisions!
integrate(\(x) x^(p-1) * sin(k*x), 0, 240000, subdivisions = 60024)
gamma(p) * sin(pi*p/2) / k^p

# upper = Inf
integrate(\(x) x^(p-1) * cos(k*x), 0, 240000, subdivisions = 60024)
gamma(p) * cos(pi*p/2) / k^p


### Transforms:
p = 1/sqrt(5); k = sqrt(3);
# upper = Inf; drastic increase in subdivisions!
integrate(\(x) x^(p-1) * sin(k*x) * log(x), 0, 240000, subdivisions = 60024)
gamma(p) * sin(pi*p/2) * (digamma(p) - log(k)) / k^p +
	+ pi/2 * gamma(p) * cos(pi*p/2) / k^p;

