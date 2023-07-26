


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
# s in (-1, 1)
# see file: Integrals.Trig.Fractions.Poly.R;
p = 2 - sqrt(2)
k = sqrt(3)
integrate(\(x) sin(k*x) / x^(p+1), 0, 200000, subdivisions=40000)
- gamma(-p) * sin(pi*p/2) * k^p;

