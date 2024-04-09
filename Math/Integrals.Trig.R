

### Integrals: Basic Trig


### I( x^(p-1) * sin(k*x) )
# Maths 505: An absolutely wonderful integral!
# https://www.youtube.com/watch?v=xBSxDLQHf3o

# Note: see Mellin Transform

p = 1/sqrt(3)
# upper = Inf
integrate(\(x) x^(p-1) * sin(x), 0, 240000, subdivisions = 20024)
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

