

### Various Radicals

### I( sqrt(b^2 + sqrt(b^4 + x^4)) ) on [0, 1]
# Maths 505: A deceivingly tough integral
# https://www.youtube.com/watch?v=nGx8j7-KHxI
# Substitution: x^2 = tan(u);

integrate(\(x) sqrt(1 + sqrt(1 + x^4)), 0, 1)
v = sqrt(tan(pi/8));
log((1+v)/(1-v))/4 + v/(1 - v^4) + atan(v)/2;

###
lim = sqrt(5)
integrate(\(x) sqrt(1 + sqrt(1 + x^4)), 0, lim)
v = sqrt(tan(atan(lim^2)/2));
log((1+v)/(1-v))/4 + v/(1 - v^4) + atan(v)/2;

### Gen:
b = 5^(1/4)
integrate(\(x) sqrt(b^2 + sqrt(b^4 + x^4)), 0, 1)
v = sqrt(tan(atan(1/b^2)/2));
b^2 * (log((1+v)/(1-v))/4 + v/(1 - v^4) + atan(v)/2);

# TODO: D(b)

######################
######################

### I( 1 / sqrt(x * (x^2 + b^2) * (x + sqrt(x^2 + b^2))) )
# Maths 505: A crazy nested roots integral!
# https://www.youtube.com/watch?v=vumE5v6p-mg
# Classic Substitution: x = tan(u)

lim = sqrt(5)
integrate(\(x) 1 / sqrt(x * (x^2 + 1) * (x + sqrt(x^2 + 1))), 0, lim)
2*sqrt(2) * atan(sqrt(tan(atan(lim)/2)))
2*sqrt(2) * atan(sqrt((sqrt(lim^2 + 1) - 1)/lim))

###
b = sqrt(5)
integrate(\(x) 1 / sqrt(x * (x^2 + b^2) * (x + sqrt(x^2 + b^2))), 0, 1)
2*sqrt(2) * atan(sqrt(sqrt(b^2 + 1) - b)) / b


### Diff-Type
lim = sqrt(5)
integrate(\(x) 1 / sqrt(x * (x^2 + 1) * (sqrt(x^2 + 1) - x)), 0, lim)
v = sqrt((sqrt(lim^2 + 1) - 1)/lim);
sqrt(2) * log((1 + v) / (1 - v))

###
b = sqrt(5)
integrate(\(x) 1 / sqrt(x * (x^2 + b^2) * (sqrt(x^2 + b^2) - x)), 0, 1)
v = sqrt((sqrt(b^2 + 1) - b));
sqrt(2) * log((1 + v) / (1 - v)) / b;

# Derivation:
integrate(\(x) 1 / sqrt(sin(x) * (1 - sin(x))), 0, atan(lim))
v = sqrt(tan(atan(lim)/2));
sqrt(2) * log((1 + v) / (1 - v))

