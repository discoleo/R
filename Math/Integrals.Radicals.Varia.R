

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

