

### Trig & Exp:
### Basic Functions

### Examples:

# I( x^p * TRIG(k*x) * exp(- a*x) )


#####################

### I( x^p * sin(k*x) * exp(- a*x) )

### Basic:
a = sqrt(3); k = sqrt(5);
integrate(\(x) sin(k*x) * exp(- a*x), 0, Inf)
k / (a^2 + k^2)

### Ext: x^p * ...
p = sqrt(3)
integrate(\(x) x^p * sin(x) / exp(x), 0, Inf)
gamma(p + 1) * Im(1 / (1-1i)^(p + 1))

### Gen 1:
p = sqrt(3); a = sqrt(5);
integrate(\(x) x^p * sin(x) / exp(a*x), 0, Inf)
gamma(p + 1) * Im(1 / (a - 1i)^(p + 1))

### Gen:
p = sqrt(3); a = sqrt(5) - sqrt(2); k = sqrt(3);
integrate(\(x) x^p * sin(k*x) / exp(a*x), 0, Inf)
gamma(p + 1) * Im(1 / (a - k*1i)^(p + 1))


#####################

### I( x^p * cos(k*x) * exp(- a*x) )

###
a = sqrt(3); k = 1/sqrt(5);
integrate(\(x) cos(k*x) * exp(- a*x), 0, Inf)
a / (a^2 + k^2)

###
p = sqrt(3)
integrate(\(x) x^p * cos(x) / exp(x), 0, Inf)
Re(gamma(p + 1) / (1-1i)^(p + 1))

###
p = sqrt(3); a = sqrt(5);
integrate(\(x) x^p * cos(x) / exp(a*x), 0, Inf)
gamma(p + 1) * Re(1 / (a - 1i)^(p + 1))

###
p = sqrt(3); a = sqrt(5) - sqrt(2); k = sqrt(3);
integrate(\(x) x^p * cos(k*x) / exp(a*x), 0, Inf)
gamma(p + 1) * Re(1 / (a - k*1i)^(p + 1))

### Case: p = -1
a = sqrt(5) - sqrt(2); k = sqrt(3);
integrate(\(x) (exp(-x) - cos(k*x) * exp(-a*x)) / x, 0, Inf)
log(a^2 + k^2) / 2


##########
### Simple

### I( sin(x)^n / x^n * exp(-x) )
# Maths 505: Feynman's technique is incredibly overpowered!
# https://www.youtube.com/watch?v=_d1PJiIm-lE

integrate(\(x) sin(x)^2/x^2 / exp(x), 0, Inf)
atan(2) - log(5)/4

###
integrate(\(x) sin(x)^3/x^3 / exp(x), 0, Inf)
atan(3) - 3*log(5)/8


# Derivation:
t = sqrt(2)
integrate(\(x) sin(x)^3 / exp(t*x), 0, Inf)
integrate(\(x) (3*sin(x) - sin(3*x))/4 / exp(t*x), 0, Inf)
Im(1/(t + 3i) - 3/(t + 1i)) / 4
# =>
3/4*(1/(t^2 + 1) - 1/(t^2 + 9))

#
integrate(\(x) sin(x)^3/x / exp(x), 0, Inf, rel.tol=1E-10)
# Constant of Integration = pi/4;
pi/4 - (atan(1) - 1/3*atan(1/3))*3/4
(atan(1/3) + pi/4)/4


###
t = sqrt(2)
integrate(\(x) sin(x)^5 / exp(t*x), 0, Inf)
integrate(\(x) (sin(5*x) + 20*sin(x)^3 - 5*sin(x))/16 / exp(t*x), 0, Inf)
integrate(\(x) (sin(5*x) - 5*sin(3*x) + 10*sin(x))/16 / exp(t*x), 0, Inf)
- Im(1/(t + 5i) - 5/(t + 3i) + 10/(t + 1i)) / 16
# =>
5/16*(1/(t^2 + 25) - 3/(t^2 + 9) + 2/(t^2 + 1))


###
t = sqrt(2)
integrate(\(x) sin(x)^7 / exp(t*x), 0, Inf)
Im(1/(t + 7i) - 7/(t + 5i) + 21/(t + 3i) - 35/(t + 1i)) / 4^3
# =>
-7/4^3 * (1/(t^2 + 49) - 5/(t^2 + 25) + 9/(t^2 + 9) - 5/(t^2 + 1))


###
t = sqrt(2)
integrate(\(x) sin(x)^9 / exp(t*x), 0, Inf)
- Im(1/(t + 9i) - 9/(t + 7i) + 36/(t + 5i) - 84/(t + 3i) + 126/(t + 1i)) / 4^4
9/4^4 * (1/(t^2 + 81) - 7/(t^2 + 49) + 20/(t^2 + 25) - 28/(t^2 + 9) + 14/(t^2 + 1))


### I(I(...)): I[n]( atan(...) )
# x/a * atan(x/a) - 1/2*log(x^2 + a^2)
a = sqrt(3)
integrate(\(x) 1/a * atan(x/a), 0, 1)
1/a * atan(1/a) - log(a^2 + 1)/2 + log(a)

### I[2]
a = sqrt(3)
# x^2/(2*a) * atan(x/a) - a/2 * atan(x/a) - 1/2*x*log(x^2 + a^2) + x/2;
integrate(\(x) x/a * atan(x/a) - 1/2*log(x^2 + a^2), 0, 1)
1/(2*a) * atan(1/a) - a/2 * atan(1/a) - 1/2*log(a^2 + 1) + 1/2

