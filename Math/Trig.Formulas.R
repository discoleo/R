#################
##
## Trig: Formulas
## ATAN & ACOS


### Resources

# - Other Trig formulas are also found
#   in file: Integrals.Fractions.Unity.Definite.R;
# TODO: copy/move here?

#################

### PROD( TRIG )
# 1. Michael Penn: do you know the value of this product??
#    https://www.youtube.com/watch?v=k0bvmRSEIvQ

### Prod( COS )
n = 8
id = seq(n);
cs = cos(pi * id / (2*n+1));
#
prod(cs)
2^(-n)


### Prod( SIN )
n = 7 # ODD
nh = (n - 1)/2; id = seq(nh);
sn = sin(pi * id / n);
#
prod(sn)
sqrt(n) / 2^nh;


### SUM

### Sum( n * SIN() )
n = 11 # ODD
p = 3
id = seq((n-1)/2); sgn = (-1)^((p %% 2) + 1);
#
sum(id * sin(2*pi*p*id/n)) # ==
sgn * n/4 / sin(pi*p/n)


### Sum( cos * log(sin) )
# - part of digamma;

n = 7; p = 3; # n = ODD
id = seq((n-1) / 2);
cs = cos(2*pi*p*id/n); sn = sin(pi*id/n)
#
sum(cs * log(sn))
(digamma(p/n) + Euler + log(2*n)) / 2 + pi/4 / tan(p*pi/n)


### Sum( cos * log(cos) )

n = 7; p = 3; # n = ODD
id = seq((n-1) / 2);
cs = cos(2*pi*p*id/n); cs2 = cos(pi*id/n)
#
sum(cs * log(cs2))
(digamma(p/(2*n) + 1/2) - digamma(p/(2*n))) / 4 - pi/4 / sin(pi*p/n)


### Other

###
n = 5
id = seq(2, n-1, by = 2)
x = 2^(1/n); cs = cos(id*pi/n); sn = sin(id*pi/n);
cs2 = cos(id*pi/2/n); sn2 = sin(id*pi/2/n);

sum((1 - cs) * log(2 - 2*cs))
n*log(2) + 2*log(sqrt(n) / 2^((n-1)/2)) - 2*sum(cs * log(sn2))
log(2) + log(n) - 2*sum(cs * log(sn2))
- (digamma(1/n) + Euler + pi/2 / tan(pi/n))


#################

###
Im(acos(3/2 + 0i))
Re(log((sqrt(5) + 3) * 1i/2))
2 * log((1+sqrt(5))/2)


###
Im(acos(4/3 + 0i))
log((sqrt(7) + 4)/3)

###
Im(acos(5/3 + 0i))
log(3)
# varia:
(2*digamma(1/2) - digamma(1/2 - 1/3) - digamma(1/2 + 1/3)) / 3


###
Im(acos(2 + 0i))
log(sqrt(3) + 2)


###
Im(acos(3 + 0i))
2 * log(sqrt(2) + 1)

###
Im(acos(4 + 0i))
log((sqrt(15) + 4))

###
Im(acos(sqrt(2) + 0i))
log(sqrt(2) + 1)


### Gen:
n = sqrt(7)
Im(acos(n + 0i))
log((sqrt(n^2 - 1) + n))


##################
##################

###
n = 5
# n = 7 # ODD
#
id = seq(2, n-1, by = 2)
x = 2^(1/n); cs = cos(id*pi/n); sn = sin(id*pi/n);

sum(2 * sn * atan((1 - cs) / sn))
pi/2 / sin(pi/n)

