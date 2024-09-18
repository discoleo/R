
### Trig: Formulas

# Various Trig Identities


####################

### Helper Functions

Euler   = 0.57721566490153286060651209008240243079;
Catalan = 0.915965594177219015054603514;


##################
##################

### Trig(3*x)
x = pi/7
cos(3*x) - 4*cos(x)^3 + 3*cos(x)
sin(3*x) + 4*sin(x)^3 - 3*sin(x)
#
tan(3*x) + (4*sin(x)^3 - 3*sin(x)) / (4*cos(x)^3 - 3*cos(x))
tan(3*x) + (tan(x)^3 - 3*tan(x)) / (1 - 3*tan(x)^2)
# x * (x^2 - 3) / (3*x^2 - 1)

### ATAN
x = sqrt(5)
# Note: force result between [-pi/2, pi/2];
atan((3*x - x^3) / (1 - 3*x^2))
3 * atan(x) - pi
#
atan((1 - 3*x^2) / (3*x - x^3))
3/2*pi - 3 * atan(x)

# Transforms:
x = sqrt(5)
atan((2*x^3 - 3*x + 1) / (2*x^3 - 6*x^2 + 3*x))
- 3 * atan(1 - 1/x)
3 * atan(x/(x-1)) - 3*pi/2
#
atan((2*x^3 - 3*x - 1) / (2*x^3 + 6*x^2 + 3*x))
3 * atan(x/(x+1)) - pi/2;

# (1 - x)/(1 + x)
th = pi/11;
x  = (1 - tan(th)) / (1 + tan(th));
(x^3 + 3*x^2 - 3*x - 1) / (x^3 - 3*x^2 - 3*x + 1)
tan(3*th)

# TODO: more;


### 1/3 Angle:
x = pi/7
cos(x/3) # ==
((cos(x) + 1i*sin(x))^(1/3) + (cos(x) - 1i*sin(x))^(1/3)) / 2

### Note:
# D(cos(3*x)) => - 3*sin(3*x), but masked;
x = pi/7
- 12*sin(x)*cos(x)^2 + 3*sin(x)
- 3*sin(3*x)


### Examples:
x = pi/7
integrate(\(x) -(4*sin(x)^3 - 3*sin(x)) / (4*cos(x)^3 - 3*cos(x)), 0, up=x)
integrate(\(x) -(tan(x)^3 - 3*tan(x)) / (1 - 3*tan(x)^2), 0, up=x)
integrate(\(x) -(1 - 4*x^2) / (4*x^3 - 3*x), lower = cos(x), 1)
integrate(\(x) -(x^3 - 3*x) / ((x^2 + 1) * (1 - 3*x^2)), 0, up = tan(x))
integrate(\(x) tan(3*x), 0, up=x)
- log(cos(3*x)) / 3;


### Ex 2:
integrate(\(x) (atan((2*x^3 - 3*x + 1) / (2*x^3 - 6*x^2 + 3*x)) - pi/2) / x, 0, 1/2)
integrate(\(x) (-3*atan(1 - 1/x) - 3*pi/2) / x, 0, 1/2)
integrate(\(x) 3 * (atan(x - 1) - pi/2) / x, 2, Inf)
3/8 * pi*log(2) - 3*Catalan


### Ex 3:
integrate(\(x) atan(- (x^3 - 3*x) / (1 - 3*x^2)), 0, 1)
# strange Bug:
integrate(\(x) atan(tan(3*x)) * (1 + tan(x)^2), 0, pi/6)$value +
	integrate(\(x) atan(tan(3*x)) * (1 + tan(x)^2), pi/6, pi/4)$value
pi*tan(pi/6) - pi + 3*pi/4 - 3/2 * log(2)
pi*tan(pi/6) - pi/4 - 3/2 * log(2)


### Ex 4:
integrate(\(x) atan( (x^3 + 3*x^2 - 3*x - 1) / (x^3 - 3*x^2 - 3*x + 1) ), 0, 1)
integrate(\(x) 6*atan(x) / (x+1)^2, 0, 1)$value +
	+ integrate(\(x) - 2*pi / (x+1)^2, tan(pi/6), 1)$value
integrate(\(x) 6*atan(x) / (x+1)^2 + pi*(1 - 2/(1 + tan(pi/6))), 0, 1)
pi - 2*pi/(1 + tan(pi/6)) + 3*log(2)/2 


#############

#############
### Trig(5*x)
x = pi/7
cos(5*x) - (16*cos(x)^4 - 20*cos(x)^2 + 5)*cos(x)
sin(5*x) - (16*sin(x)^4 - 20*sin(x)^2 + 5)*sin(x)

tan(5*x) - (tan(x)^4 - 10*tan(x)^2 + 5)*tan(x) /
	(1 - 10*tan(x)^2 + 5*tan(x)^4);


### Examples:
x = atan(3) / 5;
tan(x)^5 - 10*tan(x)^3 + 5*tan(x) +
	- tan(5*x) * (5*tan(x)^4 - 10*tan(x)^2 + 1)
# Polynomial:
x = tan(atan(5) / 5 + seq(0,4) * 2*pi/5);
x^5 - 10*x^3 + 5*x - 5 * (5*x^4 - 10*x^2 + 1) # = 0

#
x = tan(atan(6) / 5 + seq(0,4) * 2*pi/5);
x^5 - 10*x^3 + 5*x - 6 * (5*x^4 - 10*x^2 + 1) # = 0


### Complex
x = tan(atan(6i) / 5 + seq(0,4) * 2*pi/5) / 1i;
x^5 + 10*x^3 + 5*x - 6 * (5*x^4 + 10*x^2 + 1) # = 0
round0(x)

#
tn = 7;
x = tan(atan(tn * 1i) / 5 + seq(0,4) * 2*pi/5) / 1i;
x^5 + 10*x^3 + 5*x - tn * (5*x^4 + 10*x^2 + 1) # = 0
round0(x)

# Shift =>
tn = 7;
x = tan(atan(tn * 1i) / 5 + seq(0,4) * 2*pi/5) / 1i - tn;
x^5 - 10*(tn^2 - 1)*x^3 - 20*tn*(tn^2 - 1)*x^2 +
	- 5*(tn^2 - 1)*(3*tn^2 + 1)*x - 4*tn*(tn^4 - 1) # = 0
round0(x)

# Special Case:
x = tan(atan(3i) / 5 + seq(0,4) * 2*pi/5) / 2i - 3/2;
x^5 - 20*x^3 - 60*x^2 - 70*x - 30 # = 0
round0(x)


# TODO: Transforms


### Examples
integrate(\(x) atan((x^5 - 10*x^3 + 5*x) / (5*x^4 - 10*x^2 + 1)), 0, 1)
integrate(\(x) atan(tan(5*x)) * (1 + tan(x)^2), 0, pi/4)
pi*tan(pi/10) + pi/4 - 5/2 * log(2)

# Derivation:
# Note: with Interval splitting;
integrate(\(x) 5*x * (1 + tan(x)^2), 0, pi/4)$value +
	- pi*(tan(pi/4) - tan(pi/10));
integrate(\(x) 5 * atan(x), 0, 1)$value +
	- pi*(tan(pi/4) - tan(pi/10));
- pi*(tan(pi/4) - tan(pi/10)) + 5*pi/4 - 5/2 * log(2)
pi*tan(pi/10) + pi/4 - 5/2 * log(2)


###############

### Gen:
k = pi + 1; # k >= 2!
integrate(\(x) atan(tan(k*x)) * (1 + tan(x)^2), 0, pi/4)
pi*tan(pi/(2*k)) - pi + k*pi/4 - k/2 * log(2)


################
################

### ATAN

###
x = 1/sqrt(7)
atan(1+x) + atan(1-x)
pi/2 - atan(x^2/2)

###
x = 1/sqrt(7)
atan(1 + x*sqrt(2)) + atan(1 - x*sqrt(2))
pi/2 - atan(x^2)

###
x = 1/sqrt(7)
atan(1 + 2*x) + atan(1 - 2*x)
pi/2 - atan(2*x^2)

###
x = 1/sqrt(7); k = sqrt(5);
atan(1 + sqrt(2*k)*x) + atan(1 - sqrt(2*k)*x)
pi/2 - atan(k*x^2)


### Varia:
x = 1/sqrt(7)
atan(x - 1/2) + atan(x - 3/2)
- pi/2 + atan((4*x^2 - 8*x - 1) / (x-1) / 8)

###
x = 1/sqrt(7)
atan(x/(x-1)) + atan(x/(x+1))
- atan(2*x^2)
#
x = sqrt(3)
atan(x/(x-1)) + atan(x/(x+1))
pi - atan(2*x^2)

###
curve(atan(x/(x-1)) + atan(x/(x+1)), -6, 6, ylim = c(-3/4*pi, pi))
curve(- atan(2*x^2), add = TRUE, col = "red", lty=2)
curve(pi - atan(2*x^2), add = TRUE, col = "green", lty=2)
abline(v = c(-1,1), col = "orange", lty=2)


###############
### Complex ###

###
x = sqrt(5)
atan(x + 1i) + atan(x - 1i)
pi/2 + atan(x/2)

###
x = sqrt(5)
atan(sqrt(2) + x * 1i) + atan(sqrt(2) - x * 1i)
pi/2 + atan((x^2 + 1) / sqrt(8))


###
x = sqrt(5)
a = sqrt(2) + 1; b = sqrt(2*a); # fixed;
atan(a + x*b*1i) + atan(a - x*b*1i)
pi/2 + atan(x^2 + 1)
#
atan(a + x*b*1i) - atan(a - x*b*1i)
atan(1i * sqrt(a/2) * (x^2 + sqrt(2)) / x) - pi/2


###
x = sqrt(5)
atan((x + 1i*sqrt(3))/2) + atan((x - 1i*sqrt(3))/2)
pi/2 + atan((x^2 - 1) / (4*x))

