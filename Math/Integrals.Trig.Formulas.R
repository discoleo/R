
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
x = tan(atan(7i) / 5 + seq(0,4) * 2*pi/5) / 1i;
x^5 + 10*x^3 + 5*x - 7 * (5*x^4 + 10*x^2 + 1) # = 0
round0(x)

# TODO: Transforms

