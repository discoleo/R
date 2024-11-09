##################
##
## Integrals: Trig
## Inverse Fractions
##
## Leonard Mada


### Examples
# I( x^p / sin(x) )
# I( x^p / cos(x) )
# - on various intervals;


####################

### Helper Functions

Catalan = 0.915965594177219015054603514;


#################
#################

#################
### "Inverse" ###
#################

library(pracma)

### Refs:
# 1) Maths 505: Integral of x^2/sin(x) from zero to pi/2
#    https://www.youtube.com/watch?v=g4aKTQyETZw
# 2) Maths 505: A cool integral for Apery's constant
#   (ζ(3)): int 0 to 1 (x(1-x))/sin(πx)
#    https://www.youtube.com/watch?v=TNAtG8gXWuU


### I( 1 / sin(x) )
integrate(\(x) 1 / sin(x) - 1/x, 0, pi/2)
- log(pi/4)

###
integrate(\(x) x / sin(x), 0, pi/2)
2*Catalan

###
integrate(\(x) x^2 / sin(x), 0, pi/2)
2*pi*Catalan - 7/2*pracma::zeta(3)

###
integrate(\(x) x^3 / sin(x), 0, pi/2)
(pracma::psi(3, 3/4) - pracma::psi(3, 1/4)) / 128 + 3/2*pi^2*Catalan

###
integrate(\(x) x^4 / sin(x), 0, pi/2)
3/2*31*zeta(5) + pi^3*Catalan +
	+ pi*(pracma::psi(3, 3/4) - pracma::psi(3, 1/4)) / 64;


### cos(x)

# Note:
# - can be used to derive I( x^p * tan(x) )

### I( 1 / cos(x) )
integrate(\(x) 1/cos(x) - 1/(pi/2-x), 0, pi/2)
- log(pi/4)

### x^1
integrate(\(x) (pi/2 - x) / cos(x), 0, pi/2)
2 * Catalan

### x^2
integrate(\(x) (pi^2/4 - x^2) / cos(x), 0, pi/2)
7/2*pracma::zeta(3)

### x^3
integrate(\(x) x^2 * (pi/2 - x) / cos(x), 0, pi/2)
(pracma::psi(3, 3/4) - pracma::psi(3, 1/4)) / 128 +
	+ 7/2*pi*pracma::zeta(3);

### x^4
integrate(\(x) x^2 * (pi^2/4 - x^2) / cos(x), 0, pi/2)
- 3*31/2*pracma::zeta(5) + 35/8 * pi^2 * pracma::zeta(3);

# Variant:
integrate(\(x) x^3 * (pi/2 - x) / cos(x), 0, pi/2)
- 3*31/2 * pracma::zeta(5) + 21/8 * pi^2 * pracma::zeta(3) +
- (pracma::psi(3, 3/4) - pracma::psi(3, 1/4)) * pi / 256;


# Alternative Forms:
integrate(\(x) x / cos(x) - (pi/2) / (pi/2 - x), 0, pi/2)
- pi/2 * log(pi/4) - 2*Catalan;

integrate(\(x) x^2 / cos(x) - (pi/2)^2 / (pi/2 - x), 0, pi/2)
- 7/2*pracma::zeta(3) - (pi/2)^2 * log(pi/4);

integrate(\(x) x^3 / cos(x) - (pi/2)^3 / (pi/2 - x), 0, pi/2)
- (pracma::psi(3, 3/4) - pracma::psi(3, 1/4)) / 128 +
	- 21/4*pi*pracma::zeta(3) - (pi/2)^3 * log(pi/4);


### Note: Psi-Relations

7/4 * zeta(3) # ==
- (pracma::psi(2, 3/4) + pracma::psi(2, 1/4)) / 64;
# while:
pi^3 / 16 # ==
(pracma::psi(2, 3/4) - pracma::psi(2, 1/4)) / 64;

#
13/16 * zeta(3) # ==
- (pracma::psi(2, 2/3) + pracma::psi(2, 1/3)) / 64;
91/16 * zeta(3) # ==
- (pracma::psi(2, 5/6) + pracma::psi(2, 1/6)) / 64;

###
45/4 * zeta(4) # == 2 * pi^4 ==
(pracma::psi(3, 3/4) + pracma::psi(3, 1/4)) / 128;

###
93 * zeta(5) # ==
- (pracma::psi(4, 3/4) + pracma::psi(4, 1/4)) / 256;
# while:
5/16 * pi^5 # ==
(pracma::psi(4, 3/4) - pracma::psi(4, 1/4)) / 256;


### Composed: on [0, pi]

# - some terms cancel out when using: x*(pi - x);
# - for intuition of underlying principle, see e.g:
#   blackpenredpen: I saw this NTU entrance exam integral on Dcard
#   and I just had to integrate it
#   https://www.youtube.com/watch?v=PWGJGnji5Nk

### I( 1 / sin(x) )
integrate(\(x) 1 / sin(x) - 1/x - 1/(pi-x), 0, pi)
- 2*log(pi/2)

### I( x / sin(x) )
integrate(\(x) x / sin(x) - pi / (pi-x), 0, pi)
- pi*log(pi/2)

### I( x^2 / sin(x) )
integrate(\(x) x^2 / sin(x) - pi^2 / (pi-x), 0, pi)
- 7*pracma::zeta(3) - pi^2*log(pi/2)

###
integrate(\(x) x*(pi - x) / sin(x), 0, pi)
7*zeta(3)

###
pracma::integral(\(x) x*(pi - x)*(pi + x) / sin(x), -pi, pi)
3*7*pi*zeta(3)

###
pracma::integral(\(x) x^2*(pi - x)*(pi + x) / sin(x), 0, pi)
2*7*pi^2*zeta(3) - 3*31*zeta(5)

###
integral(\(x) x^3*(pi - x)*(pi + x) / sin(x), -pi, pi)
7*7*pi^3*zeta(3) - 15*31*pi*zeta(5)

###
integral(\(x) x^4*(pi - x)*(pi + x) / sin(x), 0, pi)
9/2*7*pi^4*zeta(3) - 39/2*31*pi^2*zeta(5) + 45/2*127*zeta(7)

###
integral(\(x) x^5*(pi - x)*(pi + x) / sin(x), -pi, pi)
11*7*pi^5*zeta(3) - 90*31*pi^3*zeta(5) + 315/2*127*pi*zeta(7)


### Other Intervals

### on [0, pi/3]

### I( 1 / sin(x) )
integrate(\(x) 1 / sin(x) - 1/x, 0, pi/3)
log(3)/2 - log(pi/2)

### I( x / sin(x) )
integrate(\(x) x / sin(x), 0, pi/3)
pi/3 * log(tan(pi/6)) + sin(pi/3) / 12^2 *
	( 3 * pracma::psi(1, 1/6) - 3 * pracma::psi(1, 5/6) +
	+ 5 * pracma::psi(1, 1/3) - 5 * pracma::psi(1, 2/3) +
	+ pracma::psi(1, 1/12) - pracma::psi(1, 11/12) +
	- pracma::psi(1, 5/12) + pracma::psi(1,  7/12) );


# Derivation:

# on [0, pi/3]
integrate(\(x) log(tan(x)), 0, pi/6)
- sin(pi/3)/(2*12^2) *
	( 3 * pracma::psi(1, 1/6) - 3 * pracma::psi(1, 5/6) +
	+ 5 * pracma::psi(1, 1/3) - 5 * pracma::psi(1, 2/3) +
	+ pracma::psi(1, 1/12) - pracma::psi(1, 11/12) +
	- pracma::psi(1, 5/12) + pracma::psi(1,  7/12) );
# by parts =>
integrate(\(x) x / (cos(x)*sin(x)), 0, pi/6)
integrate(\(x) 1/2 * x / sin(x), 0, pi/3)
pi/6 * log(tan(pi/6)) + sin(pi/3)/(2*12^2) *
	( 3 * pracma::psi(1, 1/6) - 3 * pracma::psi(1, 5/6) +
	+ 5 * pracma::psi(1, 1/3) - 5 * pracma::psi(1, 2/3) +
	+ pracma::psi(1, 1/12) - pracma::psi(1, 11/12) +
	- pracma::psi(1, 5/12) + pracma::psi(1,  7/12) );


##############

###########
### TAN ###
###########

###
integrate(\(x) 1/tan(x) - 1/x, 0, pi/2)
- log(pi/2)

###
integrate(\(x) 1/tan(x) - 1/x, 0, pi/4)
- log(2)/2 - log(pi/4)

###
k = pi/5
integrate(\(x) 1/tan(x) - 1/x, 0, k)
log(sin(k)) - log(k)


###
integrate(\(x) x / tan(x), 0, pi/2)
pi*log(2)/2
#
integrate(\(x) 2*x * (1/tan(x) - tan(x)), 0, pi/4)

# =>
integrate(\(x) x / tan(x), 0, pi/4)
pi*log(2)/8 + Catalan/2


### on [0, pi/3]
# - see file: Integrals.Log.Trig.R;
integrate(\(x) x / tan(x), 0, pi/3)
pi*log(3)/6 - sqrt(3)/(8*6^2) *
	( pracma::psi(1, 1/12) - pracma::psi(1, 11/12) +
	- pracma::psi(1, 5/12) + pracma::psi(1, 7/12) +
	- 5 * pracma::psi(1, 1/6) + 5 * pracma::psi(1, 5/6) +
	- 3 * pracma::psi(1, 2/6) + 3 * pracma::psi(1, 4/6));

###
integrate(\(x) cos(x) / sin(x)^2 - 1/x^2, 0, pi/3)
-2/sqrt(3) + 3/pi

### [solvable on any interval]
integrate(\(x) x * cos(x) / sin(x)^2 - 1/x, 0, pi/3)
log(3)/2 - log(pi/2) - pi/3 / sin(pi/3) + 1


###
integrate(function(x) x*(pi-x)*(pi+x) / tan(x), 0, pi)
- 3/2 * pi * zeta(3)

###
integrate(function(x) x*(pi-x)*(pi+x) / tan(x), 0, pi/2)
9/16 * pi * zeta(3) + 3/8 * pi^3 * log(2)

###
integrate(function(x) x^2 / tan(x), 0, pi/2)
- 7/8 * zeta(3) + 1/4 * pi^2 * log(2)

###
integrate(function(x) x^3 / tan(x), 0, pi/2)
- 9/16 * pi * zeta(3) + 1/8 * pi^3 * log(2)

###
integrate(function(x) x^4 / tan(x), 0, pi/2)
93/32 * zeta(5) - 9/16 * pi^2 * zeta(3) + 1/16 * pi^4 * log(2)


### Helper:
k = 1/5;
fk = function(x, k) log((1 - cos(k*x)) / (1 + cos(k*x))) / (2*k);
integrate(\(x) 1 / sin(k*x), pi/5, pi/3)
fk(pi/3, k) - fk(pi/5, k)


### Derived
# on [0, pi]

###
integrate(function(x) x*(pi - x) / sin(x), 0, pi)
7*zeta(3)

###
integrate(function(x) x^2*(pi - x) / sin(x), 0, pi)
1/2*7*pi*zeta(3)

###
integral(function(x) x^3*(pi - x) / sin(x), 0, pi)
3/2*7*pi^2*zeta(3) - 3*31*zeta(5)

###
integral(function(x) x^4*(pi - x) / sin(x), 0, pi)
4/2*7*pi^3*zeta(3) - 9/2*31*pi*zeta(5)

###
integral(function(x) x^5*(pi - x) / sin(x), 0, pi)
5/2*7*pi^4*zeta(3) - 30/2*31*pi^2*zeta(5) + 45/2*127*zeta(7)


### Derived
# on [0, pi/2]

###
integrate(function(x) x*(pi - x) / sin(x), 0, pi/2)
7/2*zeta(3)

###
integrate(function(x) x^2 / sin(x), 0, pi/2)
2*pi*Catalan - 7/2*zeta(3)

###
integrate(function(x) x^2*(pi - x) / sin(x), 0, pi/2)
1/2*pi^2*Catalan - 7/2*pi*zeta(3) +
	- (pracma::psi(3, 3/4) - pracma::psi(3, 1/4)) / 128;

###
integrate(function(x) x^3*(pi - x) / sin(x), 0, pi/2)
3/2*7*pi^2*zeta(3) - 3*31*zeta(5) - integrate(function(x) x*(pi - x)^3 / sin(x), 0, pi/2)$value
- 3/2*31*zeta(5) + 1/2*pi^3*Catalan +
	- pi*(pracma::psi(3, 3/4) - pracma::psi(3, 1/4)) / 128;

###
integrate(function(x) x^2*(pi - x)^2 / sin(x), 0, pi/2)
- 1/2*7*pi^2*zeta(3) + 3/2*31*zeta(5)


#################
#################

### Higher Powers

### I( x^p / sin(x)^2 )
# Maths 505: A surprisingly fascinating integral
# https://www.youtube.com/watch?v=8DWa0zIM9lY
# - see also file: Integrals.Log.Trig.R;

###
integrate(\(x) 1 / sin(x)^3 - 1/x^3 - 1/(2*x), 0, pi/2)
1/(2*pi^2/4) + log(2)/2 - log(pi/2)/2 - 1/12
2/pi^2 + log(2) - log(pi)/2 - 1/12
# Note: for arbitrary interval, see a few paragraphs below;

### I( 1 / sin(x)^2 )
integrate(\(x) 1 / sin(x)^2 - 1/x^2, 0, pi/2)
2/pi
#
integrate(\(x) 1 / sin(x)^2 - 1/x^2, 0, pi/4)
4/pi - 1

### I( x / sin(x)^2 )
integrate(\(x) x / sin(x)^2 - 1/x, 0, pi/2)
1 - log(pi/2)
#
integrate(\(x) x / sin(x)^2 - 1/x, 0, pi/4)
3/2*log(2) - pi/4 + 1 - log(pi)


### I( x^2 / sin(x)^2 )
integrate(\(x) x^2 / sin(x)^2, 0, pi/2)
pi*log(2)


### I( x^3 / sin(x)^2 )
integrate(\(x) x^3 / sin(x)^2, 0, pi/2)
3/4*pi^2*log(2) - 3*7/8 * pracma::zeta(3)


### I( x^4 / sin(x)^2 )
integrate(\(x) x^4 / sin(x)^2, 0, pi/2)
1/2*pi^3*log(2) - 9/4 * pi * pracma::zeta(3)


### I( 1 / sin(x)^3 )
integrate(\(x) 1 / sin(x)^3 - 1/x^3 - 1/(2*x), 0, pi/2)
log(2)/2 - 1/12 + (2/pi^2 - log(pi/2)/2);

# Note:
# I( 1 / sin(x)^(2*n+1) ) can be computed on arbitrary intervals;

### I( 1 / sin(x)^3 )
up = pi/5;
integrate(\(x) 1 / sin(x)^3 - 1/x^3 - 1/(2*x), 0, up)
(log( (1 - cos(up)) / (1 + cos(up)) ) +
	+ 1/(cos(up) - 1) + 1/(cos(up) + 1) + 2*log(2) - 1/3) / 4 +
	+ (1/(2*up^2) - log(up)/2);

# Lim: x -> 0
x = 1E-3; # but numerical issues!
(-log(1 - cos(x)) - 1/(cos(x)-1)) / 2 - 1/(x^2) + log(x)
log(2)/2 + 1/12


### on [0, pi/4]

### I( 1 / sin(x)^2 )
integrate(function(x) 1 / sin(x)^2 - 1/x^2, 0, pi/4)
4/pi - 1
#
integrate(function(x) 1 / sin(x)^2, pi/4, pi/2)
1;

### I( x / sin(x)^2 )
integrate(function(x) x / sin(x)^2 - 1/x, 0, pi/4)
1 - pi/4 - 1/2 * log(2) - log(pi/4)
# Note: log(pi/4) is from 1/x;
integrate(function(x) x / sin(x)^2, pi/4, pi/2)
pi/4 + log(2)/2

### I( x^2 / sin(x)^2 ) on [0, pi/4]
integrate(function(x) x^2 / sin(x)^2, 0, pi/4)
- pi^2 / 16 + pi/4 * log(2) + Catalan


### Generalisation
### Arbitrary Interval:

### Base: I( 1 / sin(k*x) )
k = pi/7;
integrate(\(x) 1 / sin(k*x) - 1/(k*x), 0, 1)
log((1-cos(k))/(1+cos(k))) / (2*k) + log(2)/k - log(k)/k;

### Base: I( 1 / tan(k*x) )
k = pi/7
integrate(\(x) 1 / tan(k*x) - 1/(k*x), 0, 1)
log(sin(k)) / k - log(k) / k;

### Base: I( 1 / sin(k*x)^3 )
k = pi/7
integrate(\(x) 1 / sin(k*x)^3 - 1/(k^3*x^3) - 1/(2*k*x), 0, 1)
(log( (1 - cos(k)) / (1 + cos(k)) ) +
	+ 1/(cos(k) - 1) + 1/(cos(k) + 1) + 2*log(2) - 1/3) / (4*k) +
	+ 1/(2*k^3) - log(k)/(2*k);

### Base: I( cos(k*x) / sin(k*x)^3 )
k = pi/7
integrate(\(x) cos(k*x) / sin(k*x)^3 - 1/(k^3*x^3), 0, 1)
- sin(k)^(-2) / (2*k) + 1/(2*k^3) + 1/(6*k);

### I( 1 / tan(k*x)^3 )
# Diff =>
k = 5*pi/7
integrate(\(x) 1 / tan(k*x)^3 - 1/(k^3*x^3) + 1/(k*x), 0, 1)
- sin(k)^(-2) / (2*k) + 1/(2*k^3) + 1/(6*k) +
	- log(sin(k)) / k + log(k) / k;


### I( x / ... )
# D[k] =>

### I( x / sin(k*x)^2 )
k = pi/7
integrate(\(x) x / sin(k*x)^2 - 1/(k^2*x), 0, 1)
- cos(k) / (k * sin(k)) + log(sin(k)) / k^2 - log(k) / k^2 + 1/k^2;

### I( x / tan(k*x)^2 )
k = pi/7
integrate(\(x) x / tan(k*x)^2 - 1/(k^2*x), 0, 1)
(log(sin(k)) - log(k) - k/tan(k) + 1) / k^2 - 1/2


### I( x * cos(k*x) / sin(k*x)^2 )
k = 3*pi/7;
integrate(\(x) x * cos(k*x) / sin(k*x)^2 - 1/(k^2*x), 0, 1)
log((1-cos(k))/(1+cos(k))) / (2*k^2) - 1 / (k*sin(k)) +
	+ log(2)/k^2 - log(k)/k^2 + 1/k^2;

### I( x * cos(k*x) / sin(k*x)^4 )
k = pi/2 + 1/11;
integrate(\(x) x * cos(k*x) / sin(k*x)^4 - 1/(k^4*x^3) - 1/(6*k^2*x), 0, 1)
- (2/sin(k) + 2*(cos(k)^2 + 1)/sin(k)^3) / (12*k) +
	+ (log( (1 - cos(k)) / (1 + cos(k)) ) +
		- 2*cos(k)/sin(k)^2 + 2*log(2) - 1/3) / (12*k^2) +
	+ 1/(2*k^4) + 1/(6*k^2) - log(k)/(6*k^2);


### I( x / sin(k*x)^4 )
k = pi/7
integrate(\(x) x / sin(k*x)^4 - 1/(k^4*x^3) - 2/(3*k^2*x), 0, 1)
- cos(k) / (3*k * sin(k)^3) - 1 / (6*k^2 * sin(k)^2) - 2*cos(k) / (3*k * sin(k)) +
	+ (12*log(sin(k)) - 12*log(k) + 9/k^2 + 13) / (18*k^2);


### Varia:

# By Parts:
integrate(\(x) x^2 / tan(x)^2, 0, pi/2)
pi*log(2) - pi^3/24

# By Parts: => 1/2 * I( x^2 * (1 + 1/tan(x)^2) )
integrate(\(x) x / tan(x), 0, pi/2)
pi*log(2)/2

# By Parts:
integrate(\(x) log(sin(x)), 0, pi/2)
- pi*log(2)/2

# - for other intervals, see file:
#   Integrals.Log.Trig.R;

