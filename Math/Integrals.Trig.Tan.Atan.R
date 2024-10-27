##################
##
## Leonard Mada
##
## Integrals: Trig
## Variants: ATAN


####################

### Helper Constants
Euler   = 0.577215664901532860606512090;
Catalan = 0.915965594177219015054603514;

#####################
#####################

### Basic Integrals

### I( atan(x) / x )
integrate(\(x) atan(x) / x, 0, 1)
Catalan

# Varia:
x = exp(pracma::lambertWp(exp(-1)) / 2 + 1/2)
# Maximum of function:
log(x) / (x^2 + 1)


### I( atan(x) / x^2 )
integrate(\(x) - (atan(x) / x^2 - 1/x), 0, 1)
pi/4 + log(2)/2 - 1


####################
####################

### Log-Combinations

### I( atan(x) * log((1-x)/(1+x)) )
# Maths 505: Another INSANE integral!
# https://www.youtube.com/watch?v=KEDEzVqlAYU

integrate(\(x) atan(x) * log((1-x)/(1+x)), 0, 1)
pi^2/16 - pi/4*log(2) - Catalan;


### I( atan(x) * log(1-x) )
integrate(\(x) atan(x) * log(1-x), 0, 1)
integrate(\(x) ((1-x)*log(1-x) + x) / (x^2+1) - pi/4, 0, 1)
5/96*pi^2 + pi*log(2)/8 - log(2)^2 / 8 + log(2)/2 - pi/4 - Catalan;

### I( atan(x) * log(1+x) )
integrate(\(x) atan(x) * log(1+x), 0, 1)
-1/96*pi^2 + 3/8*pi*log(2) - log(2)^2 / 8 + log(2)/2 - pi/4;


###
integrate(\(x) atan(x) * atan(1-x), 0, 1)
integrate(\(x) (1-2*x) * atan(1-x) / (x^2+1), 0, 1)
# TODO

#
integrate(\(x) atan(1-x) / (x^2+1), 0, 1)
1.393582/2 + # (Li2((3-1i)/5) + Li2((3+1i)/5)) / 2 +
	+ (log(5)*log(5/4)/4 - pi^2/8  + pi/4*atan(1/2) + atan(1/2)*atan(1/3))/2;

#
integrate(\(x) x * atan(1-x) / (x^2+1), 0, 1)
-0.59849679 * 1i/2i + # (Li2((3-1i)/5) - Li2((3+1i)/5)) / 2i +
	+ (2*log(2)*atan(1/2) + log(5)*(pi/4 + atan(1/3) - atan(1/2))) / 4;


# Gen 1: TODO
k = 2
integrate(function(x) atan(x) * log((k - x)/(k + x)), 0, 1)
(pi/4 - log(2)/2)*log((k-1)/(k+1)) +
	+ integrate(function(x) (x*atan(x) - log(x^2+1)/2)*(1/(k - x) + 1/(k + x)), 0, 1)$value
(pi/4 - log(2)/2)*log((k-1)/(k+1)) +
	+ integrate(function(x) k*atan(x)*(1/(k-x) - 1/(x+k)), 0, 1)$value +
	- integrate(function(x) log(x^2+1)/2*(1/(k - x) + 1/(k + x)), 0, 1)$value


####################
####################

### I( log(x) * atan(x) / x )
# 1.) Maths 505: a superb integral sprinkled with some fourier analysis
# https://www.youtube.com/watch?v=aHCLS4l65VE
# 2.a) see also: Sums.Fractions.Higher.R;
# 2.b) Flammable Maths: An AMAZING Journey of Series Evaluation!
# Calculating Euler's Sum! [ Series pi^3/32 (-1)^k/(2k+1)^3 ]

integrate(\(x) log(x) * atan(x) / x, 0, 1)
- pi^3 / 32


### I( log(x) * atan(x) )
integrate(\(x) log(x) * atan(x), 0, 1)
pi^2/48 - pi/4 + log(2)/2

# Derivation:
id = seq(0, 20000)
- sum((-1)^id * (1/(2*id+1) - 1/(2*id+2))) + pi^2/48
- sum((-1)^id / (2*id+1)) + pi^2/48 + log(2)/2
pi^2/48 - pi/4 + log(2)/2


### I( log(x) * atan(x^p) / x )
integrate(\(x) log(x) * atan(x^2) / x, 0, 1)
- pi^3 / 128

### Gen:
p = sqrt(5)
integrate(\(x) log(x) * atan(x^p) / x, 0, 1)
- pi^3 / (32 * p^2)

### Other:
integrate(\(x) log(x) * x / (x^2 + 1), 0, 1)
- pi^2/48


### Pow = 2
p = sqrt(5)
integrate(\(x) log(x)^2 * x^(p-1) / (x^(2*p) + 1), 0, 1)
pi^3 / (16 * p^3)

###
integrate(\(x) log(x)^2 / (x^2 + 1), 0, 1)
pi^3 / 16


###
integrate(\(x) log(x^2 + 1) / x, 0, 1)
pi^2/24


###############
###############

### Experiments

FUN = function(x, normalize = FALSE) {
	y = sapply(x, \(lim) {
			integrate(\(x) atan(x) * atan(lim - x), 0, lim)$value;
	});
	if(normalize) y = y / x;
	return(y);
}

#
up = 3; # up = 10; # up = 50;
curve(FUN(x), 0, up)
curve(exp(x) - 1, add=TRUE, col="red")
curve(tan(pi/2 * x/up), add=TRUE, col="green")
curve(eval(x), add=TRUE, col="purple")
curve(eval(pi/2*x), add=TRUE, col="pink")

#
up = 10; # up = 30; # up = 500;
curve(FUN(x, normalize = T), 0, up)
curve(atan(x), add=TRUE, col="pink")

# Limit: lim -> Inf
# library(Rmpfr)
lim = mpfr("1E+8", 240)
integrate(\(x) {
	x = mpfr(x, 240);
	y = atan(x) * atan(lim - x) / lim;
	as.numeric(y);
}, 0, as.numeric(lim))
pi^2 / 4

# Note: the proof is left as an exercise for the reader;

