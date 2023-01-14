########################
###
### Leonard Mada
### [the one and only]
###
### Exact Integration
### Polynomial Fractions: Unity
### Definite Integrals
###
### draft v.0.1b



### Definite Integrals

# I( 1 / (x^n + 1) )
# I( x^p / (x^n + 1) )

### Intervals:
# [0, Inf], [0, 1];

# - contains various Shortcut Formulas;
# - moved specific code to this file from file:
#   Integrals.Fractions.Unity.R;


######################

### Helper Functions

### I on [0, 1]
intUnityI01 = function(n) {
	if(n %% 2 == 0) {
		return(intUnityI01Even(n));
	}
	n2 = n %/% 2;
	sum(cos(2*seq(n2)*pi/n)*log(cos(seq(n2)*pi/n))) * 2 / n +
	+ sum(seq(n2)*sin(2*seq(n2)*pi/n)) * 2*pi / n^2;
}
intUnityI01Even = function(n) {
	if(n %% 2 == 1) return(intUnityI01(n));
	isNotM4 = n %% 4 != 0;
	id = seq(1, n, by=2);
	if(isNotM4) id = id[ - (n+2)/4];
	cs = cos(id*pi/n); sn = sin(id*pi/n);
	# Base-Formula:
	# fr = sum(2*cs*x + 2)/(x^2 + 2*cs*x + 1)) / n;
	# r = sum(cs*log(x^2 + 2*cs*x + 1)) + sum(2*atan((x - cs)/sn));
	r = sum(cs*log(cs + 1)) +
		# + sum(2*sn*atan((1 - cs)/sn)) - sum(2*sn*atan((0 - cs)/sn));
		+ sum(2*sn*atan((1 - cs)/sn)) - 2*pi*sum(id*sn)/n + pi*sum(sn);
	if(isNotM4) r = r + pi/2;
	return(r/n);
}

######################
######################

### Infinite Integrals
# - shortcut formulas are available;

n = 5
integrate(function(x) 1/(x^n + 1), lower=0, upper=Inf)
integrate(function(x) x^3/(x^n + 1), lower=0, upper=Inf)
pi / sin(pi/n) / n
#
integrate(function(x) x/(x^n + 1), lower=0, upper=Inf)
integrate(function(x) x^2/(x^n + 1), lower=0, upper=Inf)
pi / sin(2*pi/n) / n

###
n = 7
integrate(function(x) 1/(x^n + 1), lower=0, upper=Inf)
integrate(function(x) x^5/(x^n + 1), lower=0, upper=Inf)
pi / sin(pi/n) / n
#
integrate(function(x) x/(x^n + 1), lower=0, upper=Inf)
integrate(function(x) x^4/(x^n + 1), lower=0, upper=Inf)
pi / sin(2*pi/n) / n
#
integrate(function(x) x^2/(x^n + 1), lower=0, upper=Inf)
integrate(function(x) x^3/(x^n + 1), lower=0, upper=Inf)
pi / sin(3*pi/n) / n


############################
############################

### Special Finite Integrals
# - shortcut formulas are available;

###
n = 9
integrate(function(x) 1/(x^n + 1), 0, 1)
intUnityI01(n)

###
n = 15
integrate(function(x) 1/(x^n + 1), 0, 1)
intUnityI01(n)

###
n = 16
integrate(function(x) 1/(x^n + 1), 0, 1)
intUnityI01(n)


### Specific Cases
# Note: generalized formulas available;

###
n = 3
integrate(function(x) 1/(x^3 + 1), 0, 1)
2/n * cos(2*pi/3)*log(cos(pi/3)) +
	+ 2*pi*sin(2*pi/3)/n^2

# Fraction Decomposition:
1/n/(x+1) + 2/n * sum( (cos(2*pi/3)*x + 1) / (x^2 + 2*cos(2*pi/3)*x + 1) )
1/n/(x+1) + 1/n * cos(2*pi/3)*(2*x + 2*cos(2*pi/3)) / (x^2 + 2*cos(2*pi/3)*x + 1) +
	+ 2/n * sin(2*pi/3)^2 / (x^2 + 2*cos(2*pi/3)*x + 1)
# I on [0, 1] =>
2/n * cos(2*pi/3)*log(cos(pi/3)) +
	+ 2/n * sin(2*pi/3)*atan( (1 + cos(2*pi/3))/sin(2*pi/3) ) +
	# - 2/n * sin(2*pi/3)*atan( (0 + cos(2*pi/3))/sin(2*pi/3) );
	- 2/n * (pi/2 - 2*pi/3) * sin(2*pi/3);

### n = 5
# TODO: simplify
n = 5
integrate(function(x) 1/(x^n + 1), 0, 1)
2*cos(2*pi/n)*log(cos(pi/n)) / n +
	+ 2*cos(4*pi/n)*log(cos(2*pi/n)) / n +
	+ 2*pi*(sin(2*pi/n) + 2*sin(4*pi/n)) / n^2;
# useful ???
- 4*log(2)*cos(2*pi/n)/n +
	+ 4/n * (4*cos(pi/n)^4 - 5*cos(pi/n)^2 + 1)*log(cos(2*pi/n)) +
	+ 2*pi*sin(2*pi/n)*(2*cos(pi/n)^2 + 3*cos(2*pi/n)) / n^2;


###
n = 7
integrate(function(x) 1/(x^n + 1), 0, 1)
2*cos(2*pi/n)*log(cos(pi/n)) / n +
	+ 2*cos(4*pi/n)*log(cos(2*pi/n)) / n +
	+ 2*cos(6*pi/n)*log(cos(3*pi/n)) / n +
	+ 2*pi*(sin(2*pi/n) + 2*sin(4*pi/n) + 3*sin(6*pi/n)) / n^2;

#
2/n * cos(2*pi/n)*log(cos(pi/n)) +
	+ 2/n * cos(4*pi/n)*log(cos(2*pi/n)) +
	+ 2/n * cos(6*pi/n)*log(cos(3*pi/n)) +
	+ 2/n*sin(2*pi/n)*atan( (1 + cos(2*pi/n))/sin(2*pi/n) ) +
	+ 2/n*sin(4*pi/n)*atan( (1 + cos(4*pi/n))/sin(4*pi/n) ) +
	+ 2/n*sin(6*pi/n)*atan( (1 + cos(6*pi/n))/sin(6*pi/n) ) +
	- 2/n * (pi/2 - 2*pi/n) * sin(2*pi/n) +
	- 2/n * (pi/2 - 4*pi/n) * sin(4*pi/n) +
	- 2/n * (pi/2 - 6*pi/n) * sin(6*pi/n)


##################

### "Continuity" & Approximation
# - "quasi-analytic" continuity possible;

###
n = seq(3, 22)
x = sapply(n, intUnityI01)
# asymptotic to y = 1;
plot(n, x)
curve( (x + 1/3)/(x+1), add=T, col="orange")

### [old] Only Odd Powers
n = seq(3, 21, by=2)
x = sapply(n, intUnityI01)
# asymptotic to y = 1;
plot(n, x)
curve( (x + 1/3)/(x+1), add=T, col="orange")

###
n = c(seq(3,20), seq(21, 101, by=3))
x = sapply(n, intUnityI01)
# asymptotic to y = 1;
plot(n, x)
curve( (x + 1/3)/(x+1), add=T, col="orange")


######################

### Special Fractional Powers
n = 7
integrate(function(x) 1/(x^(1 - 1/(n+1)) + 1), 0, 1)
n + 1 - (n+1)*integrate(function(x) 1/(x^n + 1), 0, 1)$value;

###
n = 8
integrate(function(x) 1/(x^(1 - 1/(n+1)) + 1), 0, 1)
n + 1 - (n+1)*integrate(function(x) 1/(x^n + 1), 0, 1)$value;

