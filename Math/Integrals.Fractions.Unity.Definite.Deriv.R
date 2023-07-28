

### Derivation: Specific Cases

# Note:
# - generalized formulas available;
# - uses function intUnityI01WX, see file:
#   Integrals.Fractions.Unity.Definite.R


####################

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
# simplification:
id = 1:2; cs = cos(2*id*pi/n); ch = cos(id*pi/n);
2*sum(cs * log(ch)) / n + pi/sin(pi/n) / (2*n)

# [old]
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
# simplification:
id = 1:3; cs = cos(2*id*pi/n); ch = cos(id*pi/n);
2*sum(cs * log(ch)) / n + pi/sin(pi/n) / (2*n)

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


#######################
#######################

#######################
### x^p / (x^n + 1) ###
#######################

### I( x / (x^n + 1) )

# - explicit derivation of detailed formulas;

### Odd Powers:
n = 5
cs = cos(seq((n-1)/2)*2*pi/n);
sn = sin(seq((n-1)/2)*2*pi/n);
cs2 = 2*cs^2 - 1; sn2 = 2*sn*cs;
#
integrate(function(x) x/(x^n + 1), 0, 1)
- 1/(2*n)*log(2) - 1/n * sum( cs2*log(cs + 1) ) +
	- 2/n * sum( sn2 * atan((1 + cs)/sn) ) +
	+ 2/n * sum( sn2 * (pi/2 - 2*seq((n-1)/2)*pi/n) );
# alternate:
csH = cos(seq((n-1)/2)*pi/n);
- 2/n * sum( cs2*log(csH) ) +
	- 2/n * sum( sn2 * atan((1 + cs)/sn) ) +
	+ pi/n * sum(sn2) - 4*pi/n^2 * sum(sn2*seq((n-1)/2));
# simplification:
- 2/n * sum( cs2*log(csH) ) +
	- 2/n * sum( sn2 * atan((1 + cs)/sn) ) +
	# sign of pi/n / sin() depends on p:
	+ pi/n * sum(sn2) + pi/n / sin(2*pi/n);


# Indefinite Integral:
- 1/n*log(x+1) - 1/n * sum( cs2*log(x^2 + 2*cs*x + 1) ) +
	- 2/n * sum( sn2 * atan((x - cs)/sn) )

# Fraction Decomposition:
n = 5
x = 3^(1/4) # test value;
cs = cos(seq((n-1)/2)*2*pi/n)
cs2 = 2*cs^2 - 1; # cos(2*...);
#
x/(x^n + 1)
x/n/(x+1) + 2/n * sum( (cs*x^2 + x) / (x^2 + 2*cs*x + 1) )
- 1/n/(x+1) - 2/n * sum( (cs2*x + cs) / (x^2 + 2*cs*x + 1) )
- 1/n/(x+1) - 2/n * sum( (cs2*x + cs*cs2) / (x^2 + 2*cs*x + 1) ) +
	- 2/n * sum( (cs - cs*cs2) / (x^2 + 2*cs*x + 1) )


### Odd Powers:
n = 7
p = 3; # any: ODD or EVEN integer;
cs = cos(seq((n-1)/2)*2*pi/n);
sn = sin(seq((n-1)/2)*2*pi/n);
cs2 = cos(seq((n-1)/2)*2*pi*(p+1)/n);
sn2 = sin(seq((n-1)/2)*2*pi*(p+1)/n);
csH = cos(seq((n-1)/2)*pi/n);
sign = if(p %% 2 == 0) - 1 else 1;
#
integrate(function(x) x^p/(x^n + 1), 0, 1)
# simplification:
- sign * 2/n * sum( cs2*log(csH) ) +
	+ pi/(2*n) / sin((p+1)*pi/n);

# Derivation:
(- 2/n * sum( cs2*log(csH) ) + pi/n * sum(sn2) +
	- 2/n * sum( sn2 * atan((1 + cs)/sn) ) ) * sign +
	+ pi/n / sin((p+1)*pi/n);
- sign * 2/n * sum( cs2*log(csH) ) +
	+ pi/(2*n) / sin((p+1)*pi/n);


### Using Digamma Function:
n = 7
p = sqrt(3)
- pi/(2*n) / sin((p+1)*pi/n) + integrate(function(x) x^p/(x^n + 1), 0, 1)$value;
- pi/(2*n) / sin((p+1)*pi/n) + integrate(function(x) x^p/(x^n + 1), 1, Inf)$value;
int.FrU01(n, p=p) - pi/(2*n) / sin((p+1)*pi/n);


###
n = 11
p = 3
integrate(function(x) x^p/(x^n + 1), 0, 1)
intUnityI01WX(n, p)


### Even Powers

###
n = 14
p = 2
integrate(function(x) x^p/(x^n + 1), 0, 1)
intUnityI01WX(n, p)

###
n = 14
p = 3
integrate(function(x) x^p/(x^n + 1), 0, 1)
intUnityI01WX(n, p)

###
n = 14
p = 4
integrate(function(x) x^p/(x^n + 1), 0, 1)
intUnityI01WX(n, p)


###
n = 16
p = 4
integrate(function(x) x^p/(x^n + 1), 0, 1)
intUnityI01WX(n, p)

