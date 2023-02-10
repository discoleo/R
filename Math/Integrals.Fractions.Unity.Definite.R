########################
###
### Leonard Mada
### [the one and only]
###
### Exact Integration
### Polynomial Fractions: Unity
### Definite Integrals
###
### draft v.0.1i



### Definite Integrals

# I( 1 / (x^n + 1) )
# I( x^p / (x^n + 1) )

### Intervals:
# [0, Inf], [0, 1];

# - contains various Shortcut Formulas;
# - moved specific code to this file from file:
#   Integrals.Fractions.Unity.R;

### Sections:
# A.) Helper Formulas:
# - Trigonometric identities;
# B.) Definite Integrals


######################

### Helper Functions

### I on [0, 1]
intUnityI01 = function(n) {
	if(n %% 2 == 0) {
		return(intUnityI01Even(n));
	}
	n2 = n %/% 2;
	id = seq(n2);
	sum(cos(2*pi*id/n)*log(cos(pi*id/n))) * 2 / n +
		# simplification:
		# + sum(id*sin(2*pi*id/n)) * 2*pi / n^2;
		+ pi/(2*n) / sin(pi/n);
}
intUnityI01Even = function(n) {
	if(n %% 2 == 1) return(intUnityI01(n));
	isNotM4 = n %% 4 != 0;
	id = seq(1, n, by=2);
	if(isNotM4) id = id[ - (n+2)/4]; # cs == 0;
	cs = cos(id*pi/n); # sn = sin(id*pi/n);
	# Base-Formula:
	# fr = sum(2*cs*x + 2)/(x^2 + 2*cs*x + 1)) / n;
	# r = sum(cs*log(x^2 + 2*cs*x + 1)) + sum(2*atan((x - cs)/sn));
	r = sum(cs*log(cs + 1)) +
		# + sum(2*sn*atan((1 - cs)/sn)) - sum(2*sn*atan((0 - cs)/sn));
		# simplification: sum(2*sn*atan((1 - cs)/sn));
		# simplification: - 2*pi*sum(id*sn)/n;
		+ pi/2 / sin(pi/n);
	return(r/n);
}
### x^p:
intUnityI01WX = function(m, p=1) {
	if(n %% 2 == 0) {
		return(intUnityI01EvenWX(n, p=p));
	}
	n2 = n %/% 2;
	id = seq(n2);
	cs = cos(2*pi*id*(p+1)/n); csH = cos(pi*id/n);
	sign = if(p %% 2 == 0) 1 else -1;
	pi/(2*n) / sin(pi*(p+1)/n) + sign * sum(cs*log(csH)) * 2 / n;
}
intUnityI01EvenWX = function(n, p = 1) {
	if(n %% 2 == 1) return(intUnityI01WX(n, p=p));
	id = seq(1, n, by=2);
	cs = cos(id*pi*(p+1)/n); csH = cos(pi*id/n);
	sign = if(p %% 2 == 0) 1 else -1;
	r = pi/2 / sin(pi*(p+1)/n) + sign * sum(cs*log(csH + 1));
	return(r/n);
}

# I( 1 / (x^(n/2) + 1) )
intUnityI01Half = function(n) {
	if(n %% 2 == 0) warning("Not yet implemented!");
	cs = cos(seq((n-1)/2)*2*pi/n);
	sn = sin(seq((n-1)/2)*2*pi/n);
	cs2 = 2*cs^2 - 1; sn2 = 2*sn*cs;
	#
	int = - 1/(2*n)*log(2) - 1/n * sum( cs2*log(cs + 1) ) +
		- 2/n * sum( sn2 * atan((1 + cs)/sn) ) +
		+ 2*pi/n * sum( sn2 * (1/2 - 2*seq((n-1)/2)/n) );
	return(2*int);
}
# I( 1 / (x^(n/p) + 1) )
# - works for p < n; added proper sign;
intUnityI01POld = function(n, p=3) {
	if(n %% 2 == 0) warning("Not yet implemented!");
	cs = cos(seq((n-1)/2)*2*pi/n);
	sn = sin(seq((n-1)/2)*2*pi/n);
	csp = cos(p*seq((n-1)/2)*2*pi/n);
	snp = sin(p*seq((n-1)/2)*2*pi/n);
	#
	int = - 1/(2*n)*log(2) - 1/n * sum( csp*log(cs + 1) ) +
		- 2/n * sum( snp * atan((1 + cs)/sn) ) +
		+ 2*pi/n * sum( snp * (1/2 - 2*seq((n-1)/2)/n) );
	if(p %% 2 == 1) int = - int;
	return(p*int);
}
# simplification:
intUnityI01P = function(n, p=3) {
	if(n %% 2 == 0) warning("Not yet implemented!");
	cs = cos(seq((n-1)/2)*2*pi/n);
	sn = sin(seq((n-1)/2)*2*pi/n);
	csp = cos(p*seq((n-1)/2)*2*pi/n);
	snp = sin(p*seq((n-1)/2)*2*pi/n);
	#
	int = - 1/(2*n)*log(2) - 1/n * sum( csp*log(cs + 1) ) +
		- 2/n * sum( snp * atan((1 + cs)/sn) ) +
		+ 2*pi/n * sum( snp/2 );
	if(p %% 2 == 1) int = - int;
	int = int + pi/sin(pi*p/n)/n;
	return(p*int);
}

### Plot
# sum( cos(pi*seq(...)*p/n) * log(cos(pi*seq(...)/n)) )
plot.logcos = function(n, dx = 1/8, len=129, points=TRUE,
		type = "l", col.px = "red", title=TRUE, ...) {
	p = seq(-1 + dx, n - 1 - dx, length.out=len);
	r = sapply(p, function(p) {
		integrate(function(x) x^p / (x^n + 1), 0, 1)$value;
	})
	r = r - pi / sin(pi*(p+1)/n) / (2*n);
	plot(p, r, type=type, ...);
	if(title) title(paste0("n = ", n));
	# Explicit Points
	if(points) {
		isEven = (n %% 2 == 1)
		if(isEven) {
			# even seq:
			id = 2 * seq(1, (n-1)/2);
		} else {
			id = seq(1, n, by=2);
			# csH = 1 + cos(pi*id/n);
		}
		logcsH = log(cos(pi*id/(2*n)));
		px = seq(0, n-2, by=1);
		py = sapply(px, function(p) {
			sign = if(p %% 2 == 0) 1 else -1;
			sign * sum(cos(pi*id*(p+1)/n) * logcsH);
		})
		points(px, 2*py / n, col=col.px);
		return(py);
	}
	invisible()
}


####################
####################

#################
### Section A ###
#################

### Useful Formulas
# - Trigonometric identities;

### n = ODD
n = 11
p = 3
#
id = seq((n-1)/2)
sum(id * sin(2*pi*p*id/n)) # ==
(-1)^((p %% 2) + 1) * n/4/sin(pi*p/n)

### n = ODD
n = 11
p = 4
#
id = seq((n-1)/2)
sum(id * sin(2*pi*p*id/n)) # ==
(-1)^((p %% 2) + 1) * n/4/sin(pi*p/n)

### n = 4*k + 2
k = 3
p = 1
#
n = 4*k + 2;
id = seq(1, n, by=2);
sum(id * sin(pi*p*id/n)) # ==
(-1)^((p %% 2) + 1) * n/2/sin(pi*p/n)

### n = 4*k
# - same formula;
k = 3
p = 1
#
n = 4*k;
id = seq(1, n, by=2);
sum(id * sin(pi*p*id/n)) # ==
(-1)^((p %% 2) + 1) * n/2/sin(pi*p/n)


##############

### sum(sn)
n = 10 # EVEN!
id = seq(1, n, by=2);
sn = sin(pi*id/n);
sum(sn) - 1/sin(pi/n) # = 0


######################
######################

###################
###  Section B  ###
###  Integrals  ###
###################

### Infinite Integrals
# - shortcut formulas are available;

### I( x^k/(x^n + 1) )

###
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


### Non-Standard Powers

###
n = sqrt(5)
#
integrate(function(x) 1/(x^n + 1), lower=0, upper=Inf)
pi / sin(pi/n) / n
#
integrate(function(x) x/(x^n + 1), lower=0, upper=Inf)
pi / sin(2*pi/n) / n


###
n = 3; p = sqrt(2);
#
integrate(function(x) x^p/(x^n + 1), lower=0, upper=Inf)
pi / sin(pi*(p+1)/n) / n

###
n = 5; p = sqrt(3);
integrate(function(x) x^p/(x^n + 1), lower=0, upper=Inf)
pi / sin(pi*(p+1)/n) / n


###
n = sqrt(11); p = sqrt(3);
integrate(function(x) x^p/(x^n + 1), lower=0, upper=Inf)
pi / sin(pi*(p+1)/n) / n


############################
############################

### Special Finite Integrals
# - on [0, 1];
# - shortcut formulas are available;

###
n = 9
integrate(function(x) 1/(x^n + 1), 0, 1)
intUnityI01(n)

###
n = 14
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

###
n = 18
integrate(function(x) 1/(x^n + 1), 0, 1)
intUnityI01(n)


### Derivation: Specific Cases
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

###
plot.logcos(9)
plot.logcos(10)


######################

### Special Fractional Powers
n = 7
integrate(function(x) 1/(x^(1 - 1/(n+1)) + 1), 0, 1)
n + 1 - (n+1)*integrate(function(x) 1/(x^n + 1), 0, 1)$value;

###
n = 8
integrate(function(x) 1/(x^(1 - 1/(n+1)) + 1), 0, 1)
n + 1 - (n+1)*integrate(function(x) 1/(x^n + 1), 0, 1)$value;



#######################
#######################

#######################
### x^p / (x^n + 1) ###
#######################

### I( x / (x^n + 1) )

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


### TODO:
n = 7
p = sqrt(3)
- pi/(2*n) / sin((p+1)*pi/n) + integrate(function(x) x^p/(x^n + 1), 0, 1)$value;
- pi/(2*n) / sin((p+1)*pi/n) + integrate(function(x) x^p/(x^n + 1), 1, Inf)$value;
# TODO:
# - formula to compute: 0.06055876;


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


#####################

#############
### Other ###
#############

n = 8
integrate(function(x) 1/(1-x^n)^(1+1/n), lower=0, upper=1/2^(1/n))
# == 1

n = 9
integrate(function(x) 1/(1-x^n)^(1+1/n), lower=0, upper=1/2^(1/n))
# == 1

n = 10
integrate(function(x) 1/(1-x^n)^(1+1/n), lower=0, upper=1/2^(1/n))
# == 1

n = 11
integrate(function(x) 1/(1-x^n)^(1+1/n), lower=0, upper=1/2^(1/n))
# == 1


####################

###
p = 1 - sqrt(2)
n = sqrt(11)
k = sqrt(3)
integrate(function(x) x^p / (x^n+1)^(1/k), lower=0, upper=Inf)
gamma((p+1)/n)*gamma(1/k - (p+1)/n) / gamma(1/k) / n


###
p = sqrt(5) - 2
n = sqrt(11)
k = sqrt(3)
integrate(function(x) x^p / (x^n+1)^(1/k), lower=0, upper=Inf)
IInf = gamma((p+1)/n)*gamma(1/k - (p+1)/n) / gamma(1/k) / n
print(IInf)
#
- IInf/2 + integrate(function(x) x^p / (x^n+1)^(1/k), lower=0, upper=1)$value
- IInf/2 + integrate(function(x) x^p / (x^n+1)^(1/k), lower=1, upper=Inf)$value
# TODO: find formula for -0.3248496;
