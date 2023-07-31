########################
###
### Leonard Mada
### [the one and only]
###
### Integrals: Logarithms
### Log-Fractions: Old Derivations


###########
### OLD ###
###########

# - everything is old!
# - the generalized formulas are in:
#   Integrals.Log.Fractions.R;


##############

##############
### Powers ###

### [old]
integrate(function(x) log(x)/(x^2 + 1)^(3/2), 0, Inf)
- log(2)

###
integrate(function(x) log(x)/(x^2 + 1)^2, 0, Inf)
- pi / 4

###
integrate(function(x) log(x)/(x^2 + 1)^(5/2), 0, Inf)
- 2/3 * log(2) - 1/3

###
integrate(function(x) log(x)/(x^2 + 1)^3, 0, Inf)
- pi / 4

###
k = 4
integrate(function(x) log(x)/(x^2 + 1)^k, 0, Inf)
- 23/30 * gamma(1/2)*gamma(k-1/2)/gamma(k)

###
k = 5
integrate(function(x) log(x)/(x^2 + 1)^k, 0, Inf)
- 88/(3*5*7) * gamma(1/2)*gamma(k-1/2)/gamma(k)


#######################

################
### Radicals ###
################

### x^p * log(x) / (x^n + 1)^(1/k)

# TODO:
# - log-Term in result-formula is NOT yet generalized!
#   (only n, p = Integers)
# - cleanup;
# - full generalization *IS* available;

### k = 2
n = 6
integrate(function(x) log(x)/(x^n + 1)^(1/2), 0, Inf)
- (pi*cos(2*pi/n)/sin(2*pi/n) - 2*log(1/2)) * gamma(1/n) * gamma(1/2 - 1/n) / gamma(1/2) / n^2

#
n = 6; p = 1
integrate(function(x) x^p * log(x)/(x^n + 1)^(1/2), 0, Inf)
- (pi*cos(2*(p+1)*pi/n)/sin(2*(p+1)*pi/n) + 2*log(1/2)) *
	gamma((p+1)/n) * gamma(1/2 - (p+1)/n) / gamma(1/2) / n^2


###
n = 8
integrate(function(x) log(x)/(x^n + 1)^(1/2), 0, Inf)
- (pi*cos(2*pi/n)/sin(2*pi/n) - 2*sqrt(2)*log(tan(pi/8))) *
	gamma(1/n) * gamma(1/2 - 1/n) / gamma(1/2) / n^2
# log-Term: still NOT generalized!
p = 1
integrate(function(x) x^p * log(x)/(x^n + 1)^(1/2), 0, Inf)
- (pi*cos(2*pi*(p+1)/n)/sin(2*pi*(p+1)/n) - 2*sqrt(2)*log(tan(pi*(p+1)/8))) *
	gamma((p+1)/n) * gamma(1/2 - (p+1)/n) / gamma(1/2) / n^2


### n = EVEN
n = 10
p = 0; # p = Integer!
integrate(function(x) x^p * log(x)/(x^n + 1)^(1/2), 0, Inf)
id = seq(1, n, by=2);
- (pi*cos(2*pi*(p+1)/n)/sin(2*pi*(p+1)/n) - 4*sum(cos(2*pi*(p+1)*id/n)*log(cos(pi*id/(2*n))))) *
	gamma((p+1)/n) * gamma(1/2 - (p+1)/n) / gamma(1/2) / n^2


### n = ODD
n = 7
p = 0; # p = Integer!
integrate(function(x) x^p * log(x)/(x^n + 1)^(1/2), 0, Inf)
id = seq(2, n, by=2);
- (pi*cos(2*pi*(p+1)/n)/sin(2*pi*(p+1)/n) - 4*sum(cos(2*pi*(p+1)*id/n)*log(cos(pi*id/(2*n))))) *
	gamma((p+1)/n) * gamma(1/2 - (p+1)/n) / gamma(1/2) / n^2


#############
### k = 3 ###
#############

### n = ODD
n = 7
p = 0; # p = Integer!
k = 3;
integrate(function(x) x^p * log(x)/(x^n + 1)^(1/k), 0, Inf)
# see generalized Formulas;

###
n = 4
p = 0; # p = Integer!
k = 3;
integrate(function(x) x^p * log(x) / (x^n + 1)^(1/k), 0, Inf)
id = seq(2, n, by=2); th = pi*(p+1)/n;
(pi*cos(pi/n-pi/3)/sin(pi/n) + 3*log(3)/2 +
	+ (2*(sqrt(3)+1)*acosh(26) - 4*(3+sqrt(3))*atanh(1/sqrt(3))) / 4) *
	gamma((p+1)/n) * gamma(1/k - (p+1)/n) / gamma(1/k) / n^2


###
n = 9
p = 0; # p = Integer!
k = 3; # TODO: ???
integrate(function(x) x^p * log(x) / (x^n + 1)^(1/k), 0, Inf)
id = seq(2, n, by=2); th = pi*(p+1)/n;
- gamma((p+1)/n) * gamma(1/k - (p+1)/n) / gamma(1/k) / n^2 *
	(pi/sin(2*pi/9)/2 - 2*(cos(pi/9)*log(2*cos(5*pi/9) + 1) +
		+ cos(2*pi/9)*log(cos(4*pi/9)) + cos(4*pi/9)*log(cos(pi/9))) +
	- 2*log(2)*(cos(2*pi/9) + cos(4*pi/9)))


### Other
k = sqrt(5); n = 2*k;
integrate(function(x) log(x)/(x^n + 1)^(1/k), 0, Inf)
# == 0

### [k == 2]
k = 2; p = sqrt(3);
n = 4*(p+1);
integrate(function(x) x^p * log(x)/(x^n + 1)^(1/k), 0, Inf)
# == 0

### [k == 3]
n = 10;
k = 3; p = 2/3;
integrate(function(x) x^p * log(x)/(x^n + 1)^(1/k), 0, Inf)
# == 0

### [k == 3]
k = 3;
p = 5/6; n = 2*k*(p+1);
integrate(function(x) x^p * log(x)/(x^n + 1)^(1/k), 0, Inf)
# == 0

