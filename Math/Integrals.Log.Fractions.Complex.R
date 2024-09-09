#########################
##
## Log: Complex Fractions


### Values

Euler   = 0.57721566490153286060651209008240243079;
Catalan = 0.915965594177219015054603514;


################
################

### Base:
p = sqrt(2); n = sqrt(5); k = sqrt(3);
integrate(\(x) x^p / (x^n + 1)^k, 0, Inf)
gamma((p+1)/n) * gamma(k - (p+1)/n) / gamma(k) / n

###
integrate(\(x) x^p * log(x^n + 1) / (x^n + 1)^k, 0, Inf)
gamma((p+1)/n) * gamma(k - (p+1)/n) *
	(digamma(k) - digamma(k - (p+1)/n)) / gamma(k) / n

###
p = sqrt(2); n = sqrt(5); k = sqrt(3);
b = sqrt(2)
integrate(\(x) x^p * log(x^n + b^n) / (x^n + b^n)^k, 0, Inf)
gamma((p+1)/n) * gamma(k - (p+1)/n) *
	(digamma(k) - digamma(k - (p+1)/n) + n*log(b)) /
	(gamma(k) * n * b^(n*k - p - 1))

# Special Case: k = 1
p = 1/sqrt(2); n = sqrt(5);
b = sqrt(2)
integrate(\(x) x^p * log(x^n + b^n) / (x^n + b^n), 0, Inf)
- gamma((p+1)/n) * gamma(1 - (p+1)/n) *
	(digamma(1 - (p+1)/n) + Euler - n*log(b)) /
	(n * b^(n - p - 1))


###########

###
p = sqrt(2); n = sqrt(5); k = sqrt(3);
a = sqrt(3) - sqrt(2); b = sqrt(2)
integrate(\(x) x^p * log(x^n + a^n) / (x^n + b^n)^k, 0, Inf)
# TODO


#################
#################

###
p = - 1/sqrt(3);
integrate(\(x) x^p * (pi/2 - atan(x)) / (x^2 + 1), 0, Inf)
integrate(\(x) x^p * Re(log((x+1i)/(x-1i)) / 2i) / (x^2 + 1), 0, Inf)
- gamma((p+1)/2) * gamma(1 - (p+1)/2) *
	(digamma(1 - (p+1)/2) + Euler) / 4i +
+ (gamma((p+2)/2) * gamma(1 - (p+2)/2) *
	(digamma(-p/2) + Euler) +
	- 1i*(gamma((p+1)/2) * gamma(1 - (p+1)/2) *
	(digamma(1 - (p+1)/2) + Euler))) / 4 +
	- gamma(p+1) * gamma(-p) / 2 *
	((digamma(-p) + Euler - log(1i)) / (-1i)^p +
	(digamma(-p) + Euler - log(-1i)) / 1i^p );

###
p = - 1/sqrt(3);
integrate(\(x) Re(x^p * log(x - 1i) / (x^2 + 1)), 0, Inf)
integrate(\(x) Im(x^p * log(x - 1i) / (x^2 + 1)), 0, Inf)
(gamma((p+2)/2) * gamma(1 - (p+2)/2) *
	(digamma(-p/2) + Euler) +
	- 1i*(gamma((p+1)/2) * gamma(1 - (p+1)/2) *
	(digamma(1 - (p+1)/2) + Euler))) / 4i +
	- gamma(p+1) * gamma(-p) / 2i *
	((digamma(-p) + Euler - log(1i)) / (-1i)^p +
	(digamma(-p) + Euler - log(-1i)) / 1i^p );


###
n = 2; k = 1; # !!!
p = - 1/sqrt(3);
integrate(\(x) Re(x^p * log(x - 1i) / (x + 1i)), 0, Inf)
integrate(\(x) Im(x^p * log(x - 1i) / (x + 1i)), 0, Inf)
(gamma((p+1 + k)/n) * gamma(k - (p+1 + k)/n) *
	(digamma(k) - digamma(k - (p+1 + k)/n)) +
	- 1i*(gamma((p+1)/n) * gamma(k - (p+1)/n) *
	(digamma(k) - digamma(k - (p+1)/n)))) / gamma(k) / n +
	+ gamma(p+1) * gamma(-p) *
	(digamma(-p) + Euler - log(1i)) / 1i^(- p);


###
n = 2; k = 1; # !!!
p = - 1/sqrt(3);
integrate(\(x) Re(x^p * log(x^2 + 1) / (x + 1i)), 0, Inf)
integrate(\(x) Im(x^p * log(x^2 + 1) / (x + 1i)), 0, Inf)
(gamma((p+1 + k)/n) * gamma(k - (p+1 + k)/n)*
	(digamma(k) - digamma(k - (p+1 + k)/n)) +
	- 1i*(gamma((p+1)/n) * gamma(k - (p+1)/n)*
	(digamma(k) - digamma(k - (p+1)/n)))) / gamma(k) / n

# TODO: log(x - 1i) / (x + 1i)


###
n = 4; k = 1; # !!!
p = 1/sqrt(3);
integrate(\(x) Re(x^p * log(x^4 + 1) / (x^2 + 1i)), 0, Inf)
integrate(\(x) Im(x^p * log(x^4 + 1) / (x^2 + 1i)), 0, Inf)
(gamma((p+1 + 2*k)/n) * gamma(k - (p+1 + 2*k)/n)*
	(digamma(k) - digamma(k - (p+1 + 2*k)/n)) +
	- 1i*(gamma((p+1)/n) * gamma(k - (p+1)/n)*
	(digamma(k) - digamma(k - (p+1)/n)))) / gamma(k) / n

