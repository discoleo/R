#########################
##
## Log: Complex Fractions


### Base:
p = sqrt(2); n = sqrt(5); k = sqrt(3);
integrate(\(x) x^p / (x^n + 1)^k, 0, Inf)
gamma((p+1)/n) * gamma(k - (p+1)/n) / gamma(k) / n

###
integrate(\(x) x^p * log(x^n + 1) / (x^n + 1)^k, 0, Inf)
gamma((p+1)/n) * gamma(k - (p+1)/n)*
	(digamma(k) - digamma(k - (p+1)/n)) / gamma(k) / n


#################
#################

###
n = 4; k = 1; # !!!
p = 1/sqrt(3);
integrate(\(x) Re(x^p * log(x^4 + 1) / (x^2 + 1i)), 0, Inf)
integrate(\(x) Im(x^p * log(x^4 + 1) / (x^2 + 1i)), 0, Inf)
(gamma((p+1 + 2*k)/n) * gamma(k - (p+1 + 2*k)/n)*
	(digamma(k) - digamma(k - (p+1 + 2*k)/n)) +
	- 1i*(gamma((p+1)/n) * gamma(k - (p+1)/n)*
	(digamma(k) - digamma(k - (p+1)/n)))) / gamma(k) / n

