

### Derivations


### Limits:

#
p = 1/sqrt(3); n = sqrt(5);
q = 1E-5;
(gamma((p+1)/n) * (digamma(q) - digamma((p+1)/n+q)) / gamma((p+1)/n+q) +
- gamma((p+2)/n) * (digamma(q) - digamma((p+2)/n+q)) / gamma((p+2)/n+q)) * gamma(q) / n;
#
(1/2*(pracma::psi(1, (p+1)/n) - pracma::psi(1, (p+2)/n)) +
	+ Euler * (digamma((p+1)/n) - digamma((p+2)/n)) +
	+ 1/2 * (digamma((p+1)/n)^2 - digamma((p+2)/n)^2)) / n +
	# TODO: ugly limit;
	- (pracma::psi(1, (p+1)/n) * gamma((p+2)/n+q) * gamma((p+1)/n) +
		- pracma::psi(1, (p+2)/n) * gamma((p+1)/n+q) * gamma((p+2)/n)) *
	q * gamma(q) / gamma((p+1)/n+q) / gamma((p+2)/n+q) / n;


# Derivation:
((digamma(q) - (q*pracma::psi(1, (p+1)/n) + digamma((p+1)/n))) * gamma((p+2)/n+q) * gamma((p+1)/n) +
- (digamma(q) - (q*pracma::psi(1, (p+2)/n) + digamma((p+2)/n))) * gamma((p+1)/n+q) * gamma((p+2)/n) ) *
	gamma(q) / gamma((p+1)/n+q) / gamma((p+2)/n+q) / n;
#
(digamma(q) * (gamma((p+2)/n+q) * gamma((p+1)/n) - gamma((p+1)/n+q) * gamma((p+2)/n)) +
- digamma((p+1)/n) * gamma((p+2)/n+q) * gamma((p+1)/n) +
+ digamma((p+2)/n) * gamma((p+1)/n+q) * gamma((p+2)/n) ) *
	gamma(q) / gamma((p+1)/n+q) / gamma((p+2)/n+q) / n +
- (pracma::psi(1, (p+1)/n) * gamma((p+2)/n+q) * gamma((p+1)/n) +
	- pracma::psi(1, (p+2)/n) * gamma((p+1)/n+q) * gamma((p+2)/n)) *
	q * gamma(q) / gamma((p+1)/n+q) / gamma((p+2)/n+q) / n;


# Sub 2:
(pracma::psi(1, (p+1)/n) * gamma((p+2)/n+q) * gamma((p+1)/n) +
	- pracma::psi(1, (p+2)/n) * gamma((p+1)/n+q) * gamma((p+2)/n)) *
	q * gamma(q) / gamma((p+1)/n+q) / gamma((p+2)/n+q) / n;
# TODO


# Sub 1:
(digamma(q) * (gamma((p+2)/n+q) * gamma((p+1)/n) - gamma((p+1)/n+q) * gamma((p+2)/n)) +
- digamma((p+1)/n) * gamma((p+2)/n+q) * gamma((p+1)/n) +
+ digamma((p+2)/n) * gamma((p+1)/n+q) * gamma((p+2)/n) ) *
	gamma(q) / gamma((p+1)/n+q) / gamma((p+2)/n+q) / n;
#
(digamma(q) * (gamma((p+2)/n+q) * gamma((p+1)/n) +
	- gamma((p+1)/n+q) * gamma((p+2)/n)) +
+ (digamma((p+2)/n) - digamma((p+1)/n)) * gamma((p+1)/n) * gamma((p+2)/n) ) *
	gamma(q) / gamma((p+1)/n+q) / gamma((p+2)/n+q) / n;
#
((- 1/q - Euler + pi^2/6 * q) * (gamma((p+2)/n+q) * gamma((p+1)/n) +
	- gamma((p+1)/n+q) * gamma((p+2)/n)) +
+ (digamma((p+2)/n) - digamma((p+1)/n)) * gamma((p+1)/n) * gamma((p+2)/n) ) *
	gamma(q) / gamma((p+1)/n+q) / gamma((p+2)/n+q) / n;
#
((digamma((p+2)/n) - digamma((p+1)/n)) * gamma((p+1)/n) * gamma((p+2)/n) +
	- (Euler + 1/q) * (gamma((p+2)/n+q) * gamma((p+1)/n) +
	- gamma((p+1)/n+q) * gamma((p+2)/n))
) * gamma(q) / gamma((p+1)/n+q) / gamma((p+2)/n+q) / n;
#
(1/2*(pracma::psi(1, (p+1)/n) - pracma::psi(1, (p+2)/n)) +
	+ Euler * (digamma((p+1)/n) - digamma((p+2)/n)) +
	+ 1/2 * (digamma((p+1)/n)^2 - digamma((p+2)/n)^2)) / n;


p = 1; n = 4;
#
((digamma((p+2)/n) - digamma((p+1)/n)) * gamma((p+1)/n) * gamma((p+2)/n) +
	- (Euler + 1/q) * (gamma((p+2)/n+q) * gamma((p+1)/n) +
	- gamma((p+1)/n+q) * gamma((p+2)/n))
) * gamma(q);
#
(1/2*(pracma::psi(1, (p+1)/n) - pracma::psi(1, (p+2)/n)) +
	+ Euler * (digamma((p+1)/n) - digamma((p+2)/n)) +
	+ 1/2 * (digamma((p+1)/n)^2 - digamma((p+2)/n)^2) ) * gamma((p+1)/n) * gamma((p+2)/n);


#########

### [old]

###
p = 1/sqrt(3); n = sqrt(5);
integrate(\(x) x^p * (1 - x) * log(1-x^n) / (1 - x^n), 0, 1)
# Limit:
q = 1E-4;
(gamma((p+1)/n) * (digamma(q) - digamma((p+1)/n+q)) / gamma((p+1)/n+q) +
- gamma((p+2)/n) * (digamma(q) - digamma((p+2)/n+q)) / gamma((p+2)/n+q)) * gamma(q) / n;
#
- ((pracma::psi(1, q) - pracma::psi(1, (p+1)/n+q) +
	- (digamma(q) - digamma((p+1)/n+q)) * digamma((p+1)/n+q)) *
	gamma((p+1)/n) / gamma((p+1)/n+q) +
- (pracma::psi(1, q) - pracma::psi(1, (p+2)/n+q) +
	- (digamma(q) - digamma((p+2)/n+q)) * digamma((p+2)/n+q)) *
	gamma((p+2)/n) / gamma((p+2)/n+q)
) * gamma(q) / digamma(q) / n;
#
(pracma::psi(1, q) * (gamma((p+1)/n) / gamma((p+1)/n+q) - gamma((p+2)/n) / gamma((p+2)/n+q)) +
- digamma(q) * (digamma((p+1)/n+q)*gamma((p+1)/n) / gamma((p+1)/n+q) +
	- digamma((p+2)/n+q)*gamma((p+2)/n) / gamma((p+2)/n+q)) +
(- pracma::psi(1, (p+1)/n+q) + digamma((p+1)/n+q) * digamma((p+1)/n+q)) *
	gamma((p+1)/n) / gamma((p+1)/n+q) +
- (- pracma::psi(1, (p+2)/n+q) + digamma((p+2)/n+q) * digamma((p+2)/n+q)) *
	gamma((p+2)/n) / gamma((p+2)/n+q)
) / n;
#
(pracma::psi(1, q) * (gamma((p+1)/n) / gamma((p+1)/n+q) - gamma((p+2)/n) / gamma((p+2)/n+q)) +
- digamma(q) * (digamma((p+1)/n+q)*gamma((p+1)/n) / gamma((p+1)/n+q) +
	- digamma((p+2)/n+q)*gamma((p+2)/n) / gamma((p+2)/n+q)) +
	+ pracma::psi(1, (p+2)/n) - pracma::psi(1, (p+1)/n) +
	+ digamma((p+1)/n)^2 - digamma((p+2)/n)^2) / n;

#
pracma::psi(1, (p+2)/n) - pracma::psi(1, (p+1)/n) +
	+ digamma((p+1)/n)^2 - digamma((p+2)/n)^2;

#
pracma::psi(1, q) * (gamma((p+1)/n) / gamma((p+1)/n+q) - gamma((p+2)/n) / gamma((p+2)/n+q)) +
- digamma(q) * (digamma((p+1)/n+q)*gamma((p+1)/n) / gamma((p+1)/n+q) +
	- digamma((p+2)/n+q)*gamma((p+2)/n) / gamma((p+2)/n+q));
(pracma::psi(1, q) - digamma(q) * digamma((p+1)/n+q)) * gamma((p+1)/n) / gamma((p+1)/n+q) +
- (pracma::psi(1, q) - digamma(q) * digamma((p+2)/n+q)) * gamma((p+2)/n) / gamma((p+2)/n+q);
(pracma::psi(1, q) - digamma(q) * (q*pracma::psi(1, (p+1)/n) + digamma((p+1)/n))) *
	gamma((p+1)/n) / gamma((p+1)/n+q) +
- (pracma::psi(1, q) - digamma(q) * (q*pracma::psi(1, (p+2)/n) + digamma((p+2)/n))) *
	gamma((p+2)/n) / gamma((p+2)/n+q);
#
( (pracma::psi(1, q) - digamma(q) * digamma((p+1)/n)) * gamma((p+2)/n+q) * gamma((p+1)/n) +
- (pracma::psi(1, q) - digamma(q) * digamma((p+2)/n)) * gamma((p+1)/n+q) * gamma((p+2)/n)
) / (gamma((p+1)/n) * gamma((p+2)/n)) + # / gamma((p+1)/n+q) / gamma((p+2)/n+q) +
+ pracma::psi(1, (p+1)/n) - pracma::psi(1, (p+2)/n);
#
((digamma((p+2)/n) - digamma((p+1)/n)) * digamma(q) * gamma((p+1)/n) * gamma((p+2)/n) +
pracma::psi(1, q) * (gamma((p+2)/n+q) * gamma((p+1)/n) - gamma((p+1)/n+q) * gamma((p+2)/n))
) / (gamma((p+1)/n) * gamma((p+2)/n)) +
+ pracma::psi(1, (p+1)/n) - pracma::psi(1, (p+2)/n);

