########################
##
## Leonard Mada
## [the one and only]
##
## Integrals:
##   Log + Fraction Unity
##   Base-Formulas
##
## draft v.0.1a


### Base-Formulas:
# - for easy Reference;
# - are based on formulas in file:
#   Integrals.Fractions.Unity.Definite.R;
#   (but this is a better reference)


#########################

### Helper:
p = sqrt(3); n = sqrt(7); k = sqrt(5);
integrate(\(x) x^p / (x^n + 1)^k, 0, Inf)
gamma((p+1)/n) * gamma(k - (p+1)/n) / gamma(k) / n;

# D(p)
integrate(\(x) x^p * log(x) / (x^n + 1)^k, 0, Inf)
gamma((p+1)/n) * gamma(k - (p+1)/n) * (digamma((p+1)/n) - digamma(k - (p+1)/n)) / gamma(k) / n^2
# D2(p)
integrate(\(x) x^p * log(x)^2 / (x^n + 1)^k, 0, Inf)
gamma((p+1)/n) * gamma(k - (p+1)/n) / gamma(k) / n^3 *
	(pracma::psi(1, (p+1)/n) + pracma::psi(1, k - (p+1)/n) +
	+ (digamma((p+1)/n) - digamma(k - (p+1)/n))^2)
# D3(p)
integrate(\(x) x^p * log(x)^3 / (x^n + 1)^k, 0, Inf)
gamma((p+1)/n) * gamma(k - (p+1)/n) / gamma(k) / n^4 *
	(pracma::psi(2, (p+1)/n) - pracma::psi(2, k - (p+1)/n) +
	3*(pracma::psi(1, (p+1)/n) + pracma::psi(1, k - (p+1)/n)) *
	(digamma((p+1)/n) - digamma(k - (p+1)/n)) +
	+ (digamma((p+1)/n) - digamma(k - (p+1)/n))^3)


### D(k)
p = sqrt(3); n = sqrt(7); k = sqrt(5);
integrate(\(x) x^p * log(x^n + 1) / (x^n + 1)^k, 0, Inf)
gamma((p+1)/n) * gamma(k - (p+1)/n) / gamma(k) / n *
	(digamma(k) - digamma(k - (p+1)/n));

# Case: p = 0
integrate(\(x) log(x^n + 1) / (x^n + 1)^k, 0, Inf)
gamma(1/n) * gamma(k - 1/n) / gamma(k) / n *
	(digamma(k) - digamma(k - 1/n))


### D2(k)
integrate(\(x) x^p * log(x^n + 1)^2 / (x^n + 1)^k, 0, Inf)
gamma((p+1)/n) * gamma(k - (p+1)/n) / gamma(k) / n *
	(pracma::psi(1, k - (p+1)/n) - pracma::psi(1, k) +
	+ (digamma(k) - digamma(k - (p+1)/n))^2);

# Case: p = 0
integrate(\(x) log(x^n + 1)^2 / (x^n + 1)^k, 0, Inf)
gamma(1/n) * gamma(k - 1/n) / gamma(k) / n *
	(pracma::psi(1, k - 1/n) - pracma::psi(1, k) +
	+ (digamma(k) - digamma(k - 1/n))^2)

