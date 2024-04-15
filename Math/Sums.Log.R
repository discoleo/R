

### Series with Log

# TODO

Euler   = 0.57721566490153286060651209008240243079;
Catalan = 0.915965594177219015054603514;


#################
#################

id = seq(1200000)

### Sum( (-1)^n * sin(pi*b*n) * log(n) / n )
b = 1/3;
- sum((-1)^id * sin(b*id*pi) * log(id) / id)
(pi/2 - pi*b/2)*log(pi) - pi*log(gamma(1/2 + b/2)) +
	- pi*b/2 * (log(2) + Euler) - pi/2*log(cos(pi*b/2));


#################

### Sum( (-1)^n * log(2*n+1) / (2*n+1) )
id = seq(1200000)
sum((-1)^id * log(2*id+1) / (2*id+1))
pi/4 * log(pi) - pi*log(gamma(3/4)) - pi/4 * Euler


# Transformed:
id = seq(2400000)
sum((-1)^id * log(id+1/2) / (id+1/2))
pi*log(gamma(1/4)/gamma(3/4)/gamma(1/2)) + pi/2 * digamma(1/2) - digamma(1/2) - Euler

pi/2 * log(pi) + 2*log(2) - 2*pi*log(gamma(3/4)) +
	- log(2)*(digamma((1/2+1)/2) - digamma(1/2/2)) / 2 - pi/2 * Euler;
pi/2 * log(pi/2) + 2*log(2) - 2*pi*log(gamma(3/4)) - pi/2 * Euler;
pi/2 * log(pi) + 2*log(2) - 2*pi*log(gamma(3/4)) - pi/2 * (log(2) + Euler);

