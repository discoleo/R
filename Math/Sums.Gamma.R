
### Sums: Gamma Function


### Sum( gamma(n)^2 / gamma(2*n) )
# Maths 505: A seemingly impossible infinite series
# https://www.youtube.com/watch?v=xM4z0ncARlw
# - transformed to Beta-function;

id = 1:40

###
sum( gamma(id)^2 / gamma(2*id) )
2*pi / sqrt(27)

###
sum( gamma(id) * gamma(id+1) / gamma(2*id+1) )
pi / sqrt(27)

### Sum( gamma(n)^2 / gamma(2*n+1) )
sum( gamma(id)^2 / gamma(2*id+1) )
pracma::psi(1,1)/3
pi^2/18


###################
###################

# Maths 505: An awesome calculus result I cooked up
# https://www.youtube.com/watch?v=WoxHJyOOma4


### Sum( (-1) gamma(n + 1/4) / gamma(n + 7/4) )
id = seq(0, 80)
sum((-1)^id * gamma(id + 1/4) / gamma(id + 7/4))
2*sqrt(pi)
#
id = mpfr(seq(0, 81), 240)
sum(rep(c(1,-1), length(id)/2) * gamma(id + 1/4) / gamma(id + 7/4))
2*sqrt(pi)

# Based on:
integrate(\(x) sqrt(1 - x^4) / (x^4 + 1), 0, 1)
integrate(\(x) sqrt(1/tan(x) - tan(x)) / 2, 0, pi/4)
pi/4


### Gen 1:

### Sum( (-1) gamma(n/2 + 1/8) / gamma(n/2 + 11/8) )
id = mpfr(seq(0, 401), 240)
sum(rep(c(1,-1), length(id)/2) * gamma(id/2 + 1/8) / gamma(id/2 + 11/8))
sqrt(2) * beta(5/8, 1/4) / gamma(5/4);
4*sqrt(2) * gamma(5/8) / gamma(7/8);

# Based on:
integrate(\(x) (1 - x^8)^(1/4) / (x^4 + 1), 0, 1)
integrate(\(x) (1/tan(x)^2 - tan(x)^2)^(1/4) / 2, 0, pi/4)
integrate(\(x) sqrt(2)/4 * (cos(x) / sin(x)^2)^(1/4), 0, pi/2)
sqrt(2)/8 * beta(5/8, 1/4)

