

constEuler = 0.57721566490153286060651209008240243079


zeta = function(n) {
	(1 + 1/(2^(n-1) - 1)) * pracma::eta(n);
}


### sum( (zeta(n) - 1) / n )
# Maths 505: A Beautiful Riemann Zeta Sum
# https://www.youtube.com/watch?v=MJn4Fwv3_tU

n = 100
sum(sapply(seq(2, n), \(n) (pracma::zeta(n) - 1) / n))
1 - constEuler


###
n = 200
sum(sapply(seq(2, n), \(n) pracma::zeta(n) / n))
log(n)



### sum( zeta(2*n) * x^(2*n) )
# Michael Penn: Some identities involving the Riemann-Zeta function.
# https://www.youtube.com/watch?v=2W2Ghi9idxM

n = 100
sum(sapply(seq(n), \(n) pracma::zeta(2*n) / 4^n))
1/2

#
sum(sapply(seq(n), \(n) pracma::zeta(2*n) / 2^n))
(1 - pi/sqrt(2)/tan(pi/sqrt(2))) / 2

#
x = 1/sqrt(3)
sum(sapply(seq(n), \(n) pracma::zeta(2*n) * x^(2*n)))
(1 - pi*x/tan(pi*x)) / 2


### sum( (zeta(2*n) - 1) * x^(2*n) )
sum(sapply(seq(n), \(n) pracma::zeta(2*n) - 1))
3/4

#
x = 1/sqrt(3)
sum(sapply(seq(n), \(n) (pracma::zeta(2*n) - 1) * x^(2*n)) )
(3 - pi*x/tan(pi*x)) / 2 - 1/(1 - x^2)


### sum( zeta(4*n) - 1 )
# Michael Penn: More identities involving the Riemann-Zeta function!
# https://www.youtube.com/watch?v=bRdGQKwusiE

n = 500
sum(sapply(seq(n), \(n) zeta(4*n) - 1))
7/8 - pi/4 * (exp(2*pi) + 1) / (exp(2*pi) - 1)

### sum( zeta(6*n) - 1 )
n = 500
sum(sapply(seq(n), \(n) zeta(6*n) - 1))
x = cos(2*pi/3) + 1i*c(1,-1)*sin(2*pi/3)
1 + 1/4 - 1/3 - sum(pi*x/tan(pi*x)) / 6


###
# Michael Penn: Another Riemann-Zeta function identity.
# https://www.youtube.com/watch?v=AmWvTN6cPz8

n = 500
sum(sapply(seq(n), \(n) pracma::zeta(2*n) / ((n+1)*(2*n+1))) )
1/2
#
n = 2000
sum(sapply(seq(n), \(n) zeta(2*n) / ((n+1)*(2*n+1))) )
1/2


###
# Michael Penn: More Riemann Zeta function identities!!
# https://www.youtube.com/watch?v=JwxgwXUruRM

n = 500
sum(sapply(seq(n), \(n) zeta(2*n) - zeta(2*n+1)) )
1/2

###
n = 500
sum(sapply(seq(n), \(n) zeta(2*n+1) - 1) )
1/4

###
n = 1000
sum(sapply(seq(n), \(n) zeta(n+1) - 1) )
1

