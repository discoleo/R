

### Constants
Catalan = 0.915965594177219015054603514;


##########################

### sum( 1 / (j^2 + 1)^p )

###
j = seq(0, 10000)
sum( 1 / (j^2 + 1)^2 )
pi^2/4 / cosh(pi)^2 + pi/4 / tanh(pi) + 1/2
# TODO


###
j = seq(10000)
sum( j / (j^2 + 1)^2 )
(pracma::psi(1, 1 + 1i) - pracma::psi(1, 1 - 1i)) * 1i/4
(pracma::psi(1, 1i) - pracma::psi(1, - 1i)) * 1i/4


###
k = sqrt(5)
j = seq(10000)
sum( j / (j^2 + k^2)^2 )
(pracma::psi(1, k*1i) - pracma::psi(1, - k*1i)) * 1i / (4*k)
# TODO: simplify ?

# Derivation:
# - generate Polygamma function;
k = sqrt(3)
j = seq(10000)
sum( j / (j^2 + k^2)^2 )
sum( (j + k*1i)^2 / (j^2 + k^2)^2 - (j - k*1i)^2 / (j^2 + k^2)^2) / 4i / k
sum( 1 / (j - k*1i)^2 - 1 / (j + k*1i)^2) / 4i / k
sum( 1/j^2 - 1 / (j + k*1i)^2 - (1/j^2 - 1 / (j - k*1i)^2) ) / 4i / k
(pracma::psi(1, k*1i) - pracma::psi(1, - k*1i)) * 1i/ (4*k)


### Psi Order 1:
k = sqrt(3)
(pracma::psi(1, k*1i) + pracma::psi(1, 1 - k*1i))
pi^2 / sin(pi*k*1i)^2
- 4*pi^2 / (exp(pi*k) - exp(-pi*k))^2


###
# Kölbig, K. S. (1996). The polygamma function psi^k(x) for x=1/4 and x=3/4.
# J. Comput. Appl. Math. 75 (1): 43–46.
# doi:10.1016/S0377-0427(96)00055-6

(pracma::psi(1, 1/4) - pracma::psi(1, 3/4))
16 * Catalan
#
pracma::psi(1, 1/4)
pi^2 + 8*Catalan
#
pracma::psi(1, 3/4)
pi^2 - 8*Catalan


### Psi Order 2:

### Psi(2, 1/3)
# see file: Integrals.AnalyticCont.R;
pracma::psi(2, 1/3)
- 26 * pracma::zeta(3) - pi^3 * cos(pi/3) / sin(pi/3)^3
#
pracma::psi(2, 2/3)
- 26 * pracma::zeta(3) + pi^3 * cos(pi/3) / sin(pi/3)^3

