



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

