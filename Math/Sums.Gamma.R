
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
pi^2/6 - (pracma::psi(1,3) - pracma::psi(1,1)/3) - 5/4

