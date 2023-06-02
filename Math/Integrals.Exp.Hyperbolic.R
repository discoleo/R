



### I( x^2 / cosh(x^2)^2 )
# Maths 505: This EPIC integral is the best thing you'll see today
# https://www.youtube.com/watch?v=aQEJevhQjTQ

integrate(\(x) x^2 / cosh(x^2)^2, 0, Inf)
sqrt(pi/8)*(1 - sqrt(2)) * pracma::zeta(1/2)


###
p = sqrt(5)
integrate(\(x) x^p / cosh(x^2)^2, 0, Inf)
gamma(p/2 + 1/2) * pracma::zeta(p/2 - 1/2) * (1 - 2^(3/2 - p/2)) / 2^(p/2 - 1/2)

# TODO: more generalizations

