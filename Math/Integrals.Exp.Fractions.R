

# Note:
# erf(1) = 2 * (pnorm(1, sd = 1/sqrt(2)) - 0.5);
# erf(1) = 2 * (pnorm(sqrt(2)) - 0.5);
# I(exp(-x^2), 0, 1) = erf(1) * sqrt(pi)/2;


################

###
# Maths 505: A crazy yet perfect integral
# https://www.youtube.com/watch?v=regnh-HOR0c
# - using Feynman's technique;

### I( exp(-x^2) / (x^2 + 1) ) on [0, Inf]
integrate(\(x) exp(-x^2) / (x^2 + 1), 0, Inf)
exp(1) * (pi/2 - gamma(1/2) * integrate(\(x) exp(-x^2), 0, 1)$value)
exp(1) * (pi/2 - gamma(1/2) * (pnorm(sqrt(2)) - 0.5)*sqrt(pi))


###
integrate(\(x) exp(-x^3) / (x^3 + 1), 0, Inf)
exp(1) * (pi/3/sin(pi/3) - gamma(1/3) * integrate(\(x) x * exp(-x^3), 0, 1)$value)
exp(1) * (pi/3/sin(pi/3) - gamma(1/3)/2 * integrate(\(x) exp(-x^1.5), 0, 1)$value)


###
integrate(\(x) exp(-x^4) / (x^4 + 1), 0, Inf)
exp(1) * (pi/4/sin(pi/4) - gamma(1/4) * integrate(\(x) x^2 * exp(-x^4), 0, 1)$value)
exp(1) * (pi/4/sin(pi/4) - gamma(1/4)/3 * integrate(\(x) exp(-x^(4/3)), 0, 1)$value)

