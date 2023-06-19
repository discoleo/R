

### Integrals: Log & Exp


Euler   = 0.57721566490153286060651209008240243079;
Catalan = 0.915965594177219015054603514;
# Note:
# Catalan = - I(log(x)/(x^2 + 1), lower=0, upper=1)


### I( exp(-x) * log(x)^2 )
# 1.) Maths 505: Another awesome integral with a beautiful result
# https://www.youtube.com/watch?v=wUY6TJxfFTk
# 2.) Dr Peyam: A mathematically stunning formula
# https://www.youtube.com/watch?v=hmYYYhQF1RM
# 3.) Michael Penn: Euler's other constant
# https://www.youtube.com/watch?v=UEqI9GKYozU
# [the basic proof]


###
integrate(\(x) exp(-x) * log(x), 0, Inf)
- Euler

###
n = sqrt(3)
integrate(\(x) exp(-n*x) * log(x), 0, Inf)
- Euler/n - log(n)/n


### Pow = 2
integrate(\(x) exp(-x) * log(x)^2, 0, Inf)
pi^2/6 + Euler^2

### Pow = 3
integrate(\(x) exp(-x) * log(x)^3, 0, Inf)
- 2 * pracma::zeta(3) - Euler*(pi^2/2 + Euler^2)
# psi2(1)*gamma(1) + 2*psi1(1)*d(gamma)(1) + psi(1)*d2(gamma)(1)
- 2* pracma::zeta(3) - 2*pi^2/6 * Euler - Euler*(pi^2/6 + Euler^2)

