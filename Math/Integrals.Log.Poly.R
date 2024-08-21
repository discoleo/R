

### Integrals: Log(Poly)


################
################

### I( log(x^2 + 2*sin(a)*x + 1) / x ) on [0,1]
# Maths 505: A beautifully unusual integral calculation
# https://www.youtube.com/watch?v=Na5kpJY0s9g

a = asin(sqrt(3)/4)
integrate(\(x) log(x^2 - 2*sin(a)*x + 1) / x, 0, 1)
- a^2/2 - pi/2*a + pi^2 / 24

###
a = asin(1/pi)
integrate(\(x) log(x^2 - 2*sin(a)*x + 1) / x, 0, 1)
- a^2/2 - pi/2*a + pi^2 / 24

###
a = asin(1/pi)
integrate(\(x) log(x^2 + 2*sin(a)*x + 1) / x, 0, 1)
- a^2/2 + pi/2*a + pi^2 / 24

###
a = asin(pi + 0i)
integrate(\(x) log(x^2 + 2*Re(sin(a))*x + 1) / x, 0, 1)
- a^2/2 + pi/2*a + pi^2 / 24
