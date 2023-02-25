

# qncubed3:  TWO WAYS to destroy this INSANE integral feat. @maths_505
# https://www.youtube.com/watch?v=Ry5q4NsZDx0

integrate(function(x) sin(x)^2 / (x^2*(x^2+1)), 0, Inf)
pi/4*(1 + exp(-2))


### Helper:
k = sqrt(3)
integrate(function(x) cos(k*x) / (x^2 + 1), 0, Inf, subdivisions=4096*2, rel.tol=1E-6)
pi/2*exp(-k)


###
k = sqrt(3)
b = sqrt(5)
integrate(function(x) cos(k*x) / (x^2 + b^2), 0, Inf, subdivisions=4096*2, rel.tol=1E-6)
pi/2*exp(-k*b)/b


### TODO
k = 1
integrate(function(x) sin(k*x) / (x^2 + 1), 0, Inf, subdivisions=4096*2, rel.tol=1E-5)


### Gen 1:
k = sqrt(3)
integrate(function(x) sin(k*x)^2 / (x^2 + 1), 0, Inf, subdivisions=4096, rel.tol=1E-6)
pi/4*(1 - exp(-2*k))


### Gen 2:
k = sqrt(3)
integrate(function(x) sin(k*x)^2 / (x^2*(x^2 + 1)), 0, Inf, subdivisions=1024, rel.tol=1E-8)
pi/4*(2*k - 1 + exp(-2*k))

# f = (1 - exp(2i*k*z)) / (z^2 * (z^2 + 1));
# Res at 1i = pi*(1 - exp(-2*k));


########################

### Numerical Stability:

k = sqrt(3)
start = 0; steps = 1000
r = sapply(seq(start, steps), function(id) {
	startI = 2*pi*id/k;
	endI = 2*pi*(id+1)/k;
	integrate(function(x) cos(k*x) / (x^2 + 1), startI, endI, subdivisions=4096*2, rel.tol=1E-8)$value
})
r = sum(r)
# FAILS miserably!
r + integrate(function(x) cos(k*x) / (x^2 + 1), 2*pi*(steps+1)/k, Inf, subdivisions=4096*2, rel.tol=1E-8)$value
# residual: should be < 1E-8;
print(r, 12)
print(pi/2*exp(-k), 12)


######################
######################

### "Inverse"

Catalan = 0.915965594177219015054603514;


###
integrate(function(x) x / sin(x), 0, pi/2)
2*Catalan

###
integrate(function(x) x^2 / sin(x), 0, pi/2)
2*pi*Catalan - 7/2*pracma::zeta(3)

###
integrate(function(x) x*(pi - x) / sin(x), 0, pi)
7*pracma::zeta(3)

###
pracma::integral(function(x) x*(pi - x)*(pi + x) / sin(x), -pi, pi)
21*pi * pracma::zeta(3)

