


###
# Maths 505: This wacky integral has a beautiful result
# https://www.youtube.com/watch?v=a2ZPqB2Syfo
# Note: series expansion + Dirichlet kernel;

###
integrate(\(x) sin(x + sin(2*x)) * exp(cos(2*x)) / sin(x), 0, pi/2)
integrate(\(x) Im(exp(1i*x + exp(2i*x))) / sin(x), 0, pi/2)
pi/2 * exp(1)


###
integrate(\(x) cos(x + sin(2*x)) * exp(cos(2*x)) / cos(x), 0, pi/2)
pi/2 * exp(-1)


### Dirichlet kernel:
k = 7; x = sqrt(5)
sin((2*k+1)*x) / sin(x) # ==
2*sum(cos(2*seq(k)*x)) + 1;

