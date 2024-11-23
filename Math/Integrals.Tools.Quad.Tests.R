


####################

### Helper Functions

source("Integrals.Tools.Quad.R")


####################

#############
### Tests ###
#############

### Test 1:
f = \(x) sin(x) / (x + sin(3*x)^5 + cos(5*x)^5)

# Interval:
lim = c(0, 3);

n = 81; # precision too low;
# n = 181; # n = 201; # n = 281; # as good as base-R;
tmp = intw(lim, n=n)

print(integrate(f, lim[1], lim[2]), 12)
print(integrate(f, lim[1], lim[2], rel.tol=1E-12), 12)
print(integrate(\(x) { x = mpfr(x, 240); as.numeric(f(x)); },
	lim[1], lim[2], rel.tol=1E-12), 12)
sum(tmp$w * f(tmp$x))
tmp = intw(lim, n=281); sum(tmp$w * f(tmp$x))

