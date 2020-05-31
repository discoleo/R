
### Regression Variants
###
### Leonard Mada
###
### draft v 0.1


# - various variations of the classic regression;
# - minimise the actual distance to the regression line;


############################

### helper functions

solve.reg = function(x, y) {
	x.m = mean(x)
	y.m = mean(y)
	dx = x - x.m
	dy = y - y.m

	a = sum(x * dy) / sum(x * dx)
	b = y.m - a * x.m

	c(a, b)
}
solve.dist = function(x, y) {
	x.m = mean(x)
	y.m = mean(y)
	dx = x - x.m
	dy = y - y.m

	b2 = sum(x.m*y.m - x*y)
	b1 = sum(y^2 - y.m^2 - x^2 + x.m^2) / b2

	a = -b1/2 + sqrt(b1^2/4 + 1)
	# TODO: - sqrt() when liniar coeff < 0;
	b = mean(y) - a * mean(x)
	return(c(a, b))
}
dist.line = function(x, y, a, b) {
	sum((y - a*x - b)^2) / (a^2 + 1)
}
dist.line2 = function(x, y, a, b) {
	# Test: calculations of dist.line()
	c = y + x/a
	y.l = (b + a^2 * c) / (a^2 + 1)
	x.l = (y.l - b) / a
	sum((x - x.l)^2 + (y - y.l)^2)
}
### Random Generator
rliniar.gen = function(n, b1=3, b1.var=10, method=1) {
	x = 1:n
	if(method == 1) {
		# Method 1:
		x.var = b1.var/x # b1.var = 10
	} else if(method == 2) {
		# Method 2:
		x.var = b1.var * sqrt(x) # b1.var = 0.4
	} else {
		# Method 3:
		x.var = b1.var * log(x+0.5) # b1.var = ...
	}
	#
	y = x * (b1 + runif(n, -x.var, x.var))
	#
	y = y + runif(n, 0, 5)
	x = x + runif(n, 0, 2)
	#
	return(data.frame(x=x, y=y))
}

############

x.df = rliniar.gen(100, b1.var=100, method=1)
x = x.df$x
y = x.df$y


plot(x, y)


###########

### Ordinary liniar regression
r = solve.reg(x, y)
a = r[1]
b = r[2]
r
dist.line(x, y, a, b)

curve(a * x + b, add=T, col="red")

### Liniar Distance to line
# - fails with non-constant variance;
# - may be useful with polynomial regression;
r = solve.dist(x, y)
a = r[1]
b = r[2]
r
dist.line(x, y, a, b)

curve(a * x + b, add=T, col="blue")


##############
### Comparison
summary(lm(y ~ x))



####################



##############

### Derivation
### Test
a * sum((y - a*x - b)^2) + (a^2 + 1)*sum(x * (y - a*x - b))
sum(a*(dy - a*dx)^2 + (a^2 + 1) * x * (dy - a*dx))
sum(a*dy^2 + a^3*dx^2 - 2*a^2*x*y + 2*a^2*x.m*y.m - a*(a^2 + 1)*x*dx + (a^2 + 1) *x*dy)
sum(a*dy^2 - (a^2-1)*x*y + 2*a^2*x.m*y.m - a*(x^2 - x.m^2) - (a^2 + 1)*x.m*y.m)
sum(a*dy^2 - (a^2-1)*x*y + (a^2-1)*x.m*y.m - a*(x^2 - x.m^2))
sum(a*y^2 - a*y.m^2 - (a^2-1)*x*y + (a^2-1)*x.m*y.m - a*(x^2 - x.m^2))
sum(a^2*(x.m*y.m - x*y) + a*(y^2 - y.m^2 - x^2 + x.m^2) - (x.m*y.m - x*y))

### TODO:
# - variant for polynomial regression;

