
### Regression Variants
###
### Leonard Mada
###
### draft v 0.1c


# - various variations of the classic regression;
# - minimise the square of the actual distance to the regression line;
# - added effects of scaling:
#  -- different scales on the 2 axes can degrade the performance;


############################

# dir.create(tempdir())

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
rliniar.gen = function(n, b1=3, b1.var=10, method=1, clip.neg=TRUE) {
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
	x.var = ifelse(x.var < 0, -x.var, x.var)
	y = x * (b1 + runif(n, -x.var, x.var))
	#
	y = y + runif(n, 0, 5)
	x = x + runif(n, 0, 2)
	#
	if(clip.neg) {
		y[y < 0] = 0
	}
	return(data.frame(x=x, y=y))
}
### Test the regressions
plot.test = function(n, b1.var=100, method=1, clip=TRUE, scale=TRUE) {
	x.df = rliniar.gen(n, b1.var=b1.var, method=method, clip.neg=clip)
	### Scale
	if(scale) {
		x = scale(x.df[,1])
		y = scale(x.df[,2])
	} else {
		x = x.df[,1]
		y = x.df[,2]
	}
	# Start plot
	plot(x, y)

	### Ordinary liniar regression
	r = solve.reg(x, y)
	print(r)
	print(dist.line(x, y, r[1], r[2]))
	curve(r[1] * x + r[2], add=T, col="red")

	### Liniar Distance to line
	r = solve.dist(x, y)
	print(r)
	print(dist.line(x, y, r[1], r[2]))
	curve(r[1] * x + r[2], add=T, col="blue")
	#
	invisible(x.df)
}
plot.all.reg = function(x, y, scale=FALSE, col=c("red", "blue")) {
	if(scale) {
		### Scale
		x = scale(x)
		y = scale(y)
	}
	plot(x, y)
	### Ordinary liniar regression
	r = solve.reg(x, y)
	print(r)
	d = dist.line(x, y, r[1], r[2])
	print(d)

	curve(r[1] * x + r[2], add=T, col=col[1])

	### Liniar Distance to line
	# - fails with non-constant variance / different metrics on x/y;
	# - may be useful with polynomial regression;
	r = solve.dist(x, y)
	print(r)
	d = dist.line(x, y, r[1], r[2])
	print(d)

	curve(r[1] * x + r[2], add=T, col=col[2])
}

##################
##################

b1.var = 0.4 # 100
x.df = rliniar.gen(100, b1.var=b1.var, method=1)
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

### Effects of Scaling

### Scale
x.sc = scale(x)
y.sc = scale(y)


plot(x.sc, y.sc)

### Ordinary liniar regression
r = solve.reg(x.sc, y.sc)

a = r[1]
b = r[2]
r
dist.line(x.sc, y.sc, a, b)

curve(a * x + b, add=T, col="red")

### Liniar Distance to line
r = solve.dist(x.sc, y.sc)

a = r[1]
b = r[2]
r
dist.line(x.sc, y.sc, a, b)

curve(a * x + b, add=T, col="blue")


#######

plot.test(200, -0.3, method=3)

##############
### Comparison
summary(lm(y ~ x))



####################

###
n = 100
x = 1:n
x = sample(x, n, replace=TRUE)
y = c(x + runif(n, 0, 2), x + runif(n, 0, 5), x - runif(n, 5, 30))
x = c(x,x,x)
x = x + runif(n, -2, 2)

plot.all.reg(x, y)

plot.all.reg(x, y, scale=T)


###
n = 100
x = 1:n
x = sample(x, n, replace=TRUE)
x.b = (1 + runif(n, -1/x, 1/x)) * x
y = c(x.b + runif(n, 10, 15), x.b + runif(n, -5, 5), x.b - runif(n, 10, 15))
x = c(x,x,x)
x = x + runif(n, -2, 2)

plot.all.reg(x, y)

plot.all.reg(x, y, scale=T)


####################
####################

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

