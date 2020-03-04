
### Dynamical Systems

### Lab 1: Discrete Dynamical Systems
### Ch 14: Iterated Maps, p 347 [pdf 355]
### Ch 14.6: Examples ["Python Prigrams"], p 326 [pdf 385]
###
### Leonard Mada

# [draft v0.1]
# TODO:
# - Logistic Map;
# - Bifurcation Diagrams;


# helper function: plot window
plot.init = function(x.lim, y.lim, labels) {
	plot.new()
	plot.window(xlim = x.lim, ylim = y.lim)
	axis(1)
	axis(2)
	title(xlab=labels[1], ylab=labels[2])
}

#######################

### Example 14.a.) Tent Map

plot.tent = function(X, X1, X2, mu = 2) {
	plot.init(x.lim, y.lim, c("x", "T(x)"))
	lines(X1, mu*X1)
	lines(X2, mu*(1-X2))
	lines(X, X, col="blue") # bisector x = y
}


### actual trajectories

###
orbit = function(x0, iter, print=FALSE, mu = 2) {
	x = x0
	print(x)
	inputs = c()
	outputs = c()

	for(i in seq(2, iter)) {
		inputs = c(inputs, x)
		inputs = c(inputs, x)
		outputs = c(outputs, x)
		if(x <= 1/2) {
			x = mu * x
		} else if(x > 1/2) {
			x = mu - mu * x
		}
		outputs = c(outputs, x)
		if(print) {
			print(x)
		}
	}
	return(data.frame(inputs, outputs))
}

###
mu = 2 # 1.7

# Plot the tent function and line y=x.
x.lim = c(0, 1)
y.lim = c(0, mu/2)
X1 = c(0, 0.5)
X2 = c(0.5, 1)
X = seq(0, 1, length.out=200)

plot.tent(X, X1, X2, mu=mu)

###

x = 1/5
iter = 20
col.iter = "red"
x.df = orbit(x, iter, mu=mu)
lines(x.df, lwd=2, col=col.iter)

# + epsilon
x = 1/5 + 0.001
iter = 20
# col2rgb("pink")
col.iter = rgb(255, 100, 100, max = 255, alpha = 100, names = "pink50")
x.df = orbit(x, iter, mu=mu)
lines(x.df, lwd=2, col=col.iter)

x = 1/7 # 1/5 + 1/3
iter = 30
col.iter = "green"
x.df = orbit(x, iter)
lines(x.df, lwd=2, col=col.iter)

x = 1/8
iter = 30
col.iter = "purple"
x.df = orbit(x, iter)
lines(x.df, lwd=2, col=col.iter)


x = 1/11
iter = 30
col.iter = "darkgreen"
x.df = orbit(x, iter)
lines(x.df, lwd=2, col=col.iter)

# non-periodic
plot.tent(X, X1, X2)
#
x = 1/pi
iter = 50
col.iter = "darkgreen"
x.df = orbit(x, iter)
lines(x.df, lwd=2, col=col.iter)



#######################

### Example 14.b.) Bifurcation Diagram

# TODO

