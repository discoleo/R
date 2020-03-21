
### Dynamical Systems

### Lab 1: Discrete Dynamical Systems
### Ch 14: Iterated Maps, p 347 [pdf 356]
### Ch 14.1: Tent Map
### Ch 14.6: Examples ["Python Programs"], p 326 [pdf 385]
###
### Leonard Mada

# [draft v0.2.1]
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

# transparent color
# useful when the orbits overlap
t_col = function(colors, percent = 50, name = NULL) {
	#     colors = color names
	#    percent = % transparency
	#       name = an optional name for the color
	
	t.col = c()
	for(color in colors) {
		rgb.val = col2rgb(color)
		## Make new color using input color as base and alpha set by transparency
		t.col.new = rgb(rgb.val[1], rgb.val[2], rgb.val[3],
			max = 255, alpha = (100 - percent) * 255 / 100,
			names = paste(color, name, sep=""))

		## Save the color
		t.col = c(t.col, t.col.new)
	}
	
	invisible(t.col)
}

#######################

### Tent Map
# Compute
tent.f = function(x, mu=2, iter=0, digits=NA) {
	r = x
	iter = iter + 1
	for(i in 1:iter) {
		r = ifelse(r <= 1/2, r*mu, mu*(1-r))
		if( ! is.na(digits)) {
			# for numeric stability, when possible;
			r = round(r, digits)
		}
	}
	return(r)
}
all.tent = function(x0, mu, iter, digits=NA) {
	sapply(0:iter, function(iter) tent.f(x0, mu=mu, iter, digits))
}
# Plot
tent.plot = function(x1=0, x2=1, mu = 2) {
	plot.init(x.lim, y.lim, c("x", "T(x)"))
	x.mid = (x1 + x2) / 2
	X1 = c(x1, x.mid)
	X2 = c(x.mid, x2)
	X = c(x1, x2)
	lines(X1, mu*X1)
	lines(X2, mu*(1-X2))
	lines(X, X, col="blue") # bisector x = y
}
# composed tent function
ctent.plot = function(n=2, x1=0, x2=1, mu = 2, x.lim=c(x1, x2), y.lim=c(0, mu/2)) {
	plot.init(x.lim, y.lim, c("x", "T(x)"))
	#
	div = 2*n
	i = 0:div
	x.mid = ((div-i)*x1 + i*x2) / div
	X1 = x.mid
	f.val = mu*X1
	# TODO: proper & generalization
	# c(1/(2*mu), 1/2, 1-1/(2*mu))
	pow = c(1/(2*mu), 1/2, 1-1/(2*mu))
	pow = c(0, pow, 1)
	f.val = rep(mu/2, 2*n+1)
	f.val[2*(0:n) + 1] = 0
	#
	lines(pow, f.val)
	#
	X = c(x1, x2)
	lines(X, X, col="blue") # bisector x = y
}

ctent.plot(2, mu=3/2)


### actual trajectories
orbit = function(x0, iter, mu = 2, print=FALSE, digits=NA) {
	x = x0
	print(x)
	inputs = c()
	outputs = c()

	for(i in seq(2, iter)) {
		inputs = c(inputs, x, x)
		# inputs = c(inputs, x)
		outputs = c(outputs, x)
		if(x <= 1/2) {
			x = mu * x
		} else if(x > 1/2) {
			x = mu - mu * x
		}
		if( ! is.na(digits)) {
			x = round(x, digits)
		}
		outputs = c(outputs, x)
		if(print) {
			print(x)
		}
	}
	return(data.frame(inputs, outputs))
}

##########################

### Example 1: [pdf p 358]

x.lim = c(0, 1)
x1 = 0
x2 = 1

#############
### E 1.1 ###
mu = 1/2
y.lim = c(0, 0.3 + mu/2)

### Init
# Plot the tent function and line y=x
tent.plot(x1, x2, mu=mu)

# E 1.1.a
x = 1/4
iter = 20
col.iter = "red"
x.df = orbit(x, iter, mu=mu)
lines(x.df, lwd=2, col=col.iter)

# E 1.1.b
x = 1/2
iter = 20
col.iter = t_col("green", 60)
x.df = orbit(x, iter, mu=mu)
lines(x.df, lwd=2, col=col.iter)

# E 1.1.c
x = 3/4 # solution to exercise is incorrect
iter = 20
col.iter = t_col("purple", 60)
x.df = orbit(x, iter, mu=mu)
lines(x.df, lwd=2, col=col.iter)

# E 1.1.d (extra)
# mu = 1/2 has NO fixed points, except (0, 0)
x = 1/3
iter = 20
col.iter = t_col("darkgreen", 60)
x.df = orbit(x, iter, mu=mu)
lines(x.df, lwd=2, col=col.iter)
#
x = 2/3
iter = 20
col.iter = t_col("darkgreen", 60)
x.df = orbit(x, iter, mu=mu)
lines(x.df, lwd=2, col=col.iter)

legend(0.8, 0.5, legend = c(3/4, 1/2, 1/4, "", "2/3", "1/3"),
	fill=c(t_col("purple", 60), t_col("green", 60), "red", "transparent", "darkgreen", "darkgreen"))


#############
### E 1.2 ###
mu = 1
y.lim = c(0, 0.5 + mu/2)

### Init
# Plot the tent function and line y=x
tent.plot(x1, x2, mu=mu)

# E 1.2.a
x = 1/3 # fixed point
iter = 20
col.iter = "red"
x.df = orbit(x, iter, mu=mu)
lines(x.df, lwd=2, col=col.iter)

# E 1.2.b
x = 2/3
iter = 20
col.iter = "green"
x.df = orbit(x, iter, mu=mu)
lines(x.df, lwd=2, col=col.iter)

# E 1.2.c extra
x = 4/5
iter = 20
col.iter = "darkgreen"
x.df = orbit(x, iter, mu=mu)
lines(x.df, lwd=2, col=col.iter)


#############
### E 1.3 ###
mu = 3/2
y.lim = c(0, 0.25 + mu/2)

### Init
# Plot the tent function and line y=x
tent.plot(x1, x2, mu=mu)

# E 1.3.a
x = 3/5 # fixed point
iter = 20
col.iter = "red"
x.df = orbit(x, iter, mu=mu)
lines(x.df, lwd=2, col=col.iter)

# E 1.3.b
x = 6/13 # 4-cycle / period = 2
iter = 20
col.iter = "red"
x.df = orbit(x, iter, mu=mu)
lines(x.df, lwd=2, col=col.iter)
# 4-cycle
tent.f(6/13, mu)
tent.f(6/13, mu, 4)

# E 1.3.c
x = 1/3 # does NOT enter a cycle; (but resembles a 16-cycle)
iter = 30
x.df = orbit(x, iter, mu=mu)
len = length(x.df[,1])
len1 = 7
len2 = 7
lenc = len - len1 - len2
col.iter = c(rep(t_col("darkgreen"), len1), rep(t_col("green"), len2), rep(t_col("orange"), lenc))
# lines(x.df, lwd=2, col=col.iter)
invisible(lapply(2:len, function(id) lines(x.df[c(id-1, id),], lwd=2, col=col.iter[id]) ))
# 16-cycle
tent.f(1/3, mu, 14)
tent.f(1/3, mu, 30)


#############
### E 1.4 ###
mu = 2
y.lim = c(0, mu/2)

### Init
# Plot the tent function and line y=x
tent.plot(x1, x2, mu=mu)

# E 1.4.a
x = 1/3 # enters fixed point 2/3
iter = 20
col.iter = "red"
x.df = orbit(x, iter, mu=mu)
lines(x.df, lwd=2, col=col.iter)

# E 1.4.b
x = 1/5 # enters 4-cycle
iter = 20
col.iter = t_col("orange")
x.df = orbit(x, iter, mu=mu)
lines(x.df, lwd=2, col=col.iter)

# E 1.4.c
x = 1/7 # enters 6-cycle
iter = 20
col.iter = t_col("springgreen")
x.df = orbit(x, iter, mu=mu)
lines(x.df, lwd=2, col=col.iter)

# E 1.4.c "extra"
x = 1/9 # enters 8-cycle
iter = 30
col.iter = t_col("darkgreen")
x.df = orbit(x, iter, mu=mu)
lines(x.df, lwd=2, col=col.iter)

# E 1.4.d
x = 1/11 # enters 10-cycle
iter = 20
col.iter = t_col("pink")
x.df = orbit(x, iter, mu=mu)
lines(x.df, lwd=2, col=col.iter)

# E 1.4.e "extra"
tent.plot(x1, x2, mu=mu)
#
x = 1/49 # enters 10-cycle
iter = 30
col.iter = t_col("red")
x.df = orbit(x, iter, mu=mu)
lines(x.df, lwd=2, col=col.iter)


tent.plot(x1, x2, mu=mu)
#
x = 31/32/5 # enters 10-cycle
iter = 30
x.df = orbit(x, iter, mu=mu)
len = length(x.df[,1])
len1 = 7
len2 = 7
lenc = len - len1 - len2
col.iter = c(rep(t_col("red"), len1), rep(t_col("orange"), len2), rep(t_col("blue"), lenc))
# lines(x.df, lwd=2, col=col.iter)
invisible(lapply(2:len, function(id) lines(x.df[c(id-1, id),], lwd=2, col=col.iter[id]) ))


##########################
### Example 2: [pdf p 359]

mu = 2
sapply(0:20, function(iter) tent.f(2/10, mu=mu, iter, 2))
sapply(0:23, function(iter) tent.f(21/100, mu=mu, iter, 2))
sapply(0:306, function(iter) tent.f(201/1000, mu=mu, iter, 3))

iter = 24
data.frame("x1"=all.tent(2/10, mu, iter, 2),
	"x2"=all.tent(21/100, mu, iter, 2),
	"x3"=all.tent(201/1000, mu, iter, 3),
	"x4"=all.tent(2001/10000, mu, iter, 4))


##########################
### Example 3: [pdf p 361]

# already covered in Example 1;


##########################
### Experiment
### for Definition 1
### [pdf p 364]

orbit.plot = function(x, iter, mu, col.iter, digits=NA) {
	x.df = orbit(x, iter, mu=mu, digits=digits)
	lines(x.df, lwd=2, col=col.iter)
	invisible(x.df)
}

mu = 3/2
y.lim = c(0, mu/2)

tent.plot(x1, x2, mu=mu)
iter = 10
col.iter = c(t_col("blue"), t_col(c("red", "green", "red", "pink"), 75))
col.iter = rev(col.iter)
digits = c(14, 14, 14, 14, 2)
for(i in 5:1) {
	iter.all = i + 3
	x.df = orbit.plot(1/mu^i/2, iter.all, mu, col.iter[i])
	x.dim = dim(x.df)
	# (to avoid numerical stability problems)
	# last.val = x.df[ x.dim[1], 1]
	# orbit.plot(last.val, 1, mu, col.iter[i], 2)
}
# unfortunately, the point (1/2, mu/2) is repelling
orbit.plot(1/2, iter, mu, col.iter[5], 2)

###
mu = 3/2
y.lim = c(0, mu/2)

tent.plot(x1, x2, mu=mu)
iter = 60
col.iter = t_col("red")
orbit.plot(1/2, iter, mu, col.iter)


###########################
###########################

###########################
### Example 14.a.) Tent Map

x.lim = c(0, 1)
x1 = 0
x2 = 1

mu = 2 # 1.7
y.lim = c(0, mu/2)

### Init
# Plot the tent function and line y=x
tent.plot(x1, x2, mu=mu)

x = 1/5
iter = 20
col.iter = "red"
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
tent.plot(x1, x2, mu=mu)
#
x = 1/pi
iter = 50
col.iter = "darkgreen"
x.df = orbit(x, iter)
lines(x.df, lwd=2, col=col.iter)

###
tent.plot(x1, x2, mu=mu)
#
x = 1/5
iter = 20
col.iter = "red"
x.df = orbit(x, iter, mu=mu)
lines(x.df, lwd=2, col=col.iter)
# + epsilon
x = 1/5 + 0.001
iter = 400
col.iter = rgb(255, 100, 100, max = 255, alpha = 100, names = "pink50")
x.df = orbit(x, iter, mu=mu, digits=3)
lines(x.df, lwd=2, col=col.iter)

sapply(0:40, function(iter) tent.f(21/100, mu=mu, iter, 2))
sapply(0:306, function(iter) tent.f(201/1000, mu=mu, iter, 3))
# can safely round to 3 digits;

#######################

### Example 14.b.) Bifurcation Diagram

# TODO

