
### Dynamical Systems

### Lab 1: Discrete Dynamical Systems
### Helper Functions
### cover:
### Ch 14: Iterated Maps, p 347 [pdf 356]
### Ch 14.1: Tent Map
### Ch 14.6: Examples ["Python Programs"], p 326 [pdf 385]
###
### Leonard Mada

# [draft v0.2.5]
# TODO:
# - Logistic Map;
# - Bifurcation Diagrams;


### Helper functions

# plot window
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
	
	return(t.col)
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
vtent.f = function(x, mu=2, iter=0, digits=NA) {
	r = x
	iter = iter + 1
	for(i in 1:iter) {
		isLess = (r <= 1/2);
		r[isLess] = mu * r[isLess];
		r[ ! isLess ] = mu * (1 - r[ ! isLess ])
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
tent.plot = function(x1=0, x2=1, mu = 2, x.lim=c(x1, x2), y.lim=c(0, mu/2)) {
	plot.init(c(x1, x2), y.lim, c("x", "T(x)"))
	x.mid = (x1 + x2) / 2
	X1 = c(x1, x.mid)
	X2 = c(x.mid, x2)
	X = c(x1, x2)
	lines(X1, mu*X1)
	lines(X2, mu*(1-X2))
	lines(X, X, col="blue") # bisector x = y
}
# composed tent function
ctent.plot = function(iter=1, mu = 2, x1=0, x2=1, x.lim=c(x1, x2), y.lim=c(0, mu/2)) {
	plot.init(x.lim, y.lim, c("x", "T(x)"))
	# TODO: proper & generalization
	# n=2: c(1/(2*mu), 1/2, 1-1/(2*mu))
	i = 1:iter
	pow = 1/(2 * mu^i)
	if(iter > 1) {
		tmp = c()
		for(ipos in 1:(iter-1)) {
			tmp = c(tmp, 2*head(pow, -ipos) - tail(pow, -ipos))
		}
		# exclude non-applicable values
		tmp = tmp[tmp <= 1/2]
		#
		tmp3 = c()
		if(iter > 2) {
			# TODO: full
			tmp3 = 2*(pow[1] - pow[2]) + pow[3]
		}
		pow = c(tmp, tmp3, pow)
	}
	pow = c(0, pow, 1/2, 1-pow, 1)
	pow = sort(pow)
	f.val = vtent.f(pow, iter = iter, mu=mu) # TODO: iter ???
	#
	lines(pow, f.val)
	#
	X = c(x1, x2)
	lines(X, X, col="blue") # bisector x = y
	return(pow)
}

mu = 3/2
curve(vtent.f(x, 2, mu=mu), from=0, to=1, n=256)
curve(vtent.f(x, 3, mu=mu), from=0, to=1, n=256)

mu = 1.7
iter = 3
vlines = ctent.plot(iter, mu=mu)
curve(vtent.f(x, iter, mu=mu), from=0, to=1, n=256)
abline(v = vlines, col=c("red", "blue"))


# Cycles
plot.cycles = function(x, y=x, p.text="C", col="red", p.cex=1.5, pch=19, jitter=0.04) {
	if(is.null(dim(jitter))) {
		x.jitter = 0
		y.jitter = jitter
	} else {
		x.jitter = jitter[,1]
		y.jitter = jitter[,2]
	}
	# Cycles
	if(is.null(dim(x))) {
		id.all = 0
	} else {
		id.all = 1:(dim(x)[2])
	}
	for(cyc.id in id.all) {
		cyc.x = if(cyc.id == 0) {x;} else {x$p[x$n == cyc.id];}
		cyc.y = if(cyc.id == 0) {y;} else {y$p[y$n == cyc.id];}
		cyc.col = if(length(col) == 1) {col;} else {col[cyc.id];}
		points(cyc.x, cyc.y, pch=pch, cex=p.cex, col=cyc.col)
		len = length(cyc.x)
		for(id in 1:len) {
			text(cyc.x[id] - x.jitter, cyc.y[id] - y.jitter, p.text, col=cyc.col)
		}
	}
}
# Fixed Points
plot.fixed = function(x, y=x, mu.v=NULL, iter=1, p.text="F", col="blue", p.cex=1.5, pch=19, jitter=0.04,
		conn.col=t_col("darkgreen"), ...) {
	points(x, y, pch=pch, cex=p.cex, col=col)
	len = length(x)
	isTextVector = (length(p.text) > 1)
	for(id in 1:len) {
		p.text.var = ifelse(isTextVector, p.text[id], p.text)
		text(x[id], y[id] - jitter, p.text.var, col=col)
	}
	#
	if( ! is.null(mu.v)) {
		len = length(mu.v)
		for(id in 1:len) {
			if(is.na(mu.v[id])) {
				next;
			}
			x.df = orbit(x[id], iter, mu=mu.v[id])
			lines(head(x.df, -1), lwd=2, col=conn.col, ...)
		}
	}
}


#######################

### actual trajectories
orbit = function(x0, iter=1, mu = 2, print=FALSE, digits=NA) {
	x = x0
	print(x)
	inputs = c()
	outputs = c()

	for(i in seq(0, iter)) {
		inputs = c(inputs, x, x)
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

################

### Fixed Points

special.points = function(mu, iter=2) {
	if(iter == 1) {
		x = mu / (1 + mu^2)
		x = c(x, x * mu)
		x.df = data.frame(p=x, n=1)
	} else if(iter == 2) {
		# C2-1
		x = mu / (1 + mu^3)
		x1 = c(x, x * mu, x * mu^2)
		# C2-2
		if(mu != 1) {
			x = c(1 / (mu^3 - 1))
			x2 = c(x * (mu^2 - mu), x * (mu^3-mu^2), x * (mu^3-mu))
			x = c(x1, x2)
			id = rep(c(1:iter), rep(iter + 1, iter))
			x.df = data.frame(p=x, n=id)
		} else {
			x.df = data.frame(p=x1, n=1)
		}
	} else {
		print("Not yet implemented!")
	}
	return(x.df)
}

########################################

### Quadratic Function
quad.f = function(x, c=1/4, iter=0, digits=NA) {
	r = x
	iter = iter + 1
	for(i in 1:iter) {
		r = r^2 + c
		if( ! is.na(digits)) {
			# for numeric stability, when possible;
			r = round(r, digits)
		}
	}
	return(r)
}

# Plot
quad.plot = function(x1=0, x2=1, c=1/4, x.lim=c(x1, x2), y.lim=c(c, max(x.lim)^2 + c)) {
	plot.init(x.lim, y.lim, c("x", "P(x)"))
	x.mid = (x1 + x2) / 2
	X = c(x1, x2)
	f = function(x) x^2 + c;
	curve(f(x), from=x1, to=x2, n=256, add=TRUE)
	
	lines(X, X, col="blue") # bisector x = y
}

quad.orbit = function(x0, iter=20, c=1/4, digits=NA, ymin=0, print=FALSE) {
	x = x0
	print(x)
	inputs = c(x) # x
	outputs = c(ymin) # y

	for(i in seq(2, iter)) {
		if(i > 2) {
			inputs = c(inputs, x, x)
			outputs = c(outputs, x)
		} else {
			# skip bisector
			inputs = c(inputs, x)
		}
		x = x^2 + c
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

###########

solve.p3 = function(b) {
	# x^3 + b2*x^2 + b1*x + b0 = 0
	b = b / b[1]
	c = -b[3] / 3
	d = -b[4] / 2
	if(b[2] != 0) {
		# x = x0 - s
		s = b[2]/3
		c = c - s^2 + 2*b[2]*s/3
		d = d + (s^3 - b[2]*s^2 + b[3]*s)/2
	} else {
		s = 0
	}
	det = d^2 - c^3
	det = ifelse(det < 0, complex(re=0,im=sqrt(-det)), sqrt(det))
	#
	if(Im(det) != 0) {
		r = (d + det)^(1/3) + (d - det)^(1/3)
	} else {
		r = ifelse(d + det >= 0, (d + det)^(1/3), -(-d-det)^(1/3))
			+ ifelse(d - det >= 0, (d - det)^(1/3), -(-d+det)^(1/3))
	}
	return(Re(r - s))
}

