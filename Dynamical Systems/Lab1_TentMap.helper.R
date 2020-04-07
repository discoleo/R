
### Dynamical Systems

### Lab 1: Discrete Dynamical Systems
### Helper Functions
### cover:
### Ch 14: Iterated Maps, p 347 [pdf 356]
### Ch 14.1: Tent Map
### Ch 14.6: Examples ["Python Programs"], p 326 [pdf 385]
###
### Leonard Mada

# [draft v0.3.0]
# TODO:
# - generalization of formulas;
# - generalization of code/functions;
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

test.col = function(n = 10, offset=0) {
	label = paste("0:", n, " + offset: ", offset, sep="", collapse="")
	plot.init(c(0, n), c(0, n), c(label, label))
	all.col = colors()
	for(id in 0:(n*n - 1)) {
		rect(id %% n, id %/% n,
			 1 + id %% n, 1 + id %/% n, col = all.col[id + 1 + offset])
	}
}
# test.col(offset=match("coral4", colors()))

palette.gen = function(n=3, isTransparent=FALSE) {
	col = c("orange", "darkgreen", "purple", "coral4", "pink")
	if(isTransparent) {
		col = t_col(col)
	}
	return(col)
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
ctent.plot = function(iter=1, mu = 2, x1=0, x2=1,
		x.lim=c(x1, x2), y.lim=c(0, mu/2), addLines=TRUE) {
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
	if(addLines) {
		abline(v = pow, col=c("red", "blue"))
	}
	return(pow)
}

plot.empiric.tent = function(mu=2, iter=1, vlines=NULL, col=c("red", "blue"), n=256) {
	curve(vtent.f(x, iter, mu=mu), from=0, to=1, n=n) # empiric
	if( ! is.null(vlines)) {
		abline(v = vlines, col=col)
	}
}

mu = 3/2
curve(vtent.f(x, 2, mu=mu), from=0, to=1, n=256)
curve(vtent.f(x, 3, mu=mu), from=0, to=1, n=256)

mu = 1.7
iter = 3
vlines = ctent.plot(iter, mu=mu)
curve(vtent.f(x, iter, mu=mu), from=0, to=1, n=256)
abline(v = vlines, col=c("red", "blue"))


### Cycles
# only Points
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
		id.all = unique(x$n)
		# may break col[cyc.id]
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
# plot ALL Cycles
plot.all.cycles = function(mu, iter=2, x=NULL, p.text=NULL,
		col=palette.gen(iter+1, transparent=TRUE), p.col="green", p.cex=1.5, lwd=2) {
	# Cycles
	if(is.null(x)) {
		x = special.points(mu=mu, iter)
	}
	iter.mod = iter + 1
	if(is.null(p.text)) {
		p.text = paste("C", iter.mod, sep="", collapse="")
	}
	plot.cycles(x, p.text=p.text, p.cex=p.cex, col=p.col)
	id = 0
	# TODO: apply
	for(x.p in x$p) {
		x.df = orbit(x.p, iter, mu=c)
		lines(x.df, lwd=lwd, col=col[(id %/% iter.mod) + 1])
		id = id + 1
	}
	return(x)
}
# Fixed Points + evolution to F
plot.e.fixed = function(x, mu, iter=2, p.text=c("eF", "F"), col="red", p.cex=1.5, lwd=2) {
	plot.fixed(x, p.text=p.text, p.cex=p.cex, col=col)
	for(x.p in x) {
		x.df = orbit(x.p, iter, mu=mu)
		lines(x.df, lwd=lwd, col=col)
	}
}
plot.fixed = function(x, y=x, mu.v=NULL, iter=1, p.text="F", col="blue", p.cex=1.5, pch=19, jitter=0.04,
		conn.col=t_col("darkgreen"), ...) {
	points(x, y, pch=pch, cex=p.cex, col=col)
	len = length(x)
	isTextVector = (length(p.text) > 1)
	for(id in 1:len) {
		if(isTextVector && is.na(p.text[id])) {next;}
		p.text.var = if(isTextVector) p.text[id] else p.text;
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

### Cycle Points

# TODO: generalization of formulas;
special.points = function(mu, iter=2) {
	# Cn-1
	id = 0:iter
	iter_n = iter + 1
	x = mu / (1 + mu^iter_n)
	x1 = x * mu^id
	if(iter == 1) {
		return(data.frame(p=x1, n=1))
	}
	
	# other C's
	if(mu != 1) {
		div = 1 / (mu^iter_n - 1)
		iter_h = tail(id, -1) + 1
		mu_p = mu^iter_h
		x2 = div * (mu * mu_p - mu_p)
		x2 = c(x2, div * mu_p[iter])
	}
	# TODO
	if(iter == 2) {
		# C3-2
		if(mu != 1) {
			x = 1 / (mu^3 - 1)
			x2 = c(x * (mu^2 - mu), x * (mu^3-mu^2), x * (mu^3-mu))
			x = c(x1, x2)
			id = rep(c(1:iter), rep(iter + 1, iter))
			x.df = data.frame(p=x, n=id)
		} else {
			x.df = data.frame(p=x1, n=1)
		}
	} else if(iter == 3) {
		# C4: 2,3,4
		if(mu != 1) {
			x = 1 / (mu^4 - 1)
			m.coeff = c(mu^3 - mu^2, mu^4 - mu^3)
			x2 = c(x * (mu^2 - mu), x * m.coeff, x* (mu^4 - mu))
			# degenerate Cyc 2: x * (mu^4-mu^2)
			x3 = c(x * (mu^3-mu), x * (mu^4-mu^2))
			# x = (c^3 - c^2 + c)/(c^4+1)
			m.coeff = c(m.coeff, sum(m.coeff)) + mu
			div = 1/(mu^4 + 1)
			x4 = div * m.coeff
			x4 = sort(c(x4, mu * x4[1]))
			###
			x = c(x1, x2, x3, x3, x4)
			iter_n = iter + 1
			id = rep(c(1:iter_n), rep(iter_n, iter_n))
			# id = rep(c(1:3), rep(iter + 1, 3))
			x.df = data.frame(p=x, n=id)
		} else {
			# degenerate Cyc 2: x * (mu^4-mu^2)
			x3 = c(x * (mu^3-mu), x * (mu^4-mu^2))
			id = rep(1:2, rep(4,2))
			x.df = data.frame(p=c(x1, x3, x3), n=1)
		}
	} else {
		x.df = data.frame(p=x1, n=1)
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

clean = function(x, rm.x, digits=10) {
	x.rounded = round(x, digits)
	rem = round(rm.x, digits)
	doRemove = x.rounded %in% rem
	x = sort(x[ ! doRemove ])
	return(x)
}

solve.p3 = function(b) {
	# x^3 + b2*x^2 + b1*x + b0 = 0
	b = b / b[1]
	c = -b[3] / 3
	d = -b[4] / 2
	if(b[2] != 0) {
		# x = x0 - s
		s = b[2]/3
		c = c + s^2 # = c - s^2 + 2*b[2]*s/3 = - b[3] / 3;
		d = d - s^3 + b[3]*s/2 # d + (s^3 - b[2]*s^2 + b[3]*s)/2
	} else {
		s = 0
	}
	det = d^2 - c^3
	det = ifelse(det < 0, complex(re=0, im=sqrt(-det)), sqrt(det))
	#
	if(Im(det) != 0) {
		r = (d + det)^(1/3) + (d - det)^(1/3)
	} else {
		r = ifelse(d + det >= 0, (d + det)^(1/3), -(-d-det)^(1/3)) +
			ifelse(d - det >= 0, (d - det)^(1/3), -(-d+det)^(1/3))
	}
	return(Re(r - s))
}

#####################
#####################

png.name = function(name) {
	last = tail(name, 1)
	len = nchar(last) + 1
	if(len < 3 || substr(last, len - 3, len) != "png") {
		name = c(name, ".png")
	}
	paste(name, sep="", collapse="")
}

to.png = function(name="test", file=png.name(name), bg="white", pointsize=15, ...) {
	png(file=file, bg=bg, pointsize=pointsize, ...)
}

