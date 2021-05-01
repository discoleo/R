########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions


#######################

library(polynom)
library(pracma)

### helper Functions

rootn = function(r, n) {
	if(n %% 2 == 0) return(r^(1/n)); # TODO: complex?
	ifelse( (Im(r) == 0 & Re(r) < 0), - (-r)^(1/n), r^(1/n) )
}
unity = function(n=3, all=TRUE) {
	m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
	if(all) {
		m = m^(0:(n-1))
	}
	return(m)
}
roots.f = function(K, s, n=length(s)) {
	# roots for Class 1 polynomials;
	# s = includes s0;
	m = unity(n=n, all=T)
	k = rootn(K, n);
	order1 = n - 1;
	r = sapply(seq(n), function(id) sum(s * (k*m[id])^(0:order1)))
	return(r)
}
roots.cl2.f = function(s, n = length(s)) {
	# roots for Class 2 polynomials;
	m = unity(n+1, all=T)[-1]; # exclude 1;
	r = sapply(seq(n), function(id) sum(s * m[id]^(0:n)))
	r = round0(r)
}
### Multiplication
mult.p = function(p1, p2) {
	p.m = outer(p1, p2)
    p = as.vector(tapply(p.m, row(p.m) + col(p.m), sum))
	return(p)
}
### Multi-variable Multiplication
mult.pm = function(p1, p2) {
	# P1
	split.df = function(p.df) {
		p.l = lapply(colnames(p.df), function(name) p.df[,name]);
		names(p.l) = colnames(p.df);
		return(p.l);
	}
	if(is.data.frame(p1)) {
		p1 = split.df(p1);
	}
	p1.b0 = p1$coeff;
	p1 = p1[ ! names(p1) %in% "coeff"];
	# P2
	if(missing(p2)) {
		p2.b0 = p1.b0;
		p2 = p1;
	} else {
		if(is.data.frame(p2)) {
			p2 = split.df(p2);
		}
		p2.b0 = p2$coeff;
		p2 = p2[ ! names(p2) %in% "coeff"];
	}
	# helper
	prod.b0 = function(p1, p2=p1) outer(p1, p2);
	# Adjust Vars
	vars = unique(c(names(p1), names(p2)));
	len  = length(p1[[1]])
	for(v in vars[ ! vars %in% names(p1)]) p1[[v]] = rep(0, len);
	len  = length(p2[[1]])
	for(v in vars[ ! vars %in% names(p2)]) p2[[v]] = rep(0, len);
	# print(p1); print(p2);
	# Multiply
	p.m = lapply(vars, function(name) outer(p1[[name]], p2[[name]], function(i, j) i+j))
	p.b0 = prod.b0(p1.b0, p2.b0);
	p.l = lapply(p.m, as.vector);
	p.v = do.call(cbind, p.l)
	p.v = cbind(p.v, b0=as.vector(p.b0))
	p.r = aggregate(b0~., p.v, sum);
	colnames(p.r) = c(vars, "coeff");
	return(p.r);
}
pow.p = function(p, n=2) {
	if(n == 1) return(p);
	if(is.double(n) && (n == round(n))) n = as.integer(n);
	if( ! is.integer(n)) stop("n must be integer!")
	# Multiply
	# TODO: vectorize: as.integer(intToBits(3));
	p.r = NULL;
	p.pow = p;
	while (n > 0) {
		print(n)
		if (n %% 2 == 1) {
			if(is.null(p.r)) p.r = p.pow else p.r = mult.pm(p.r, p.pow);
		}
		if(n == 1) break;
        p.pow = mult.pm(p.pow);
        n = n %/% 2;
    }
	x.name = names(p)[1];
	id = order(p.r[,x.name], decreasing=TRUE);
	p.r = p.r[id,];
	return(p.r);
}

# round to 0
round0 = function(m, tol=1E-7) {
	m[abs(Re(m)) < tol & abs(Im(m)) < tol] = 0
	isZero = (Re(m) != 0) & (abs(Re(m)) < tol)
	if(any(isZero)) {
		m[isZero] = complex(re=0, im=Im(m[isZero]))
	}
	isZero = (Im(m) != 0) & (abs(Im(m)) < tol)
	if(any(isZero)) {
		m[isZero] = Re(m[isZero])
	}
	return(m)
}
round0.p = function(p, tol=1E-7) {
	p = round0(as.vector(p), tol=tol)
	class(p) = "polynomial"
	return(p)
}

### Other
perm.gen = function(x) {
	len = length(x)
	id = seq(len)
	id.m = outer(id, id, function(i, j) ((i+j+1) %% len + 1))
	p.m = x[id.m]
	dim(p.m) = dim(id.m)
	p.m
}

### Solvers:

### Simple systems:
solve.En = function(x, max.perm=0, n=4, duplicates=FALSE) {
	id = 1:length(x)
	if(max.perm == 1) {
		id.gr = matrix(id, nrow=1)
	} else {
		id.l = rep(list(id), n)
		id.gr = expand.grid(id.l)
		if( ! duplicates) {
			isDuplic = apply(id.gr, 1, function(id.val) any(duplicated(id.val)))
			id.gr = id.gr[ ! isDuplic , ]
		}
		if(max.perm > 0) {
			max.perm = min(max.perm, nrow(id.gr));
			id.gr = head(id.gr, n=max.perm);
		}
	}
	sol = if(n == 4) cbind(
			x1=x[id.gr[,1]], x2=x[id.gr[,2]], x3=x[id.gr[,3]], x4=x[id.gr[,4]])
		else cbind(x[id.gr[,1]], x[id.gr[,2]], x[id.gr[,3]])
	sol.names = if(n >= 4) paste0("x", seq(n)) else
		if(n == 3) c("x", "y", "z") else c("x", "y");
	colnames(sol) = sol.names;
	return(sol);
}
solve.EnAll = function(m, max.perm=0, n=4) {
	# generates ncol(m) * (nrow(m)!) root combinations/permutations!
	l = lapply(seq(ncol(m)), function(id) solve.En(as.vector(m[,id]), max.perm=max.perm, n=n));
	do.call(rbind, l)
}

### decomposed polynomial systems
solve.S = function(S, R, b=0) {
	# generic solver (based on existing S = x+y+z)
	b2 = if(length(b) > 1) b[2] else 0; # Ext A2;
	b3 = if(length(b) > 2) b[3] else 0; # Ext A3;
	x = sapply(S, function(x) roots(c(1, -x, R[2] - b2*x, - R[3] + b3*x)))
	len = length(S)
	S = matrix(S, ncol=len, nrow=3, byrow=T)
	yz = R[3]/x - b3
	yz.s = S - x
	# TODO: robust (when necessary)
	# Note: this simple method is NOT robust!
	yz.d = sqrt(yz.s^2 - 4*yz)
	y = (yz.s + yz.d) / 2
	z = yz.s - y
	cbind(as.vector(x), as.vector(y), as.vector(z))
}

solve.mS = function(S, b=0) {
	# generic solver (based on existing S = x+y+z)
	# S = cbind(S, E2, E3)
	b2 = if(length(b) > 1) b[2] else 0; # Ext A2;
	b3 = if(length(b) > 2) b[3] else 0; # Ext A3;
	x = sapply(seq(nrow(S)), function(id) roots(c(1, -S[id, 1], S[id, 2] - b2*S[id, 1], - S[id, 3] + b3*S[id, 1])))
	E2 = S[,2]; E3 = S[,3]; S = S[,1]; len = length(S)
	S  = matrix(S, ncol=len, nrow=3, byrow=T)
	E3 = matrix(E3, ncol=len, nrow=3, byrow=T)
	yz = E3/x - b3
	yz.s = S - x
	# TODO: robust (when necessary)
	yz.d = sqrt(yz.s^2 - 4*yz)
	y = (yz.s + yz.d) / 2
	z = yz.s - y
	cbind(as.vector(x), as.vector(y), as.vector(z))
}

### Print

# Print multi-variable Poly
print.monome = function(name, p) {
	v = p[,name];
	v.r = rep("", length(v));
	v.r[v > 1] = paste0(name, "^", v[v > 1]);
	v.r[v == 1] = name;
	return(v.r);
}
print.p = function(p) {
	id.coeff = match("coeff", colnames(p));
	coeff = p[,id.coeff]; p = p[, - id.coeff];
	p.str = sapply(colnames(p), print.monome, p=p);
	paste.nonempty = function(str, collapse="*") {
		str = str[nchar(str) > 0]
		paste(str, collapse=collapse)
	}
	p.str = apply(p.str, 1, paste.nonempty);
	coeff.str = as.character(coeff);
	hasCoeff = (coeff != 1);
	p.str[hasCoeff] = paste(coeff.str[hasCoeff], p.str[hasCoeff], sep = "*");
	return(paste(p.str, collapse=" + "));
}



#######################
#######################

#############
### Tests ###
#############

### Multi-variable Multiplication

# (x^3 + b1*x - R)^3
p = list(
	x = c(3,1,0),
	b1 = c(0,1,0),
	R = c(0,0,1),
	coeff = c(1,1,-1)
)

### Test
mult.pm(p)

p.v = pow.p(p, 3)
p.v

print.p(p.v[,c(2,3,4,1)])

