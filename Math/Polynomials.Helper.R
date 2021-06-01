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

### Round to 0
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

### Root
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

sort.sol = function(sol, useRe=TRUE, ncol=1) {
	if(useRe) {
		id = order(abs(sol[,ncol]), Re(sol[,ncol]));
	} else {
		id = order(abs(sol[,ncol]));
	}
	return(sol[id,]);
}

### Polynomials

### Multiplication
mult.p = function(p1, p2) {
	p.m = outer(p1, p2)
    p = as.vector(tapply(p.m, row(p.m) + col(p.m), sum))
	return(p)
}
### Multi-variable Multiplication
mult.all.pm = function(p) {
	if( ! is.list(p)) stop("p must be a list of polynomials");
	len = length(p);
	pR = p[[1]];
	for(id in seq(2, len)) {
		p2 = p[[id]];
		if(is.numeric(p2)) {
			pR = mult.sc.pm(pR, p2);
		} else {
			pR = mult.pm(pR, p2);
		}
	}
	return(pR);
}
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
pow.pm = function(p, n=2) {
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
mult.sc.pm = function(p, s, coeff.name="coeff") {
	# Multiplication by scalar
	if(is.data.frame(p)) {
		p[ , coeff.name] = p[ , coeff.name] * s;
	} else if(is.list(p)) {
		p[[coeff.name]] = p[[coeff.name]] * s;
	} else stop("p must be a polynomial!")
	return(p);
}
reduce.pm = function(p) {
	# remove Monomes with coeff == 0;
	id = which(p$coeff != 0)
	if(is.data.frame(p)) {
		return(p[id, ]);
	}
	p = lapply(p, function(m) m[id]);
	return(p);
}
align.pm = function(p1, p2, align.names=TRUE) {
	p1 = reduce.pm(p1); p2 = reduce.pm(p2);
	n1 = names(p1); n2 = names(p2);
	### Coefficients
	n1 = n1[ ! n1 %in% "coeff"];
	n2 = n2[ ! n2 %in% "coeff"];
	xc = intersect(n1, n2);
	xall = union(n1, n2);
	pad.pm = function(p, vnew) {
		if(is.data.frame(p)) {
			p[, vnew] = 0;
		} else {
			len = length(p$coeff);
			zero = rep(0, len);
			for(nn in vnew) {
				p[[nn]] = zero;
			}
		}
		return(p);
	}
	# add missing variables
	if(length(xall) != length(xc)) {
		### p1
		n1new = n2[ ! n2 %in% n1];
		p1 = pad.pm(p1, n1new);
		### p2
		n2new = n1[ ! n1 %in% n2];
		p2 = pad.pm(p2, n2new);
	}
	#
	if(align.names) {
		id = match(names(p1), names(p2));
		list(p1=p1, p2=p2[id]);
	} else {
		list(p1=p1, p2=p2);
	}
}
add.pm = function(p1, p2) {
	l = align.pm(p1, p2);
	p1 = l[[1]]; p2 = l[[2]];
	n1 = names(p1); n2 = names(p2);
	### to DF
	id = match(n2, n1);
	p = rbind(as.data.frame(p1), as.data.frame(p2)[,id]);
	### Sum
	p.r = aggregate(coeff~., p, sum);
	return(reduce.pm(p.r));
}
diff.pm = function(p1, p2) {
	p2$coeff = - p2$coeff;
	return(add.pm(p1, p2));
}
sort.pm = function(p, sort.coeff=1, xn=NULL) {
	pP = p[, - which(names(p) == "coeff")];
	pow.tot = sapply(seq(nrow(p)), function(id) sum(pP[id, ]));
	pow.max = sapply(seq(nrow(p)), function(id) max(pP[id, ]));
	if(length(sort.coeff) == 1) {
		id = order(abs(p$coeff), pow.tot, pow.max);
	} else {
		coeff.df = data.frame(abs(p$coeff), -pow.tot, -pow.max);
		if( ! is.null(xn)) coeff.df$x = -pP[,xn];
		coeff.df = coeff.df[, sort.coeff];
		id = do.call(order, coeff.df)
	}
	return(p[id,])
}
eval.pm = function(p, x) {
	pP = p[, - which(names(p) == "coeff")];
	eval.p = function(id) {
		idx = which(pP[id,] != 0);
		prod(x[idx]^pP[id, idx], p$coeff[id]);
	}
	sum(sapply(seq(nrow(p)), eval.p))
}

#############
### Other ###

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

#########

### Print

# Print multi-variable Poly
print.monome = function(name, p) {
	v = p[,name];
	v.r = rep("", length(v));
	v.r[v > 1] = paste0(name, "^", v[v > 1]);
	v.r[v == 1] = name;
	return(v.r);
}
print.p = function(p, leading=1, order=TRUE, sort.order=TRUE) {
	### Var order
	if( ! is.numeric(leading)) leading = match(leading, names(p));
	if( ! is.na(leading)) {
		if(order) p = p[order(p[, leading], decreasing=sort.order), ];
		p = cbind(p[,-leading], p[,leading, drop=FALSE]);
	}
	###
	id.coeff = match("coeff", colnames(p));
	coeff = p[,id.coeff]; p = p[, - id.coeff];
	p.str = sapply(colnames(p), print.monome, p=p);
	# print(p.str)
	paste.nonempty = function(str, collapse="*") {
		str = str[nchar(str) > 0]
		paste(str, collapse=collapse)
	}
	if( ! is.null(dim(p.str))) p.str = apply(p.str, 1, paste.nonempty)
	else p.str = paste.nonempty(p.str);
	sign.str = ifelse(coeff > 0, " + ", " - ");
	sign.str[1] = if(coeff[1] > 0) "" else "- ";
	coeff.str = as.character(abs(coeff));
	hasCoeff = (abs(coeff) != 1);
	p.str[hasCoeff] = paste(coeff.str[hasCoeff], p.str[hasCoeff], sep = "*");
	return(paste(sign.str, p.str, sep="", collapse=""));
}
toCoeff = function(p, x="x") {
	idx = match(x, names(p));
	if(idx < 0) stop(paste0("No variable ", x));
	px = p[,x]; p = p[, - idx];
	str = tapply(seq(nrow(p)), px, function(nr) print.p(p[nr,], leading=NA))
	str[nchar(str) == 0] = "1";
	# missing powers
	x.all = seq(0, max(px));
	p.all = rep("0", length(x.all));
	p.all[1 + sort(unique(px))] = str;
	return(p.all)
}
print.coeff = function(p, x="x") {
	p = rev(toCoeff(p, x));
	sapply(p, function(p) cat(paste(p, ",\n", sep="")));
	invisible(p);
}

### Poly Calculations
perm2 = function(n, p=c(1,1)) {
	# all 2 permutations
	n1 = n - 1;
	len = n*n1/2;
	m = matrix(0, nrow=len, ncol=n);
	ioff = 0;
	for(i in seq(1, n1)) {
		m[seq(ioff+i, ioff+n1), i] = p[1];
		ioff = ioff + n1 - i;
	}
	nC = unlist(sapply(seq(2, n), function(id) seq(id, n)));
	for(i in seq(1, len)) {
		m[i, nC[i]] = p[2];
	}
	return(m)
}
perm3 = function(n, p=c(1,1,1)) {
	p2 = perm2(n, p=p[1:2]);
	if(p[1] != p[2]) p2 = rbind(p2, array(rev(p2), dim(p2)));
	if(any(p[3] == p[1:2])) {
		# TODO
	} else {
		p2 = t(p2);
		id = sapply(seq(ncol(p2)), function(id) which(p2[,id] == 0));
		putVar = function(pos=c(T,F)) {
			p21 = p2;
			id1 = id[rep(pos, length(id) %/% 2)];
			p21[seq(0, ncol(p2)-1)*n + id1] = p[3];
			return(p21);
		}
		p2 = cbind(putVar(c(T,F)), putVar(c(F,T)));
		p2 = t(p2);
	}
	return(p2);
}
prod.perm.poly = function(n, pow=c(1,1)) {
	# all 2 permutations
	m = perm2(n, p=pow);
	xn = paste0("x", seq(n));
	toPoly = function(mr) {
		id = which(mr != 0);
		pP = list(c(pow[1],0), c(0,pow[2]));
		names(pP) = xn[id];
		pP$coeff = c(1,1);
		return(pP);
	}
	if(n == 2) return(toPoly(m[1,]));
	pR = mult.pm(toPoly(m[1,]), toPoly(m[2,]));
	for(id in seq(3, nrow(m))) {
		pR = mult.pm(pR, toPoly(m[id,]));
	}
	return(pR)
}
perm.poly = function(n, p=c(1,1)) {
	m = as.data.frame(perm2(n, p=p))
	names(m) = paste0("x", seq(n));
	m$coeff = rep(1, nrow(m))
	return(m);
}
sym.poly = function(p, var="x") {
	n = length(p)
	l = list(p);
	l = rep(l, n);
	pP = 0; # TODO!
	names(pP) = paste0(var, seq(n));
	pP$coeff = 1;
	return(pP);
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

p.v = pow.pm(p, 3)
p.v

print.p(p.v[,c(2,3,4,1)])

################
### Permutations
n = 4
p = prod.perm.poly(n)
p = sort.pm(p);
p

print.p(p)

apply(perm3(4, p=c(3,2,1)), 1, sum)
table(duplicated(perm3(4, p=c(3,2,1))))

