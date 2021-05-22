###################
###
### Rigidity Theory
###
### Leonard Mada
###
### draft v.0.1g


### Rigidity Theory
### Global Rigidity



###################

library(pracma)

### Helper functions

### generate 1D Framework
rpx1D = function(n, lim=c(1, 20), replace=TRUE) {
	sample(seq(lim[1], lim[2]), 10, replace=replace)
}
# connect the points
rpx1D.con = function(r, p=NULL) {
	# does NOT check if a framework is already possible;
	len = length(r);
	p = if(missing(p) || is.null(p)) NULL
		else if(length(p) >= 2) p[1:2]
		else c(p, 1-p);
	b = sample(c(1,-1), len, TRUE, prob=p);
	s = sum(b*r);
	if(s != 0) {
		b = c(b, -sign(s));
		r = c(r, abs(s));
	}
	return(list(p=r, b=b));
}

### Framework Realizations
### Problem:
# - massive redundancy!
# - the reductions simplify greatly the subproblems,
#   but create massive redundancy between the subproblems;
solve.fr1D = function(f) {
	l = test(f, 2);
	if(l$sum == 1) {
		l$SAT = FALSE;
		return(l);
	}
	l = test(f, 3);
	m = as.integer(names(l$tbl));
	len1 = l$tbl[m == 1]; len2 = l$tbl[m == 2];
	if(length(len1) == 0) len1 = 0;
	if(length(len2) == 0) len2 = 0;
	l$R = list();
	if(len1 == 0 && len2 == 0) {
		len = l$tbl[1];
		# few Elements vs Many
		if(len <= 3) {
			if(len == 1) {
				isSat = (f[1] == 0); id = 1;
			} else if(len == 2) {
				isSat = (f[1] == f[2]); id = c(1, -2);
			} else {
				iSat = c(f[1]+f[2]-f[3], f[1]-f[2]+f[3], f[1]-f[2]-f[3]);
				isSat.all = (iSat == 0); isSat = any(isSat.all);
				id = list(c(1,2,-3), c(1,-2,3), c(1,-2,-3));
				id = if(isSat) unlist(id[isSat.all]) else 0;
			}
			l$SAT = isSat; l$R = list(id=id);
		} else {
			l2 = solve.fr1D(f / 3);
			if(l2$SAT) {
				l$SAT = TRUE;
				l$R = c(l$R, l2$R);
			} else {
				l$SAT = FALSE;
			}
		}
		return(l);
	}
	# TODO:
	m3 = f %% 3;
	if(len1 == 0) {
		if(len2 == 1) {
			l$SAT = FALSE;
			return(l);
		}
		if(len2 %% 2 == 0) {
			f2 = combine.diff2(f[m3 == 2], add=f[m3 == 0], p=3);
			l2 = solve.fr1D(f2);
			if(l2$SAT) {
				l$R = c(l$R, l2$R);
			}
		}
		if(len2 %% 3 == 0) {
			f2 = combine.sum3(f[m3 == 2], add=f[m3 == 0], p=3);
			l2 = solve.fr1D(f2);
			if(l2$SAT) {
				l$R = c(l$R, l2$R);
			}
		}
	} else if(len2 == 0) {
	}
	return(l)
}
combine.diff2 = function(v, add, p=3) {
	# TODO: resolve massive redundancy!
	# - creating all combinations generates massive redundancy;
	id = seq(1, length(v), by=2);
	v2 = c(add, abs(v[id] - v[id + 1])) / p;
	return(v2);
}
combine.sum3 = function(v, add, p=3) {
	# TODO: resolve massive redundancy!
	# - creating all combinations generates massive redundancy;
	id = seq(1, length(v), by=3);
	v2 = c(add, abs(v[id] + v[id + 1] + v[id + 2])) / p;
	return(v2);
}
# exploratory test to reduce redundancy
latin6 = function() {
	m = matrix(c(
		1, 1, 1, 1, 1, 1, # s3 + s3
		1, 1, 1,-1,-1,-1, # s3 - s3
		1, 1,-1, 1,-1,-1,
		1, 1,-1,-1, 1,-1,
		1, 1,-1,-1,-1, 1,
		1,-1, 1, 1,-1,-1,
		1,-1, 1,-1, 1,-1,
		1,-1, 1,-1,-1, 1,
		-1,1, 1, 1,-1,-1,
		-1,1, 1,-1, 1,-1,
		-1,1, 1,-1,-1, 1
	), nrow=6)
	return(m);
}
# Test
sort(abs(apply(latin6() * sqrt(2:7), 2, sum)))

###  Test Sat
test = function(x, p) {
	list(sum=sum(x) %% p, tbl=table(x %% p), p=p)
}
simplify = function(x) {
	if(x$sum == 0) return(0);
	# works only with p = prime number;
	# Note: NOT fully correct!
	# 2 types of simplification: %% p or %% 2;
	tbl = x$tbl %% 2;
	# correction for: (%% p)-times;
	m = as.integer(names(tbl));
	# Note: should do the sum(tbl %% p == 0);
	# but even then it is only = 0 (mod p)!
	doCorrect = (tbl == 1) & (m != 0);
	cum.sum = sum((x$tbl[doCorrect] %/% x$p) * m[doCorrect]) + tbl[m == 0];
	tbl[m == 0] = cum.sum %% 2;
	tbl[doCorrect] = (x$tbl[doCorrect] %% x$p) %% 2;
	r = sum(tbl * m) %% x$p;
	# TODO: more simplifications possible;
	# TODO: all correct combinations;
	return(list(tbl=tbl, r=r));
}
simplify.p3 = function(x) {
	# full implementation of congruence (mod 3);
	# SAT: only for (mod 3)!
	if( ! is.list(x)) {
		l = test(x, 3);
	} else l = x;
	m = as.integer(names(l$tbl));
	simple.mod = function(l, mlog) {
		s2 = l$tbl[mlog];
		l$SAT = TRUE; l$type = "Simple";
		if(length(s2) == 0) return(l);
		if(s2 %% 2 == 0) {
			l$tbl[mlog] = 0;
		} else if(s2 %% 3 == 0) {
			l$tbl[mlog] = 0;
		} else {
			# else NO realization possible!
			l$SAT = FALSE;
		}
		return(l);
	} 
	if(all(m != 1)) {
		return(simple.mod(l, m == 2));
	} else if(all(m != 2)) {
		return(simple.mod(l, m == 1));
	}
	print("Non-simple")
	# both present: = 1 (mod 3) & = 2 (mod 3);
	s1 = l$tbl[m == 1]; s2 = l$tbl[m == 2];
	l$SAT = TRUE; l$tbl0 = l$tbl;
	if((s1 %% 3 == 0 || s1 %% 2 == 0) && (s2 %% 3 == 0 || s2 %% 2 == 0)) {
		l$tbl = l$tbl[m == 0];
		return(l);
	}
	if(s1 %% 2 == 0) {
		s2 = s2 %% 2; # == 1;
		s1 = -2;
		l$tbl[m == 1] = s1; l$tbl[m == 2] = s2;
		return(l);
	}
	if(s2 %% 2 == 0) {
		s1 = s1 %% 2; # == 1;
		s2 = -2;
		l$tbl[m == 1] = s1; l$tbl[m == 2] = s2;
		return(l);
	}
	s1 = s1 %% 2; # == 1;
	s2 = s2 %% 2; # == 1;
	l$tbl[m == 1] = s1; l$tbl[m == 2] = s2;
	return(l);
}

### Analyse
gcd.v = function(v, p) {
	gcd(v, p)
}
gcd.all = function(v) {
	len = length(v) - 1;
	d = sapply(seq(1, len), function(id) gcd(tail(v,-id), v[id]));
	d = unlist(d);
	m = diag(v);
	len = length(v);
	id = expand.grid(1:len, 1:len);
	m[id[,1] > id[,2]] = d;
	return(m);
}

#####################
#####################

#####################
### 1D Frameworks ###
#####################

n = 10;
f = rpx1D(n)


### Q: Is this a valid cyclic framework?
# (exists) b[i] = {-1, 1} such that sum(b*f) = 0;

test(f, 2)$tbl;

l = simplify.p3(f)
l

### Example:
# - at least 2 realizations (for any order);
f = c(1, 3, 8, 10, 13, 14, 17, 18, 18, 20)
l = simplify.p3(f)
l

### Solution:
# {1, 3, 13, 17}: {(1+/-3), (13+/-17)}, ...; # (mod 2)-variant;
# {1, 10, 13}, {8, 14, 17, 20}:
#   "3x1" + "2-2" => {(1+10+13), (8-14), (17-20)}, ...; # 3 variants
#   "3x2" + "1+2" => {(1-10), (13+8), (14+17+20)}, ...; # 4*3 = 12 variants;
### Case 1: "3x1" + "2-2"
# C.1.1: {3, 18, 18, 24, -6, -3} => scale by +/-1/3
#   => {1, 6, 6, 8, 2, 1}; # smaller problem;
#   =>  (1-1), {6,6,8,2} => NO realization!
#   =>  (1+1), {6,6,8,2} => realization: see below;
# C.1.2: {3, 18, 18, 24, -9, -6} => scale by +/-1/3
#   = > {1, 6, 6, 8, 3, 2}; # smaller problem;
#   =>  (1-3), {6,6,8,2} => 6+6-(8+2+3-1)
18 + 18 - (1+10+13) - (20-14 + (17-8) - 3); # = 0;
# C.1.3: {3, 18, 18, 24, -12, -3} => scale by +/-1/3
#   => {1, 6, 6, 8, 4, 1}; # smaller problem;
#   =>  (1-1), {6,6,8,4} => 6+6-8-4 + (1-1)
18 + 18 - (1+10+13) - (20-8) + 3 - (17-14) # (same as previous) OR
18 + 18 - (1+10+13) - (20-8) - 3 + (17-14); # = 0;
### Case 2: "3x2" + "1+2"
# C.2.1: {3, 18, 18, 14+17+20, 8+1, 10-13} => scale by +/-1/3
#   => {1, 6, 6, 17, 3, 1}; # smaller problem;
#   =>  (1-3), {6,6,8,2} => 17 - (6+6+1+3+1)
14 + 17 + 20 - (18+18+3+8+1+13-10); # = 0;


### TODO: solve recursive problem

###
# p = prime number
# if(sum(p) %% p == 1):
# NOT realizable (for odd n)!
# for n = even => more work;
p = 2
test(f, p)
p = 3
l = test(f, p)
simplify(l)
p = 5
l = test(f, p)
simplify(l)
p = 7
l = test(f, p)
simplify(l)


### TODO:
# - design full algorithm;

### generate a valid framework
rpx1D.con(f)


########################

### Complex Analysis
gcd.all(f)

