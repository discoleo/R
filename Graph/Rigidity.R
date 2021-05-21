###################
###
### Rigidity Theory
###
### Leonard Mada
###
### draft v.0.1e-ex2


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
test = function(v, p) {
	list(sum=sum(f) %% p, tbl=table(f %% p), p=p)
}
simplify = function(x) {
	if(x$sum == 0) return(0);
	# works only with p = prime number;
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
### TODO: solve recursive problem
### Case: "3x1" + "2-2"
# e.g. {3, 18, 18, 24, -6, -3} => scale by +/-1/3
#   => {1, 6, 6, 8, 2, 1}; # smaller problem;
#   =>  (1-1), {6,6,8,2} => NO realization!
# e.g. {3, 18, 18, 24, -9, -6} => scale by +/-1/3
#   = > {1, 6, 6, 8, 3, 2}; # smaller problem;
#   =>  (1-3), {6,6,8,2} => 6+6-(8+2+3-1)
18 + 18 - (1+10+13) - (20-14 + (17-8) - 3); # = 0;
# e.g. {3, 18, 18, 24, -12, -3} => scale by +/-1/3
#   => {1, 6, 6, 8, 4, 1}; # smaller problem;
#   =>  (1-1), {6,6,8,4} => 6+6-8-4 + (1-1)
18 + 18 - (1+10+13) - (20-8) + 3 - (17-14) # (same as previous) OR
18 + 18 - (1+10+13) - (20-8) - 3 + (17-14); # = 0;


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

