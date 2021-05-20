###################
###
### Rigidity Theory
###
### Leonard Mada
###
### draft v.0.1d


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
	# 2 types of simplification: %% p or %% 2;
	tbl = x$tbl %% 2;
	# correction for: (%% p)-times;
	tbl[tbl == 1] = (x$tbl[tbl == 1] %% x$p) %% 2;
	r = sum(tbl * as.integer(names(tbl))) %% x$p;
	# TODO: more simplifications possible;
	# Note: should do the sum( tbl %% p == 0)!
	return(list(tbl=tbl, r=r));
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
gcd.all(f)

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

