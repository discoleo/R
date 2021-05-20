###################
###
### Rigidity Theory
###
### Leonard Mada
###
### draft v.0.1b-test


### Rigidity Theory
### Global Rigidity



###################

library(pracma)

### Helper functions

### generate 1D Framework
rpx1D = function(n, lim=c(1, 20), replace=TRUE) {
	sample(seq(lim[1], lim[2]), 10, replace=replace)
}
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
	list(sum=sum(f) %% p, tbl=table(f %% p))
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
p = 3
test(f, p)
p = 5
test(f, p)
p = 7
test(f, p)


### TODO:
# - design full algorithm;

