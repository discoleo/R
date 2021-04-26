########################
###
### Leonard Mada
### [the one and only]
###
### Statistics: Moments
###
### draft v.0.1b


### Harmonic Moments & other Moments


####################

moments.h10 = function(x, pow=1) {
	if(pow != 1) x = x^pow;
	xm = sum(1/x);
	xm = length(x) / xm;
	if(pow != 1) xm = xm^(1/pow);
	return(xm)
}
moments.h20 = function(x, pow=1) {
	if(pow != 1) x = x^pow;
	xinv = 1/x;
	h.m = function(id) {
		v = xinv[id] * sum(tail(xinv, -id));
	}
	len = length(x) - 1;
	xm = sum(sapply(seq(len), h.m))
	xm = 1 / xm;
	### Normalization
	xm = xm * len * (len + 1) / 2;
	xm = if(pow != 1) xm^(1/(2*pow)) else (xm^(1/2));
	return(xm)
}
moments.h21 = function(x, pow=1, type=c("Mean", "Simple")) {
	type = match.arg(type);
	if(pow != 1) x = x^pow;
	xinv = 1/x;
	h.m = function(id) {
		xinv[id] * sum((x[id] + tail(x, -id)) * tail(xinv, -id));
	}
	len = length(x) - 1;
	xm = sum(sapply(seq(len), h.m))
	xm = 1 / xm;
	### Normalization
	# TODO: proper Normalization, e.g. when pow > 1?
	xm = xm * len * (len + 1);
	if(type == "Mean") {
		xm = xm * mean(x);
		xm = if(pow != 1) xm^(1/(2*pow)) else sqrt(xm);
	} else {
		xm = xm;
		xm = if(pow != 1) xm^(1/(pow)) else xm;
	}
	return(xm)
}


### Examples

n = 200
lambda = 4
x = rpois(n, lambda)
# exclude 0
x = x[x != 0]

### Moments

### H(1, 0) = n / sum(1/x[i])
### H(2, 0) = sqrt(n / sum(1/(x[i]*x[j])))
### H(2, 1) = n / sum((x[i]+x[j]) / (x[i]*x[j])))

### Q:
# - H(2, 1): What is the proper normalization?

moments.h10(x)
moments.h20(x)
moments.h21(x)
moments.h21(x, type="Simple")


#################

### Tests
x.test = rep(4, 20)
moments.h10(x.test)
moments.h20(x.test)
moments.h21(x.test)
moments.h21(x.test, type="Simple")

x.test = rep(7, 20)
moments.h10(x.test)
moments.h20(x.test)
moments.h21(x.test)
moments.h21(x.test, type="Simple")
