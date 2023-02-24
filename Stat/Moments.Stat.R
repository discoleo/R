########################
###
### Leonard Mada
### [the one and only]
###
### Statistics: Moments
###
### draft v.0.1e


### Harmonic Moments
### Geometric Moments
### Geometric Harmonic Moments
### & other Moments


### Moments

### Harmonic
### H(1, 0) = n / sum(1/x[i])
### H(2, 0) = sqrt(n / sum(1/(x[i]*x[j])))
### H(2, 1) = n / sum((x[i]+x[j]) / (x[i]*x[j])))


### Harmonic Dispersion

# HV(1, 0) = n / sum(1 / (HM - x[i]))
# HV(1, 1) = n / sum(x[i] / (HM - x[i]))
# where HM = harmonic mean; other types of mean also possible;


### Other

# IAS: Sihao Cheng - How to quantify random fields and textures
#   in astrophysics? (2023)
# https://www.youtube.com/watch?v=vDwEMC6UgLc
# - Scatter Statistics vs Power Spectrum vs CNN


####################

### Harmonic Moments
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
		# same as H10
		xm = xm;
		xm = if(pow != 1) xm^(1/(pow)) else xm;
	}
	return(xm)
}

### Geometric Moments
gmean.Gn0 = function(x, pow=NA) {
	xm = prod(x);
	len = length(x);
	if(len > 1) xm = xm^(1/len);
	return(xm)
}
gmean.Gnm10 = function(x, pow=1) {
	len = length(x);
	if(len < 2) return(x);
	if(pow != 1) x = x^pow; # TODO
	xm = prod(x);
	xs = sum(1/x);
	xm = xm * xs / len;
	if(len > 2) xm = xm^(1/(len-1));
	return(xm)
}

### Geometric Harmonic Moments
moments.GH20 = function(x, pow=1, useLog=TRUE) {
	if(useLog) {
		return(moments.GH20Log(x, pow=pow));
	}
	if(pow != 1) x = x^pow;
	xinv = 1/x;
	h.m = function(id) {
		v = xinv[id] + tail(xinv, -id);
	}
	len = length(x) - 1;
	xm = prod(unlist(lapply(seq(len), h.m)))
	xm = 1 / xm;
	### Normalization
	npow = pow * len * (len + 1) / 2;
	xm = xm^(1/npow);
	xm = xm * 2;
	return(xm)
}
moments.GH20Log = function(x, pow=1) {
	if(pow != 1) x = x^pow;
	xinv = 1/x;
	h.m = function(id) {
		v = xinv[id] + tail(xinv, -id);
	}
	len = length(x) - 1;
	xm = unlist(lapply(seq(len), h.m));
	xm = - sum(log(xm));
	### Normalization
	npow = pow * len * (len + 1) / 2;
	xm = xm / npow;
	xm = exp(xm) * 2;
	return(xm)
}


###################

### Examples

n = 200
lambda = 4
x = rpois(n, lambda)
# exclude 0
x = x[x != 0]

### Moments

### Q:
# - H(2, 1):
#   What is the proper normalization?
#   Simple H(2, 1) = H(1, 0);

### Harmonic
moments.h10(x)
moments.h20(x)
moments.h21(x)
moments.h21(x, type="Simple") # same as H10
### Geometric
gmean.Gn0(x)
gmean.Gnm10(x)
### Geometric Harmonic
moments.GH20(x, useLog=TRUE)


### Dispersion
mm = moments.h10(x)
n / sum(1 / (mm - x))
n / sum(x / (mm - x))

#
mm = mean(x)
n / sum(1 / (mm - x))
n / sum(x / (mm - x))


#################

### Tests
x.test = rep(4, 20)
moments.h10(x.test)
moments.h20(x.test)
moments.h21(x.test)
moments.h21(x.test, type="Simple")
gmean.Gn0(x.test)
gmean.Gnm10(x.test)
moments.GH20(x.test, useLog=TRUE)

x.test = rep(7, 20)
moments.h10(x.test)
moments.h20(x.test)
moments.h21(x.test)
moments.h21(x.test, type="Simple")
gmean.Gn0(x.test)
gmean.Gnm10(x.test)
moments.GH20(x.test, useLog=TRUE)
