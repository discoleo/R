########################
###
### Leonard Mada
### [the one and only]
###
### Integrals: Trigonometric
### Fresnel Integrals
###
### draft v.0.2d



# TODO:
# - evaluate if any of the concepts are useful:
#   Giray Ökten: Number sequences for simulation - Lecture 2
#   https://www.youtube.com/watch?v=3YFVV5YqLsI


# I( sin(k*x^n) ) on [0, Inf]
# much faster than naive pracma::integral;
fresnel = function(n, iter=10000, k=1, p=0, doWeights = TRUE) {
	id  = seq(1, iter);
	id0 = c(0, id);
	ninv = 1/n;
	# Note: 1*pi vs 2*pi does NOT help!
	limI = ((2*pi)*id0)^ninv;
	int = sapply(id, function(iStart) {
		lim = limI[c(iStart, iStart + 1)];
		if(p == 0) {
			pracma::integral(function(x) sin(k*x^n), lim[1], lim[2]);
		} else {
			pracma::integral(function(x) x^p * sin(k*x^n), lim[1], lim[2]);
		}
	})
	if(doWeights && p == 0) {
		# w = (1 - id/iter^2);
		# problematic for n > 2, but better for n = 2;
		w = (1 - rev(1/limI[-1])/(n^3 - 13/2));
		# w = (1 - 1/rev(limI[-1]*iter));
		w = w / mean(w);
		# Gauss: does not work (on 2*pi);
		# xn = (head(limI, -1) + tail(limI, -1)) / 2;
		# w = 2 / (n^2 * xn^(2*n-2) * (1 - xn^2) * (1 - int^2));
		int = sum(w*int);
		
	} else int = sum(int);
	return(int);
}
# much faster, but far LESS robust;
fresnel2 = function(n, iter=10000, subdiv=16, k=1) {
	id = seq(iter);
	ninv = 1/n;
	limX = ((2*pi)*id)^ninv;
	limI = c(0, limX);
	x = sapply(id, function(id) {
		x = seq(limI[id], limI[id + 1], length.out=subdiv);
		# return(x[-c(1, subdiv)]);
		return(x);
	})
	int = sin(k*x^n);
	int[seq(1, iter*subdiv, by=subdiv)] = 0;
	int[seq(subdiv, iter*subdiv, by=subdiv)] = 0;
	dfx = diff(c(0, x));
	int = sum(dfx * int);
	return(int)
}

# may take long!
test.fresnel = function(n=2, iter=c(1E4, 2E4), doWeights=TRUE, plot=TRUE, length=50) {
	it  = seq(iter[1], iter[2], length.out=length);
	cat("Started: ");
	val = sapply(seq(length), function(id) {
		if(id %% 10 == 0) cat(id, ", ", sep="");
		tmpIter = it[id];
		fresnel(n, iter=tmpIter, doWeights=doWeights);
	})
	cat(length, "\n", sep="");
	if(plot) {
		val0 = 1/n * sin(pi/(2*n)) * gamma(1/n);
		val1 = fresnel(n, iter=iter[1], doWeights=FALSE);
		ylim = c(min(val, val1, val0), max(val0, val1, val));
		plot(it, val, ylim=ylim);
		abline(h = val0, col="green");
		abline(h = val1, col="red");
	}
	invisible(val);
}

# Note:
# - numerically very problematic!

### I( sin(k*x^n) ) on [0, Inf]
n = 2
fresnel(n, iter=20000)
1/n * sin(pi/(2*n)) * gamma(1/n)

###
n = 3
fresnel(n)
1/n * sin(pi/(2*n)) * gamma(1/n)

###
n = 4
fresnel(n)
1/n * sin(pi/(2*n)) * gamma(1/n)

### n -> Inf
# x = y^(1/p) & p -> 0;
# upper = Unf;
integrate(\(x) sin(x) / x, 0, 10000, subdivisions = 1025)
n = 2000
sin(pi/(2*n)) * gamma(1/n)
pi / 2


### Weighted:

# Types of Corrections:
# - Weighting/"Overfitting":
#   lower weight for later [0, 2*pi]-intervals,
#   as late intervals become oscillatory and converging to 0;
# - Residual: add some residual corresponding to the tail;

# takes ~ 5 min!
test.fresnel(2)
# takes ~ 5 min!
test.fresnel(3)


###############
### with k <> 1

###
n = 2
k = 3
# presumably correct:
fresnel(n, iter=20000, k=k)
1/n * sin(pi/(2*n)) * gamma(1/n) / k^(1/n)

###
n = 3
k = 3
# presumably correct:
fresnel(n, k=k)
1/n * sin(pi/(2*n)) * gamma(1/n) / k^(1/n)


####################
####################

### I( x^p * sin(x^n) )
n = 3
p = - 3/4; # - (n+1) < p <= 0;
# very slow!
# pracma::integral(\(x) x^p * sin(x^n), 0, 1000)
fresnel(n=n, p=p, k=1)
1/n * sin(pi*(p+1)/(2*n)) * gamma((p+1)/n)

###
n = 3
# - (n+1) < p <= 0;
# but numerically problematic for p < - (n + 0.5);
p = 0.5 - (n+1);
fresnel(n=n, p=p, k=1)
1/n * sin(pi*(p+1)/(2*n)) * gamma((p+1)/n)


### Log

### I( x^p * log(x) * sin(x^n) )
n = 3
p = - 3/4; # - (n+1) < p < 0;
# Upper = Inf: very slow!
pracma::integral(\(x) x^p * log(x) * sin(x^n), 0, 100)
gamma((p+1)/n) * (digamma((p+1)/n) * sin(pi*(p+1)/(2*n)) +
	pi/2 * cos(pi*(p+1)/(2*n))) / n^2


###
n = 3
p = - 2.5; # - (n+1) < p < 0;
# Upper = Inf: very slow!
pracma::integral(\(x) x^p * log(x) * sin(x^n), 0, 100)
gamma((p+1)/n)* (digamma((p+1)/n) * sin(pi*(p+1)/(2*n)) +
	pi/2 * cos(pi*(p+1)/(2*n))) / n^2


###
n = 3
p = - 3.25; # - (n+1) < p < 0;
# Upper = Inf: very slow!
pracma::integral(\(x) x^p * log(x) * sin(x^n), 0, 100)
gamma((p+1)/n)* (digamma((p+1)/n) * sin(pi*(p+1)/(2*n)) +
	pi/2 * cos(pi*(p+1)/(2*n))) / n^2


#####################

### Higher Power: sin

### I( sin(x^n)^2 / x^p )
# 1 < p < 2*n + 1;
n = 2
p = sqrt(5);
# Upper = Inf: very slow!
pracma::integral(\(x) sin(x^n)^2 / x^p, 0, 6000)
sin(pi*(1/2 - (p-1)/(2*n))) * gamma(1 - (p-1)/n) * 2^((p-1)/n - 1) / (p-1)


### works: 0 < n <= 1
# (as it depends on x^p for convergence)
# Note: n = 0 implies 1 < p < 1 => Contradiction;
n = 1/2
p = 1.4
pracma::integral(\(x) sin(x^n)^2 / x^p, 0, 1000000)
sin(pi*(1/2 - (p-1)/(2*n))) * gamma(1 - (p-1)/n) * 2^((p-1)/n - 1) / (p-1)

###
n = 2
p = 2*n + 0.45;
# Upper = Inf: very slow!
pracma::integral(\(x) sin(x^n)^2 / x^p, 0, 6000)
sin(pi*(1/2 - (p-1)/(2*n))) * gamma(1 - (p-1)/n) * 2^((p-1)/n - 1) / (p-1)


######################
######################

### Composite Fresnel-Type

###
# Maths 505: A RIDICULOUSLY AWESOME INTEGRAL!!!! int 0 to infty (sin(x^2+1/x^2))^3
# https://www.youtube.com/watch?v=CJOZoV7S2l4

### I( sin(x^2 + 1/x^2)^3 )
# upper = Inf: numerical issues;
pracma::integral(\(x) sin(x^2 + 1/x^2)^3, 0, 512*pi)
1/8 * sin(pi/4) * gamma(1/2) * (3*sin(2) + 3*cos(2) - (sin(6) + cos(6))/sqrt(3))


###
# upper = Inf: + workarounds for numerical issues;
# pracma::integral(\(x) sin(x^4 + 1/x^4)^3, 0, 256*pi)
pracma::integral(\(x) 1/4 * sin(x + 1/x)^3 / x^(3/4), 1, 65536*pi) +
	pracma::integral(\(x) 1/4 * sin(x + 1/x)^3 / x^(3/4), pi/(2*4096), 1)
# pracma::integral(\(x) sin(x^4 + 4*x^2 + 2)^3, 0, 256*pi)
pracma::integral(\(x) 1/2 * sin(x^2 + 4*x + 2)^3 / x^(1/2), 0, 1024*pi)
# TODO: sin(x)^3 = (3*sin(x) - sin(3*x))/4

# vs:
# pracma::integral(\(x) sin(x^4 - 4*x^2 - 4/x + 1/x^4)^3, 0, 256*pi)


###
# upper = Inf: numerical issues;
# pracma::integral(\(x) sin(x^4 + 4*x^2 + 2), 0, 256*pi)
pracma::integral(\(x) 1/2 * sin(x^2 + 4*x + 2) / x^(1/2), 0, 1024*pi)
pracma::integral(\(x) sin((x^2 + 2)^2)*cos(2) - cos((x^2 + 2)^2)*sin(2), 0, 256*pi)
pracma::integral(\(x) sqrt(2) * sin(4*(x^2 + 1)^2)*cos(2) +
	- sqrt(2) * cos(4*(x^2 + 1)^2)*sin(2), 0, 256*pi)
# TODO: ???


##################

### I( sin(a^2*x^2 + b^2/x^2) )
# Maths 505: A surprising integral result
# https://www.youtube.com/watch?v=P4ROv4iiK-4
# - Substitution: x = b/a / y => new I ;
#   then a * I + b * new I: t = a*x - b/x;

a = 2; b = 3;
# Note: on [0, Inf] but numerically extremely problematic!
integrate(\(x) sin(a^2*x^2 + b^2/x^2), 1/(4*pi), 4*pi, subdivisions = 529)
integrate(\(x) sin(x^2 + 2*a*b) / a, 1/(4*pi), 4*pi)
(sin(2*a*b) +  cos(2*a*b)) * sin(pi/4) * gamma(1/2) / (2*a);


# library(Rmpfr)
integrate(\(x) {
	x = mpfr(x, 240);
	y = sin(a^2*x^2 + b^2/x^2);
	as.numeric(y);
	}, 1/(32*sqrt(pi)), 2*sqrt(pi), subdivisions = 5025)$value +
integrate(\(x) {
	x = mpfr(x, 240);
	y = sin(a^2*x^2 + b^2/x^2);
	as.numeric(y);
	}, 2*sqrt(pi), 32*sqrt(pi), subdivisions = 5025)$value;
	
