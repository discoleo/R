########################
###
### Leonard Mada
### [the one and only]
###
### Integrals: Trigonometric
### Fresnel Integrals
###
### draft v.0.1c



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

###
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
p = - 3/4; # -n <= p <= 0;
# very slow!
# pracma::integral(\(x) x^p * sin(x^n), 0, 1000)
fresnel(n=n, p=p, k=1)
1/n * sin(pi*(p+1)/(2*n)) * gamma((p+1)/n)


### I( x^p * log(x) * sin(x^n) )
n = 3
p = - 3/4; # -n <= p <= 0;
# Upper = Inf: very slow!
pracma::integral(\(x) x^p * log(x) * sin(x^n), 0, 100)
gamma((p+1)/n)* (digamma((p+1)/n) * sin(pi*(p+1)/(2*n)) +
	pi/2 * cos(pi*(p+1)/(2*n))) / n^2


###
n = 3
p = - 2.5; # -n <= p <= 0;
# Upper = Inf: very slow!
pracma::integral(\(x) x^p * log(x) * sin(x^n), 0, 100)
gamma((p+1)/n)* (digamma((p+1)/n) * sin(pi*(p+1)/(2*n)) +
	pi/2 * cos(pi*(p+1)/(2*n))) / n^2

