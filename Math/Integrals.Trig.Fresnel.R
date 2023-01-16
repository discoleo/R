########################
###
### Leonard Mada
### [the one and only]
###
### Integrals: Trigonometric
### Fresnel Integrals
###
### draft v.0.1a



# much faster than naive pracma::integral;
fresnel = function(n, iter=10000, k=1) {
	id  = seq(1, iter);
	id0 = c(0, id);
	ninv = 1/n;
	limI = ((2*pi)*id0)^ninv;
	int = sapply(id, function(iStart) {
		lim = limI[c(iStart, iStart + 1)];
		pracma::integral(function(x) sin(k*x^n), lim[1], lim[2]);
	})
	int = sum(int);
	return(int);
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

