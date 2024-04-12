

### Zeta Function

### D(zeta)
dz = function(s, eps = 1E-6) {
	(pracma::zeta(s+eps) - pracma::zeta(s)) / eps;
}

### Reflection Formula

s = 1/5;
dz(1 - s) + dz(s)*gamma(s)*sin(pi*(1-s)/2) * 2^(1-s) / pi^s # ==
pracma::zeta(1-s) * (log(2*pi) + pi/2 / tan(pi*(1-s)/2)) +
	- digamma(s)*gamma(s)*pracma::zeta(s) * sin(pi*(1-s)/2) * 2^(1-s) / pi^s;
