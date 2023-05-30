

### Analytic Continuation

### Partial Zeta-Functions
# => Digamma & Polygamma Functions
### Other Partial Functions


# For an Introduction, see link to Digamma Function:
# - Lines That Connect: Extending the Harmonic Numbers to the Reals
#   https://www.youtube.com/watch?v=9p_U_o1pMKo



####################

### Helper Functions

### Simple ###

### Based on H(n)

### H(3*n) - H(n)
# - accuracy/convergence problems possible with some of the formulas;
#   (I may have had a bug in the initial code as well)
H3n = function(x, N = 30) {
	n = seq(0, N);
	sapply(x, \(x) {
		m = n + x;
		sum( (9*n+5) / ((3*n+1)*(3*n+2)*(3*n+3)) ) +
			- sum( (9*m+5) / ((3*m+1)*(3*m+2)*(3*m+3)));
	})
}
H3n0 = function(n) {
	sum(1/seq(3*n)) - sum(1/seq(n));
}
H3n0.alt = function(n) {
	n = seq(0, n-1);
	sum( (9*n+5) / ((3*n+1)*(3*n+2)*(3*n+3)) );
}


### Based on the Zeta-function
# - Zeta Partial-Sums;

### sum( 1 / n^p )
zetaSum = function(x, p=2, N = 50) {
	id = seq(1, N);
	sapply(x, \(x) sum(1/id^p) - sum(1/(x+id)^p));
}
zetaSum0 = function(n, p=2) {
	id = seq(1, n);
	sum(1/id^p);
}


### Quasi-Zeta

### sum( 1 / n^n )
zetaQSum = function(x, N = 30) {
	id = seq(1, N);
	sapply(x, \(x) {
		m = id + x;
		sum(1/id^id) - sum(1/m^m);
	})
}
zetaQSum0 = function(n) {
	id = seq(1, n);
	sum(1/id^id);
}

### sum( (-1)^n / n^n )
zetaQnegSum = function(x, N = 30) {
	id = seq(1, N);
	u  = (-1)^id;
	sapply(x, \(x) {
		m = id + x;
		b = sin(pi*m - pi/2);
		- sum(u/id^id) - sum(b/m^m);
	})
}
# non-Smooth
zetaQnegSum.nonSmooth = function(x, N = 30) {
	id = seq(1, N);
	u  = (-1)^id;
	sapply(x, \(x) {
		m = id + x;
		b = if(floor(x) %% 2 == 0) 1 else -1;
		- sum(u/id^id) + b * sum(u/m^m);
	})
}
zetaQnegSum0 = function(n) {
	id = seq(1, n);
	- sum((-1)^id / id^id);
}

### Alternating Zeta

### sum( (-1)^n / n^p )
zetaAltSum = function(x, p=1/2, N = 20000) {
	id = seq(1, N);
	u  = (-1)^id;
	sapply(x, \(x) {
		m = id + x;
		b = sin(pi*m - pi/2);
		- sum(u/id^p) - sum(b/m^p);
	})
}
zetaAltSum0 = function(n, p=1/2) {
	id = seq(1, n);
	- sum((-1)^id / id^p);
}


### Exponentials ###

### sum( n / exp(n) )
invExpSum = function(x, N = 30) {
	id = seq(0, N);
	sapply(x, \(x) sum(id/exp(id)) - sum((x+1+id)/exp(x+1+id)));
}
invExpSum0 = function(n) {
	id = seq(0, n);
	sum(id/exp(id));
}


### sum( 1 / (n*exp(n) + 1) )
invXExpSum = function(x, k=1, N = 30) {
	id = seq(0, N);
	sapply(x, function(x) {
		sum(1/(id*exp(id) + k)) - sum(1/((x+1+id)*exp(x+1+id) + k))
	});
}
invXExpSum0 = function(n, k=1) {
	id = seq(0, n);
	sum(1 / (id*exp(id) + k) );
}

#####################
#####################

################
### Examples ###
################

### Simple

### H(3*n) - H(n)

px = c(1,5,10,15);
curve(H3n(x), xlim=c(0, 30), ylim=c(0, 1.2))
points(px, sapply(px, H3n0), col="red")
abline(h = log(3), col="red")
# for derivation of log(3), see file:
# Sums.Trig.R;

###
H3n(1/2, N=10000)
2/3


##################
##################

### Zeta-based

### sum( 1 / n^p )
curve(zetaSum(x), xlim=c(0,10))
tmp = sapply(seq(1, 8), \(x) points(x, zetaSum0(x), col="red"))
#
p = 3
curve(zetaSum(x, p=p), add=T, col="blue")
tmp = sapply(seq(1, 8), \(x) points(x, zetaSum0(x, p=p), col="red"))


### Zeta(2)
# Q: Can we compute a closed form for zetaSum(3/2)?
zetaSum(1/2, N=200)
4 - pi^2/3
zetaSum(3/2, N=200)
4 - pi^2/3 + 4/9
zetaSum(5/2, N=200)
4 - pi^2/3 + 4/9 + 4/25


### Zeta(3)
zetaSum(1/2, p=3)
8 - 6*pracma::zeta(3)

###
zetaSum(1/3, p=3, N=1000)
27 - 12 * pracma::zeta(3) - pi^3 * cos(pi/3) / sin(pi/3)^3 / 2
#
zetaSum(2/3, p=3, N=1000)
27/8 - 12 * pracma::zeta(3) + pi^3 * cos(pi/3) / sin(pi/3)^3 / 2


# Derivation:
pracma::psi(2, 2/3) - pracma::psi(2, 1/3)
2 * pi^3 * cos(pi/3) / sin(pi/3)^3

#
pracma::psi(2, 1/3) / 2 + 27 + pracma::zeta(3)
zetaSum(1/3, p=3)
#
pracma::psi(2, 2/3) / 2 + 27/8 + pracma::zeta(3)
zetaSum(2/3, p=3)


### Zeta(3): zetaSum(1/4) & zetaSum(3/4)
zetaSum(1/4, p=3, N=1000)
4^3 - 27 * pracma::zeta(3) - pi^3/2 * cos(pi/4) / sin(pi/4)^3
zetaSum(3/4, p=3, N=1000)
(4/3)^3 - 27 * pracma::zeta(3) + pi^3/2 * cos(pi/4) / sin(pi/4)^3

# Derivation:
pracma::psi(2, 3/4) - pracma::psi(2, 1/4)
2 * pi^3 * cos(pi/4) / sin(pi/4)^3

#
zetaSum(1/4, p=3, N=1000) + zetaSum(3/4, p=3, N=1000)
64 + 64/27 + (2 - 64 + 64/8) * pracma::zeta(3)
#
pracma::psi(2, 3/4) / 2 + (4/3)^3 + pracma::zeta(3)
zetaSum(3/4, p=3, N=1000)


###
n = 5; # n = sqrt(11);
pracma::psi(2, 1 - 1/n) - pracma::psi(2, 1/n)
2 * pi^3 / sin(pi/n)^2 / tan(pi/n)


##############
##############

### Quasi-Zeta
# sum( 1 / n^n )

curve(zetaQSum(x), xlim=c(0, 6))
tmp = sapply(seq(1, 6), \(x) points(x, zetaQSum0(x), col="red"))

###
curve(zetaQSum(x), xlim=c(1.95, 7))
tmp = sapply(seq(2, 6), \(x) points(x, zetaQSum0(x), col="red"))

###
zetaQSum(1/2)
# TODO: ???


###############

### Quasi-Zeta
### Alternating

# sum( (-1)^n / n^n )

# TODO: Clarify if the use of sin() is the best option?
curve(zetaQnegSum(x), xlim=c(0, 6))
tmp = sapply(seq(1, 6), \(x) points(x, zetaQnegSum0(x), col="red"))

###
zetaQnegSum(20)
zetaQnegSum0(20)

###
zetaQnegSum(1/2)


### Bernoulli Integral
# Dr. Trefor Bazett: The Bernoulli Integral is ridiculous
# https://www.youtube.com/watch?v=PxyK_Tsnz10
# Dr. Peyam: Integral x^x from 0 to 1
# https://www.youtube.com/watch?v=A54_QPXdkU0


integrate(\(x) x^x, 0, 1)
zetaQnegSum(50)

###
integrate(\(x) x^(- x), 0, 1)
zetaQSum(50)


####################

### Alternating Zeta
# sum( (-1)^n / n^p )

# TODO: Clarify if the use of sin() is the best option?
# Note: large N needed for convergence of odd terms;
curve(zetaAltSum(x), xlim=c(0, 10))
tmp = sapply(seq(1, 8), \(x) points(x, zetaAltSum0(x), col="red"))


###
curve(zetaAltSum(x, p=3/2, N=5000), xlim=c(0, 10))
tmp = sapply(seq(1, 8), \(x) points(x, zetaAltSum0(x, p=3/2), col="red"))


################
################

### Exponentials

### sum( n / exp(n) )
curve(invExpSum(x), xlim=c(0,10))
tmp = sapply(seq(0, 5), \(x) points(x, invExpSum0(x), col="red"))


#####################

### sum( 1 / (n*exp(n) + 1) )
curve(invXExpSum(x), xlim=c(0,10))
tmp = sapply(seq(0, 6), \(x) points(x, invXExpSum0(x), col="red"))
#
k = 1 + 1/16
curve(invXExpSum(x, k=k), add = TRUE)
tmp = sapply(seq(0, 6), \(x) points(x, invXExpSum0(x, k=k), col="red"))


# NOT identical: still some diff?
invXExpSum0(50)
# exp(1)/2 # close, but ???
integrate(\(x) 1/(x*exp(x) + 1), 0, Inf)

