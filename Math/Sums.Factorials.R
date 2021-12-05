########################
###
### Leonard Mada
### [the one and only]
###
### Sums: Factorials
###
### draft v.0.1a


### Sums: Based on Factorials


########################

### Helper Functions

source("Polynomials.Helper.R")


sum0.fact = function(n, x=1) {
	i = seq(0, length(n)-1); sg = sign(x)^i; x = abs(x);
	sum(sg * x^n/factorial(n));
}
sum.fact = function(n, type, iter=36) {
	type = pmatch(type, c("Sum", "Harmonic", "First", "Second")); # TODO
	i = seq(0, iter);
	if(type == 1) {
		i = i[(i %% n) != 0];
		x = 1;
	} else if(type == 2) {
		i = i[(i %% n) != 0];
		x = -1;
	}
	r = sum0.fact(n=i, x=x);
	return(r);
}
sum.exp = function(n, type, x=1) {
	type = pmatch(type, c("Sum", "Harmonic", "First", "Second")); # TODO
	m = unity(n, all=FALSE);
	if(type == 1) {
		r = (exp(m*x) + exp(x/m))/2 - exp(x);
		r = r / (cos(2*pi/n) - 1);
	} else if(type == 2) {
		r = (exp(m*x) - exp(x/m))/2i;
		r = r / sin(2*pi/n);
	}
	return(r);
}


########################

###
n = 3
sum.exp(n, type="Sum")
sum.fact(n, type="Sum")
sum0.fact(c(1,2,4,5,7,8))

###
n = 3
sum.exp(n, type="Harm")
sum.fact(n, type="Harm")
sum0.fact(c(1,2,4,5,7,8), x=-1)

###
n = 3
(sum.exp(n, type="Sum") + sum.exp(n, type="Harm")) / 2
(sum.fact(n, type="Sum") + sum.fact(n, type="Harm")) / 2
sum0.fact(c(1,4,7,10))

