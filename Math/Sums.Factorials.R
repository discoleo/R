########################
###
### Leonard Mada
### [the one and only]
###
### Sums: Factorials
###
### draft v.0.1c


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

### x^j / Gamma(j + k)
# k = Offset of factorial;
sum.GammaExt = function(k, x=1, iter=20) {
	i = seq(0, iter);
	x = if(x == 1) x else x^i;
	i = i + k + 1;
	sum(x/gamma(i))
}


########################

### Basic Factorials

### By = 3
n = 3;
m = exp(seq(0, n-1) * 2i*pi/n);
id = seq(0, 10);

### Sum ( 1 / (3*j)! )
sum(exp(m)) / n;
sum( 1 / factorial(n*id) )

### Sum ( 1 / (3*j + 1)! )
sum(exp(m) * m[c(1,3,2)]) / n;
sum( 1 / factorial(n*id + 1) )

### Sum ( 1 / (3*j + 2)! )
sum(exp(m) * m[c(1,2,3)]) / n;
sum( 1 / factorial(n*id + 2) )


### Parameter: x != 1

### Sum ( x^j / (3*j)! )
x = 1/2; # x = 2; # x = 2^(1/3);
sum(exp(x*m)) / n;
sum( x^(n*id) / factorial(n*id) )


### Case: x = log(2)
# Very special Log;
x = log(2); # x = log(2)^(1/3);
sum(exp(x*m)) / n;
sum( x^(n*id) / factorial(n*id) )
#
sum(exp(x*m) * m[c(1,3,2)]) / (n*x);
sum( x^(n*id) / factorial(n*id + 1) )


### Case: x^3 = x + 1
x = pracma::roots(c(1,0,-1,-1));
x1 = x[1];
sum(exp(x1*m)) / n;
sum( x1^(n*id) / factorial(n*id) )
#
sum(exp(x1*m) - exp(m)) / n;
sum( (x1^(n*id) - 1) / factorial(n*id) );
sum( ((x1+1)^id - 1) / factorial(n*id) );

###
x = pracma::roots(c(1,0,-1,-1));
sum(exp(x[2]*m) - exp(x[3]*m)) / (2i*n);
sum( (x[2]^(n*id) - x[3]^(n*id)) / factorial(n*id) ) / 2i;


#############
### By = 5
n = 5;
m = exp(seq(0, n-1) * 2i*pi/n);
id = seq(0, 10);

### Sum ( 1 / (5*j)! )
sum(exp(m)) / n;
sum( 1 / factorial(n*id) )

### Sum ( 1 / (5*j + 1)! )
sum(exp(m) * m[c(1,5,4,3,2)]) / n;
sum( 1 / factorial(n*id + 1) )

### Sum ( 1 / (5*j + 2)! )
sum(exp(m) * m[c(1,4,2,5,3)]) / n;
sum( 1 / factorial(n*id + 2) )

### Sum ( 1 / (5*j + 3)! )
sum(exp(m) * m[c(1,3,5,2,4)]) / n;
sum( 1 / factorial(n*id + 3) )

### Sum ( 1 / (5*j + 4)! )
sum(exp(m) * m[c(1,2,3,4,5)]) / n;
sum( 1 / factorial(n*id + 4) )


##########

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


######################
######################

### Extensions:
### over Gamma

###
n = - 1/2
sum.GammaExt(n)


###
n = 1/2
sum.GammaExt(n)

###
n = 1/3
sum.GammaExt(n)

###
n = 2/3
sum.GammaExt(n)

# TODO: explore these extensions;

