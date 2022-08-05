########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S4: C2-Hetero-Symmetric
### Helper Functions
###
### draft v.0.1b


# this file:
# source("Poly.System.S4.C2.Helper.R")


####################

### Helper Functions

source("Polynomials.Helper.R")


### Test
test.S4C2 = function(sol, n, R=NULL) {
	if(length(n) == 2) n = c(n, n[2]);
	#
	err1 = apply(sol^n[1], 1, sum);
	err4 = apply(sol, 1, prod);
	x1 = sol[,1]; x2 = sol[,2]; y1 = sol[,3]; y2 = sol[,4];
	err2 = x1^n[2]*y1^n[3] + x2^n[2]*y2^n[3];
	err3 = x1*x2*(y1 + y2) + y1*y2*(x1 + x2);
	err = cbind(err1, err2, err3, err4);
	if( ! is.null(R)) err = round0(err - rep(R, each=4));
	return(err)
}

### Test Variants
test.S4C2.Var = function(sol, n, type="x1x2", R=NULL) {
	if(length(n) == 2) n = c(n, n[2]);
	#
	type = pmatch(type, c("x1x2", "x1y2"));
	if(is.na(type)) stop("Invalid type!");
	len = length(n);
	nV  = if(len == 3) c(1,1)
		else if(len == 4) c(n[4], 1)
		else n[c(4, 5)];
	FUN = if(type == 1) function(x) { x[1]^nV[1] * x[2]^nV[2] + x[3]^nV[1] * x[4]^nV[2]; }
		else if(type == 2) function(x) { x[1]^nV[1] * x[4]^nV[2] + x[2]^nV[1] * x[3]^nV[2]; }
	#
	err1 = apply(sol^n[1], 1, sum);
	err4 = apply(sol, 1, FUN);
	x1 = sol[,1]; x2 = sol[,2]; y1 = sol[,3]; y2 = sol[,4];
	err2 = x1^n[2]*y1^n[3] + x2^n[2]*y2^n[3];
	err3 = x1*x2*(y1 + y2) + y1*y2*(x1 + x2);
	err = cbind(err1, err2, err3, err4);
	if( ! is.null(R)) err = round0(err - rep(R, each=4));
	return(err)
}
