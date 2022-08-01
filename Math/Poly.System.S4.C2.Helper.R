########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S4: C2-Hetero-Symmetric
### Helper Functions
###
### draft v.0.1a


####################

### Helper Functions

source("Polynomials.Helper.R")


# Test
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
