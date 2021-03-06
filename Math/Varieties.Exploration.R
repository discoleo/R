########################
###
### Leonard Mada
### [the one and only]
###
### Varieties:
### Basic Exploration
###
### draft v.0.1b

### Introduction

### "A la recherche du mathematique perdue!"

# - a basic exploration of varieties;


####################
####################

### helper functions

library(polynom)
library(pracma)

# - the functions are in the file:
#   Polynomials.Helper.R

# - additional functions (e.g. specific solvers)
#   will be mentioned in the coresponding sections;


################################
################################


### "Triple" Elliptic Curve

# x^3 + b2*y^2 + b1*y = R
# y^3 + b2*z^2 + b1*z = R
# z^3 + b2*x^2 + b1*x = R

# - Solver: solve.Y2Y1.S3P3();
# - see file:
#   Poly.System.Hetero.Symmetric.S3.R;


variety = function(x, sVar, R=1, b=c(10,60), type="min", tol=1E-5) {
    solve.var = function(R, b2, b1) {
        sol = sort(solve.Y2Y1.S3P3(R, b=c(b2, b1), debug=F)[,1])
        sol = round0(sol, tol=tol);
        sol = Re(sol[Im(sol) == 0]);
		type = match(type, c("min", "max", "count", "min2", "max2"));
		if(is.na(type)) stop("Type NOT supported!")
		if(type == 1) {
			sol = if(length(sol) >= 1) sol[1] else NA;
		} else if(type == 2) {
			sol = if(length(sol) >= 1) sol[length(sol)] else NA;
		} else if(type == 3) {
			sol = length(sol);
		} else if(type == 4) {
			sol = if(length(sol) >= 3) sol[2] else NA;
		} else if(type == 5) {
			sol = if(length(sol) >= 4) sol[length(sol) - 1] else NA;
		}
        return(sol)
    }
	varType = match(sVar, c("R", "b1", "b2"))
	if(varType == 1) {
		sol = sapply(x, solve.var, b2=b[1], b1=b[2])
	} else if(varType == 2) {
		sol = sapply(x, solve.var, R=R, b2=b[1])
	} else if(varType == 3) {
		sol = sapply(x, solve.var, R=R, b1=b[2])
	}
    return(sol)
}
toRGB = function(col, scale=1) {
	col = col / scale;
	col[col < 0] = 1 + col[col < 0];
	rgb(col[1], col[2], col[3])
}

### Examples:

curve(variety(x, sVar="R", b=c(10,60)), from=-100, to=100, ylim=c(-10, 10))
curve(variety(x, sVar="R", b=c(10,60), type="max"), col="green", add=T)
curve(variety(x, sVar="R", b=c(10,60), type="count"), col="red", add=T)
curve(variety(x, sVar="R", b=c(10,60), type="min2"), col="orange", add=T)
curve(variety(x, sVar="R", b=c(10,60), type="max2"), col="orange", add=T)
text(-50, -2, "Varieties") # "sort of"


### Multiple
curve.v = function(b, type="max") curve(variety(x, sVar="R", b=c(10,b), type=type),
	col=toRGB(c(-b+20,60,0), 80), add=T);
by = 2;
curve(variety(x, sVar="R", b=c(10,60)), from=-100, to=100, ylim=c(-10, 10))
r = sapply(seq(30, 60, by=by), curve.v)
r = sapply(seq(30, 60, by=by), curve.v, type="min")
text(0, -1.25, "Varieties") # "sort of"


