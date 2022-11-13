########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S5: Hetero-Symmetric
### Useful Formulas
###
### draft v.0.1a


### Derivations
# - Basic derivations;
# - Numerical approaches: basic intuition;


####################

### Helper Functions

source("Polynomials.Helper.R")
source("Polynomials.Helper.EP.R")


# Solve for all initial tuples in x0
solve.all = function(FUN, x0, ..., debug=TRUE) {
	if(is.null(dim(x0))) x0 = matrix(x0, nrow=1);
	nr = nrow(x0);
	x.all = array(0, c(ncol(x0), 0));
	#
	for(id in seq(nr)) {
		xi = rbind(Re(x0[id,]), Im(x0[id,]));
		xx = multiroot(solve.S5HtMixed.Num, start=xi, ...);
		#
		x = xx$root;
		x = matrix(x, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
		if(debug) {
			cat(paste0("ID = ", id, "; Iter = ", xx$iter, "; Prec = ", xx$estim.precis));
			cat("\n   "); cat.sol(xx);
		}
		x.all = cbind(x.all, x);
	}
	x.all = t(x.all);
	colnames(x.all) = paste0("x", seq(ncol(x0)));
	rownames(x.all) = NULL;
	return(x.all);
}
# Solve using given Path:
solve.path = function(FUN, x0, path, debug=TRUE) {
	n = length(path);
	for(i in seq(n)) {
		if(debug) cat(paste0("Step: ", i, "\n"));
		x0 = solve.all(FUN, x0=x0, R=path[[i]], debug=debug);
	}
	return(x0)
}

### Print/Format
cat.sol = function(x, digits=4, sep="\n") {
	if(inherits(x, "list")) {
		x = x$root;
		x = matrix(x, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
	}
	cat("c(");
	cat(paste0(round(x, digits=digits), collapse=", ")); cat(")"); cat(sep);
}
cat.sol.m = function(x, digits=4, sep=",\n") {
	if( ! inherits(x, "matrix")) {
		return(cat.sol(x, digits=digits));
	}
	apply(x, 1, cat.sol, digits=digits, sep=sep);
	invisible();
}

#######################
#######################

### Base-System

### Numerical Solution:

# library(rootSolve)

solve.S5HtMixed.Num = function(x, R=c(0,1,0,0,1)) {
	x = matrix(x, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
	x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4]; x5 = x[5];
	s1 = x1 + x3; s2 = x2 + x4; S = s1 + s2 + x5;
	p1 = x1 * x3; p2 = x2 * x4; E5 = p1 * p2 * x5;
	ps = s1 * s2; sp = p1 + p2;
	# E2 = x1*(S - x1) + x2*(x3 + x4 + x5) + x3*(x4 + x5) + x4*x5;
	# E2 = sp + ps + x5*(S - x5);
	E3 = p1*s2 + p2*s1 + x5*(sp + ps);
	E4 = p1*p2 + x5*(p1*s2 + p2*s1);

	### E2:
	E11a = x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x1;
	#
	y = c(S, E11a, E3, E4, E5) - R;
	y = rbind(Re(y), Im(y));
	return(y);
}

# Note:
# - the double permutation (x1, x3) & (x4, x5) like c(x3, x2, x1, x5, x4),
#   and equivalent permutations (e.g. (x1, x2) & (x3, x5)),
#   are also a solution;

### Set 1:
x0 = c(-0.70449+0.64i, 0.8913, -0.70449-0.64i, 0.2589 + 1.08i, 0.2589 - 1.08i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
print(poly.calc(x)[4], 12)
x.all = c(x)


### Set 2:
x0 = c(-1.70449-0.64i, -0.8913, -1.70449+0.64i, -0.2589 + 1.08i, -0.2589 - 1.08i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
print(poly.calc(x)[4], 12)
x.all = c(x.all, x)


### Set 3:
x0 = c(0.05-0.94i, -2.8913, 0.05+0.04i, -0.2589 + 1.08i, -0.2589 - 1.08i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
print(poly.calc(x)[4], 12)
x.all = c(x.all, x)


### Set 4:
x0 = c(-1.274+0.729i, 1.23-0.47i, 0.58-0.67i, 0.178+0.725i, -0.71 - 0.31i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
# print(poly.calc(x)[4], 12)
x.all = c(x.all, x)


### Set 5:
x0 = c(-1.274-0.729i, 1.23+0.47i, 0.58+0.67i, 0.178-0.725i, -0.71 + 0.31i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
# print(poly.calc(x)[4], 12)
x.all = c(x.all, x)


### Set 6:
x0 = c(-0.8108+1.5014i, 0.7763-1.6039i, 0.6605-0.1778i, -0.5008-0.4342i, -0.1252+0.7145i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
# print(poly.calc(x)[4], 12)
x.all = c(x.all, x)


### Set 7:
x0 = c(-0.8108-1.5014i, 0.7763+1.6039i, 0.6605+0.1778i, -0.5008+0.4342i, -0.1252-0.7145i);
x0 = rbind(Re(x0), Im(x0))
xx = multiroot(solve.S5HtMixed.Num, start=x0)

x = matrix(xx$root, nr=2); xc = x[2,]; x = x[1,] + 1i * xc;
# print(poly.calc(x)[4], 12)
x.all = c(x.all, x)
x.all = matrix(x.all, nc=5, byrow=T)

round0(poly.calc(x.all)) * 27
poly.calc(apply(x.all, 1, function(x) sum(x * x[c(3,4,5,1,2)]))) * 27

###################

### Case 3:
# S = 1; E11a = 0;

### Cases 4 - 6:
# S =  1; E11a =  1;
# S =  1; E11a = -1;
# S = -1; E11a =  1;

###
x0All = list(
	# R2 = c(0,1,0,0,2)
	V01 = rbind(
	c(-0.8392+0.7287i, 1.051, -0.8392-0.7287i, 0.3137+1.2009i, 0.3137-1.2009i),
	c(-0.42-0.6586i, 0.7387+0i, -0.42+0.6586i, 0.0507+2.1058i, 0.0507-2.1058i),
	c(0.392+0.8686i, 1.4881+0i, 0.392-0.8686i, -1.136-0.4352i, -1.136+0.4352i),
	c(-0.1185+0.8347i, -0.591-0.4927i, 0.7697-0.2215i, 0.9223-1.7751i, -0.9826+1.6545i),
	c(-0.1185-0.8347i, -0.591+0.4927i, 0.7697+0.2215i, 0.9223+1.7751i, -0.9826-1.6545i),
	c(0.218+0.8128i, -0.8281-0.33i, -1.526+0.7887i, 1.5164-0.5302i, 0.6196-0.7413i),
	c(0.218-0.8128i, -0.8281+0.33i, -1.526-0.7887i, 1.5164+0.5302i, 0.6196+0.7413i)
	),
	# R2 = c(1,0,0,0,2)
	V10 = rbind(
	c(-0.7472+0.6624i, 1.3684+0i, -0.7472-0.6624i, 0.5629+1.0719i, 0.5629-1.0719i),
	c(-0.4726-0.6449i, 0.8466+0i, -0.4726+0.6449i, 0.5493+1.8422i, 0.5493-1.8422i),
	c(0.5758+1.1678i, 1.248+0i, 0.5758-1.1678i, -0.6998-0.675i, -0.6998+0.675i),
	c(-0.0203+0.9575i, -0.599-0.4604i, 0.7755-0.3705i, 1.5775-1.4408i, -0.7338+1.3143i),
	c(-0.0203-0.9575i, -0.599+0.4604i, 0.7755+0.3705i, 1.5775+1.4408i, -0.7338-1.3143i),
	c(0.2929+0.7451i, -0.803-0.3106i, -1.3451+0.6949i, 2.34-0.5178i, 0.5152-0.6116i),
	c(0.2929-0.7451i, -0.803+0.3106i, -1.3451-0.6949i, 2.34+0.5178i, 0.5152+0.6116i)
	),
	# R2 = c(1,1,0,0,2)
	V11 = rbind(
	c(-0.6816+0.6782i, 1.2058+0i, -0.6816-0.6782i, 0.5787+1.208i, 0.5787-1.208i),
	c(-0.4213-0.6133i, 0.7734+0i, -0.4213+0.6133i, 0.5346+2.094i, 0.5346-2.094i),
	c(0.4993+0.8644i, 1.7734+0i, 0.4993-0.8644i, -0.886-0.5887i, -0.886+0.5887i),
	c(-0.1608+0.9233i, -0.5632-0.4736i, 0.7705-0.2968i, 1.4773-1.6499i, -0.5238+1.4969i),
	c(-0.1608-0.9233i, -0.5632+0.4736i, 0.7705+0.2968i, 1.4773+1.6499i, -0.5238-1.4969i),
	c(0.2909+0.8266i, -0.7785-0.4025i, -1.1452+0.7978i, 2.0352-0.5726i, 0.5976-0.6494i),
	c(0.2909-0.8266i, -0.7785+0.4025i, -1.1452-0.7978i, 2.0352+0.5726i, 0.5976+0.6494i)
	),
	# R2 = c(1,-1,0,0,2)
	V1n1 = rbind(
	c(-0.815+0.6336i, 1.5611+0i, -0.815-0.6336i, 0.5344+0.9574i, 0.5344-0.9574i),
	c(-0.5474-0.6729i, 0.9565+0i, -0.5474+0.6729i, 0.5691+1.5668i, 0.5691-1.5668i),
	c(0.5646+1.6294i, 0.9276+0i, 0.5646-1.6294i, -0.5284-0.6677i, -0.5284+0.6677i),
	c(0.1043+0.9233i, -0.6431-0.4419i, 0.7535-0.4692i, 1.7136-1.1814i, -0.9284+1.1693i),
	c(0.1043-0.9233i, -0.6431+0.4419i, 0.7535+0.4692i, 1.7136+1.1814i, -0.9284-1.1693i),
	c(0.288+0.6875i, -0.7899-0.2239i, -1.5555+0.6011i, 2.5938-0.4843i, 0.4636-0.5804i),
	c(0.288-0.6875i, -0.7899+0.2239i, -1.5555-0.6011i, 2.5938+0.4843i, 0.4636+0.5804i)
	),
	# R2 = c(-1,1,0,0,2)
	Vn11 = rbind(
	c(-1.04+0.7834i, 0.9117+0i, -1.04-0.7834i, 0.0842+1.1344i, 0.0842-1.1344i),
	c(-0.4058-0.7245i, 0.7146+0i, -0.4058+0.7245i, -0.4515+1.9633i, -0.4515-1.9633i),
	c(0.2519+0.9654i, 1.0961+0i, 0.2519-0.9654i, -1.3-0.3782i, -1.3+0.3782i),
	c(-0.0945+0.7799i, -0.6323-0.508i, 0.7477-0.1733i, 0.4808-1.761i, -1.5017+1.6623i),
	c(-0.0945-0.7799i, -0.6323+0.508i, 0.7477+0.1733i, 0.4808+1.761i, -1.5017-1.6623i),
	c(0.1546+0.7818i, -0.8239-0.2526i, -2.0817+0.7679i, 1.1557-0.4064i, 0.5954-0.8906i),
	c(0.1546-0.7818i, -0.8239+0.2526i, -2.0817-0.7679i, 1.1557+0.4064i, 0.5954+0.8906i)
	),
	# R2 = c(-1,-1,0,0,2)
	Vn1n1 = rbind(
	c(-1.2386+0.5191i, 1.0386+0i, -1.2386-0.5191i, 0.2193+1.0098i, 0.2193-1.0098i),
	c(-0.5645-0.9232i, 0.7884+0i, -0.5645+0.9232i, -0.3297+1.4345i, -0.3297-1.4345i),
	c(-0.1641+1.2813i, 0.8206+0i, -0.1641-1.2813i, -0.7462-0.9506i, -0.7462+0.9506i),
	c(0.0377+0.7812i, -0.7332-0.3888i, 0.8342-0.3062i, 0.7392-1.3443i, -1.8779+1.2581i),
	c(0.0377-0.7812i, -0.7332+0.3888i, 0.8342+0.3062i, 0.7392+1.3443i, -1.8779-1.2581i),
	c(0.2276+0.7033i, -0.758-0.1166i, -2.5403+0.4807i, 1.6231-0.3825i, 0.4476-0.685i),
	c(0.2276-0.7033i, -0.758+0.1166i, -2.5403-0.4807i, 1.6231+0.3825i, 0.4476+0.685i)
	),
	# R2 = c(-5/3, -1,0,0,2.5)
	# + internal permutation: 8 solutions! (sol[2,] can be excluded)
	Vn1fn1n2f = rbind(
	c(-1.4391+0.489i, 0.9708+0i, -1.4391-0.489i, 0.1203+1.0489i, 0.1203-1.0489i),
	c(-1.4391-0.489i, 0.9708+0i, -1.4391+0.489i, 0.1203-1.0489i, 0.1203+1.0489i),
	c(-0.67+1.3553i, -0.4885+1.0588i, -0.9735-1.4403i, -0.3205-0.987i, 0.7859+0.0133i),
	c(-0.67-1.3553i, -0.4885-1.0588i, -0.9735+1.4403i, -0.3205+0.987i, 0.7859-0.0133i),
	c(0.0291+0.7954i, -0.7965-0.384i, 0.8613-0.2638i, 0.5542-1.3938i, -2.3148+1.2462i),
	c(0.0291-0.7954i, -0.7965+0.384i, 0.8613+0.2638i, 0.5542+1.3938i, -2.3148-1.2462i),
	c(0.2175+0.7363i, -0.7773-0.1054i, -3.047+0.4717i, 1.4931-0.3464i, 0.447-0.7562i),
	c(0.2175-0.7363i, -0.7773+0.1054i, -3.047-0.4717i, 1.4931+0.3464i, 0.447+0.7562i)
	)
);

### Example
# - deriving: R2 = c(-2,-1,0,0,2.5);
R2 = c(-2,-1,0,0,2.5)
path = lapply(c(-5/3, -1.8, -2), function(R1) c(R1, -1,0,0,2.5));
x.all = solve.path(solve.S5HtMixed.Num, x0All$Vn1fn1n2f, path=path, debug=T)


### Case 2:
# S = 0, E5 = 2;

R2 = c(0,1,0,0,2)
x0 = x0All$V01;

x.all = solve.all(solve.S5HtMixed.Num, x0, R=R2, debug=F)
poly.calc(apply(x.all, 1, function(x) sum(x * x[c(3,4,5,1,2)]))) * 27
# round0(poly.calc(x.all)) * 27

-12473 - 37419*x - 12473*x^2 - 10*x^3 - 10*x^4 + 27*x^5 + 80.75*x^6 + 27*x^7


### Case 6:
R2 = c(-1,1,0,0,2)
x0 = x0All$Vn11
x.all = solve.all(solve.S5HtMixed.Num, x0, R=R2, debug=F)
poly.calc(apply(x.all, 1, function(x) sum(x * x[c(3,4,5,1,2)]))) * 27

###
R2 = c(5/4, -1/3,0,0,2)
x0 = x0All$V10
x.all = solve.all(solve.S5HtMixed.Num, x0, R=R2, debug=F)
poly.calc(apply(x.all, 1, function(x) sum(x * x[c(3,4,5,1,2)]))) * 27


R2 = c(1,0,0,0,2)
-2500 + 12500*x - 12472*x^2 - 100*x^3 - 1*x^4 + 9*x^5 - 27*x^6 + 27*x^7

R2 = c(1,0,0,0,3/2)
-1406.25 + 7031.25*x - 7010.25*x^2 - 75*x^3 - 1*x^4 + 9*x^5 - 27*x^6 + 27*x^7

R2 = c(1,1/3,0,0,2)
277.185 - 1.998959*x - 12505.22*x^2 - 349.9733*x^3 - 6.351853*x^4 + 11.94444*x^5 + 0.1663215*x^6 + 27*x^7
R2 = c(1,-1/3,0,0,2)
-8048.84 + 25090.59*x - 12773*x^2 + 152.4053*x^3 + 8.401235*x^4 + 12.01852*x^5 - 54.16701*x^6 + 27*x^7

R2 = c(1,1,0,0,2)
-2564 - 25342*x - 13521*x^2 - 908.5*x^3 - 13.5*x^4 + 33.5*x^5 + 58.25*x^6 + 27*x^7
R2 = c(1,-1,0,0,2)
-27436 + 51262*x - 14399*x^2 + 721.5*x^3 + 51.5*x^4 + 35.5*x^5 - 112.75*x^6 + 27*x^7
R2 = c(-1,1,0,0,2)
-2420 - 24530*x - 11377*x^2 + 784.5*x^3 - 14.5*x^4 + 38.5*x^5 + 49.25*x^6 + 27*x^7
R2 = c(-1,-1,0,0,2)
-27692 + 48850*x - 10655*x^2 - 589.5*x^3 + 44.5*x^4 + 36.5*x^5 - 103.75*x^6 + 27*x^7

f5 = function(S, E11a, E5) {
	9*S^4 + 27*E11a^2 - 2*E11a^3*S^3/E5 - 3*E11a^4*S/E5;
}
f4 = function(S, E11a, E5) {
	- S^6 + 18*E11a^2*S^2 + 4*E11a^4*S^3/E5 - 10*E11a^3 - 21*E11a*S^4 - 3*E11a^5*S/E5;
}
f0 = function(S, E11a, E5, coeff=c(-27, 0)) {
	- 50*S^3*E5 - 10*E11a^4 - 2*S^3*E11a^5/E5 + 9*S*E11a^6/E5 - 10*E11a^3*S^2 +
	+ 12*E11a^2*S^4 + coeff[1]*E11a*E5*S + coeff[2]*E11a*S^6;
}
f2 = function(S, E11a, E5) {
	- 50*S^3*E5 - 10*E11a^4 + 9*S*E11a^6/E5 + 12*E11a^2*S^4;
}
f1 = function(S, E11a, E5, coeff=c(-14, -2)) {
	 - 375*E11a*E5*S - 2*S^3*E11a^5/E5 + coeff[1]*E11a^3*S^2 + coeff[2]*E11a*S^6;
}
f1S = function(S, E11a, E5, coeff=c(0)) {
	- 50*S^3*E5 - 2*S^3*E11a^5/E5 + 9*S*E11a^6/E5 - 375*E11a*E5*S;
}
f2S = function(S, E11a, E5, coeff=-54) {
	- 10*E11a^4 - 10*E11a^3*S^2 + 12*E11a^2*S^4 + coeff[1]*E11a*S^6;
}


27*(E11a^7 + E11b^7)*E5^2 +
	# x^6
	- 27*(E11a^6 + E11b^6)*E5^2*S^2 - (E11a*E11b)^6 +
	+ 81*E11a*E11b*(E11a^5 + E11b^5)*E5^2 + 9*(E11a*E11b)^3*(E11a^3 + E11b^3)*E5*S +
	# x^5:
	+ 9*(E11a^5 + E11b^5)*E5^2*S^4 + 27*(E11a*E11b)^2*(E11a^3 + E11b^3)*E5^2 +
	- 2*(E11a*E11b)^3*(E11a^2 + E11b^2)*E5*S^3 - 3*(E11a*E11b)^4*(E11a + E11b)*E5*S +
	# x^4: Note: 1 term from x^5 contributes as well;
	- (E11a^4 + E11b^4)*E5^2*S^6 + 18*(E11a*E11b)^2*(E11a^2 + E11b^2)*E5^2*S^2 +
	+ 4*(E11a*E11b)^4*E5*S^3 - 10*(E11a*E11b)^3*(E11a + E11b)*E5^2 +
	- 21*(E11a*E11b)*(E11a^3 + E11b^3)*E5^2*S^4 +
	# x^3:
	# TODO: + ... +
	- 50*(E11a^3 + E11b^3)*E5^3*S^3 +
	- 375*(E11a*E11b)*(E11a^2 + E11b^2)*E5^3*S + ...*(E11a*E11b)*(E11a^2 + E11b^2)*E5^2*S^6 +
	+ 12*(E11a*E11b)^2*(E11a + E11b)*E5^2*S^4 - 10*(E11a*E11b)^3*E5^2*S^2 +
	# x^2:
	+ 14*(E11a^2 + E11b^2)*E5^3*S^5 +
	- 5^5*E5^4*(E11a^2 + 3*E11a*E11b + E11b^2) +
	+ 5^5*(E11a + E11b)*E5^4*S^2 - 5^4*E5^4*S^4 # = 0


######################
######################

### Robust Derivation:


x = x.all[1,]; E3 = E4 = 0;
x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4]; x5 = x[5];
s1 = x1 + x2; p1 = x1 * x2;
s2 = x3 + x4 + x5; e2 = (x3 + x4)*x5 + x3*x4; e3 = x3*x4*x5;
S = s1 + s2; E5 = p1*e3;
E11a = x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x1;
E11b = x1*x3 + x2*x4 + x3*x5 + x4*x1 + x5*x2;
E2 = E11a + E11b;


###
### Transformed P[5] System:
s1 + s2 - S # = 0
#
s1*S - s1^2 + p1 + e2 - E2 # = 0
s1*e2 - s1*p1 + S*p1 + e3 - E3 # = 0
s1*e3 + p1*e2 - E4 # = 0
p1*e3 - E5 # = 0
x5^3 + s1*x5^2 + e2*x5 - e3 - S*x5^2 # = 0
#
s1*S - s1^2 + p1 + e2 - E2 # = 0
s1*p1*e2 - s1*p1^2 + S*p1^2 - p1*E3 + E5 # = 0
s1*E5 + p1^2*e2 - p1*E4 # = 0
p1*x5^3 + p1*s1*x5^2 + p1*e2*x5 - p1*S*x5^2 - E5 # = 0


p1 = toPoly.pm("s1*S - s1^2 + p1 + e2 - E2")
p2 = toPoly.pm("s1*p1*e2 - s1*p1^2 + S*p1^2 - p1*E3 + E5")
p3 = toPoly.pm("s1*E5 + p1^2*e2 - p1*E4")
p4 = toPoly.pm("p1*x5^3 + p1*s1*x5^2 + p1*e2*x5 - p1*S*x5^2 - E5")

pR1 = solve.lpm(p1, p4, p2, xn=c("p1", "e2"))
pR2 = solve.lpm(p1, p4, p3, xn=c("p1", "e2"))
pR1 = pR1[[2]]$Rez; pR1$coeff = - pR1$coeff;
pR2 = pR2[[2]]$Rez; pR2$coeff = - pR2$coeff;
table(pR2$s1)

tmp = gcd.pm(pR1, pR2, by="s1")
pR2 = diff.pm(pR2, mult.pm(pR1, toPoly.pm("x5^3")))

# Note: coeff a == 0!
x5^2*(S - x5)*(x5^5 - S*x5^4 + E2*x5^3 - E3*x5^2 + E4*x5 - E5)*s1^2 +
	- x5^2*(S - x5)^2*(x5^5 - S*x5^4 + E2*x5^3 - E3*x5^2 + E4*x5 - E5)*s1 +
	- E5^2 + 2*E4*E5*x5 - E4^2*x5^2 - E5*E2*S*x5^2 + E5*E2*x5^3 + E4*E2*S*x5^3 + E5*S^2*x5^3 - E4*E2*x5^4 +
	- 2*E5*S*x5^4 - E4*S^2*x5^4 + E5*x5^5 + 2*E4*S*x5^5 + E2^2*S*x5^5 - E4*x5^6 - 2*E2*S^2*x5^6 +
	+ 2*E2*S*x5^7 + S^3*x5^7 - 2*S^2*x5^8 + S*x5^9 - E2*S*x5^4*E3 - E2*x5^5*E3 + S^2*x5^5*E3 +
	- x5^7*E3 + x5^4*E3^2 # = 0

