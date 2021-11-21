########################
###
### Leonard Mada
### [the one and only]
###
### Polynomials: Helper Functions
### Tests


### fast load:
source("Polynomials.Helper.R")


### this file:
# options(warn=1); source("Polynomials.Helper.Tests.R")


####################

### Helper Functions

# Check value of a specific coefficient
checkCoeff.pm = function(p, val, pow=1, xn="x") {
	if(is.null(xn)) {
		coeff = p$coeff;
	} else {
		coeff = p$coeff[p[, xn] == pow];
	}
	if(length(coeff) != 1) stop("Wrong number of Monoms!");
	if(coeff != val) stop("Wrong value!");
	cat("Coeff: Success!\n");
	invisible(TRUE);
}
checkVal.pm = function(pval, val) {
	print(pval);
	stopifnot(pval == val);
	cat("Value: Success!\n");
	invisible(TRUE);
}
checkEmpty.pm = function(p) {
	stopifnot(nrow(p) == 0);
	cat("Empty: Success!\n");
	invisible(TRUE);
}
checkLength.pm = function(p, nr) {
	stopifnot(nrow(p) == nr);
	cat("Length: Success!\n");
	invisible(TRUE);
}
checkWarning.pm = function(wrn, len=1, txt=NULL) {
	stopifnot( ! is.null(wrn) && length(wrn) == len);
	if( ! is.null(txt)) {
		len = length(txt);
		stopifnot(names(wrn)[seq(len)] == txt);
	}
	cat("Warning: Success!\n");
	invisible(TRUE);
}


#######################
#######################

#############
### Tests ###
#############

### Section: Parser

cat("\n### Section: Parser\n\n")

### Parse:
### B0
p = toPoly.pm("1")
checkLength.pm(p, 1)
checkCoeff.pm(p, 1, x=NULL)

###
p = toPoly.pm("3")
checkLength.pm(p, 1)
checkCoeff.pm(p, 3, x=NULL)

###
p = toPoly.pm("-5")
checkLength.pm(p, 1)
checkCoeff.pm(p, -5, x=NULL)


cat("\nSection: Parse B0\n\tSuccess!\n")

##################
##################

### Multi-variable Multiplication

cat("\n### Section: Multiplication\n\n")

### x^3 + b1*x - R
pTest = data.frame(
	x = c(3,1,0),
	b1 = c(0,1,0),
	R = c(0,0,1),
	coeff = c(1,1,-1)
)
p = toPoly.pm("x^3 + b1*x - R")
pDiff = diff.pm(p, pTest)
pDiff
checkEmpty.pm(pDiff)


### TODO: p^2 - deprecate
pR = mult.pm(p)
pR
pDiff = diff.pm(pR, toPoly.pm("R^2 - 2*R*x^3 + x^6 - 2*R*x*b1 + 2*x^4*b1 + x^2*b1^2"))
checkEmpty.pm(pDiff)


### p^3
pp3 = data.frame(
	x  = c(0, 3, 1, 6, 4, 2, 9, 7, 5, 3),
	R  = c(3, 2, 2, 1, 1, 1, 0, 0, 0, 0),
	b1 = c(0, 0, 1, 0, 1, 2, 0, 1, 2, 3),
	coeff = c(-1, 3, 3, -3, -6, -3, 1, 3, 3, 1)
)
p.v = pow.pm(p, 3)
p.v
checkEmpty.pm(diff.pm(p.v, pp3))

# TODO
print.pm(p.v[,c(2,3,4,1)])


### eval 1:
x = 2; R = 2; b1 = -3;
# == 0 !
(x^3 + b1*x - R)^3
v = eval.pm(p.v, c(R, x, b1))
v
checkVal.pm(v, 0)


### eval 2:
R = 2; b1 = 3; x = -5;
# != 0!
(x^3 + b1*x - R)^3
v = eval.pm(p.v, c(R, x, b1))
v
checkVal.pm(v, -2863288)


### eval 2:
R = -2; b1 = -5; x = 2;
# != 0!
(x^3 + b1*x - R)^3
v = eval.pm(p.v, c(R, x, b1))
v
checkVal.pm(v, 0)

cat("\nSection: Multiplication\n\tSuccess!\n\n")


###################
###################

### Reduce

cat("\n### Section: Reduce\n\n")

# automatic:
p = toPoly.pm("x^3 + 0*b1*x^2 + b2*0*x + y*0 + 3*b + 2")
p
checkLength.pm(p, 3)

p = data.frame(x=3:0, y=0:3, coeff=c(1,0,0,1))
# automatic mechanism is bypassed: ?? what default ??
p = toPoly.pm(p, reduce=FALSE)
print.pm(p)
reduce0.pm(p)
checkLength.pm(p, 4)
checkLength.pm(reduce0.pm(p), 2)

pr = toPoly.pm(p, reduce=TRUE)
pr
checkLength.pm(pr, 2)

cat("\nSection: Reduce\n\tSuccess!\n\n")


###################
###################

### Advanced Parser

p1 = toPoly.pm("(x+1)^3")
p1

# (x+1)*(x+2)*...*(x+6)
pR = toPoly.pm(paste("(x+", seq(1,6), ")", collapse="*"))
pR
checkVal.pm(eval.pm(pR, -1), 0)
checkVal.pm(eval.pm(pR, -6), 0)

p2 = toPoly.pm("(x+a+b)^3")
p2 = sort.pm(p2, "x", xn2= c("a", "b"))
p2
# == 3^3
r = eval.pm(p2, c(2,4,-3))
checkVal.pm(r, 3^3)

### x^3
toPoly.pm("p1(x = x-1)")

### x^3 * (x+b)^3
toPoly.pm("p1(x = x-1) * p2(x = x-a)")

###
f = function(p1) toPoly.pm("p1(x = x-1)")
f(p1)
# (x-1)^3
f(toPoly.pm("x^3"))


### Specified Parameters
b = c(2, 5, 3)
polyGenTest = function(b, id) {
	p = toPoly.pm("x^2 + b[id]*x + 1");
	checkCoeff.pm(p, b[id], pow=1);
	print(p)
}
polyGenTest(b, 1)
polyGenTest(b, 2)
polyGenTest(b, 3)

p = toPoly.pm("x^2 + (b^2 - b + 1)[1]*x + 1")
p # ... + 3*x


### Power n
n = 3
p = toPoly.pm("x^n + b*x - R")
p


###
# [resolved] clashes with internal variable;
n = 3
m = 2
p = toPoly.pm("x^(n+m) + b*x - R")
p

f = function(m) toPoly.pm("x^(n+m) + b*x - R")
m = -1
f(0)
f(2)
f(3)

cat("\nSection: Advanced Parser\n\tSuccess!\n\n")

###############

### Eval:

# (x+1)*(x+2)*...*(x+5)
sP = paste("(x+", seq(1,5), ")", collapse="*");
pR = mult.pm(toPoly.pm(sP), toPoly.pm("a+b"));
pR
eval.pm(pR, c(-2, 1,1)) # 0
eval.pm(pR, c(0, -2,-3)) # != 0
eval.pm(pR, list(a=-6, b=6, x=2)) # 0
eval.pm(pR, list(a=-6, b=5, x=-5)) # 0
eval.pm(pR, list(a=-6, b=5, x=-6)) # != 0

cat("\nSection: Eval\n\tSuccess!\n\n")

###############

### Shift vars
p1 = toPoly.pm("a*x^3 + b*x^3 + 1")

### Test 1:
pR = shift.pm(p1, -1, "x")
pDiff = diff.pm(pR, toPoly.pm("(a+b)*(x-1)^3 + 1"))
pDiff
checkEmpty.pm(pDiff)


### Test 2:
pR = shift.pm(p1, c(-1,1), "x")
# TODO: ERROR: all possible environments tested!
# wrn = parent.frame()[["last.warning"]]; print(wrn);
# flush(stderr()); flush(stdout()); # update Warnings;
wrn = warnings(); print(wrn)
# checkWarning.pm(wrn, 1, "Same variable used!")

pDiff = diff.pm(p1, pR)
# nrow == 0 & Warning!
pDiff
checkEmpty.pm(pDiff)


### Test 3:
pR = shift.pm(p1, c(-1,1), c("a", "b"))
diff.pm(p1, pR)

### Test 4:
pR = shift.pm(p1, c(-1,2), c("a", "b"))
diff.pm(pR, toPoly.pm("(a-1)*x^3 + (b+2)*x^3 + 1"))

### Test 5:
p1 = toPoly.pm("x^3 - 1/27")
pR = shift.pm(p1, 1/3, "x")
# b0 == 0 !
pR

cat("\nSection: Shift\n\tSuccess!\n\n")


####################

### Replace vars

### with specific Value
p = toPoly.pm("(x+3)^4")
replace.pm(p, -3, xn="x")
replace.pm.numeric(p, -3, xn="x")
replace.pm.numeric(p, -2, xn="x")


p = toPoly.pm("(x+1)^3*(x-a)^3")
pR = replace.pm(p, 1, xn="a")
diff.pm(pR, toPoly.pm("(x^2-1)^3"))
# x^3 * (x+1)^3 => x^3
pR = replace.pm(p, 0, xn="a")
pR = div.pm(pR, toPoly.pm("(x+1)^3"), "x")$Rez
print.pm(pR)
# (x+a)^3 => x^3
pR = div.pm(p, toPoly.pm("(x+1)^3"), "x")$Rez
pR = replace.pm(pR, 0, xn="a")
print.pm(pR)


# x^4 + x^3 + x^2 + x + 1
p = toPoly.pm(data.frame(x=4:0, coeff=1))
p = sum.pm(p, toPoly.pm("y^2"))
p1 = replace.pm(p, unity(5, all=FALSE), xn="x")
p2 = replace.pm(p, c("x"=unity(5, all=FALSE)))
print.pm(p1)
print.pm(p2)


### Replace with character (new name)

p1 = toPoly.pm("(x+y+z)^3")
#
p2 = replaceNames.pm(p1, "y", "z")
diff.pm(p2, toPoly.pm("(x+2*y)^3"))
# Cyclic permutation:  # sequential = FALSE!
p2 = replaceNames.pm(p1, c("y","x"), xn=c("x", "y"))
diff.pm(p1, p2) # SAME!
# x => y; then y => x;
p2 = replaceNames.pm(p1, c("y","x"), xn=c("x", "y"), seq=TRUE)
diff.pm(p2, toPoly.pm("(2*x + z)^3"))
p2 = replaceNames.pm(p1, c("x", "y","x"), xn=c("z", "x", "y"), seq=TRUE)
diff.pm(p2, data.frame(x=3, coeff=3^3))


### Replace with Character / Higher Power

p1 = toPoly.pm("x^5 + c*x^3 + b0")
replace.pm(p1, "Big", "x", pow=5)

### Quintic
p1 = toPoly.pm("x^5 - 5*K*x^3 - 5*(K^2 + K)*x^2 - 5*K^3*x - K^4 - 6*K^3 + 5*K^2 - K")
# fractional powers:
r = toPoly.pm("K^(4/5) + K^(3/5) + K^(1/5)")
# - we just found a root of a non-trivial quintic!
replace.pm(p1, r, "x", pow=1)

### All roots
rootTest.pm = function(id, n=5) {
	m = unity(n, all=FALSE);
	mj = m^id;
	r = toPoly.pm("K^(4/5)*mj()^4 + K^(3/5)*mj()^3 + K^(1/5)*mj()");
	return(r)
}
replaceRoot.pm = function(id, n=5) {
	p = replace.pm(p1, rootTest.pm(id, n=n), "x", pow=1);
	p = round0.pm(p);
	p = reduce.pm(p);
	cat("\n")
	return(p)
}
#
n = 5; # fixed!
# - we just found ALL 5 roots!
lapply(seq(0, 4), replaceRoot.pm, n=n)

# more formal
r = toPoly.pm("k^4*m^4 + k^3*m^3 + k*m")
# m^5 = 1; # m = roots of unity;
pR = p1;
pR = replace.pm(pR, r, xn="x")
pR = replace.pm(pR, "K", xn="k", pow=5)
pR = replace.pm(pR, 1, xn="m", pow=5)
pR # the roots worked!

### multiple Powers
Class1Root = function(n, s=NULL) {
	sn = paste0("s", seq(n-1,1));
	r.str = paste0(sn,
		"*k^", seq(n-1, 1),
		"*m^", seq(n-1, 1), collapse="+");
	r = toPoly.pm(r.str);
	if( ! is.null(s))
		r = replace.pm(r, s, sn, pow=1);
	return(r)
}
Class1Poly = function(s, n=5) {
	# Note: non-efficient algorithm!
	if(length(s) >= n) stop("s0 not yet implemented!")
	pFactor = toPoly.pm("x - r");
	r = Class1Root(n);
	# c(s4, s3, s2, s1)
	sn = paste0("s", seq(n-1,1));
	plst = lapply(seq(0, n-1), function(id) {
		pF = replace.pm(pFactor, r, "r", pow=1);
		isValue = ! is.na(s);
		pF = replace.pm(pF, s[isValue], sn[isValue], pow=1);
		pF$m = (pF$m * id) %% n;
		return(pF);
	})
	p = mult.lpm(plst);
	p = reduce.radicals(p, n=n);
	return(p)
}
reduce.radicals = function(p, n=5) {
	# check 2 Vars:
	p = replace.pm(p, c("K", "M"), c("k", "m"), pow=c(n,n));
	p = replace.pm(p, 1, "M", pow=1);
	if("m" %in% names(p))
		p = replace.pm(p, data.frame(m=seq(0, n-2), coeff=-1), "m", pow=n-1);
	if(nrow(p) == 0) return(p); # Coeffs canceled!
	p = toPoly.pm(p); p = sort.pm(p, "x", xn2="K");
	return(p);
}

n = 5;
s = c(1, 0,-2,1)
p = Class1Poly(s, n=n)
print.pm(p, lead="x")
# Check:
pT = replace.pm(p, Class1Root(n, s), "x")
pT = reduce.radicals(pT, n=n)
pT
#
print.pm(replace.pm(p, 2, "K"), lead="x")


###
n = 5
s = c(1,0,NA,1)
p = Class1Poly(s, n=n)
print.pm(p, lead="x")
K = 3
pR = replace.pm(p, c(K,-K), c("K","s2"))
print.pm(pR, lead="x")
pR = replace.pm(p, c(K=K, s2=-K))
print.pm(pR, lead="x")
#
err = eval.pm(pR, sum(c(1,0,-K,1)*rootn(K^(4:1), n)))
round0(err)


########################

cat("\n### Section: Replace Vars\n\n")

###
p = toPoly.pm("(x*y + 2)^4 + b3*(x*y + 1)^3")
pR = replace.pm.character.pm(p, "xy", toPoly.pm("x*y"))
pR = sortColumns.pm(pR)
pR
replace.pm(pR, c(xy=-2)) # -b3
replace.pm(pR, c(xy=-1)) # 1
replace.pm(pR, toPoly.pm("b3 - 1"), "xy") # 2*b3^4 + ...
replace.pm(pR, toPoly.pm("-b3 - 2"), "xy") # 0*b3^4 - 3*b3^3 - ...


### Monomial: 2*x*y => "xy"
p = toPoly.pm("(2*x*y + 2)^4 + b3*(6*x*y + 5)^3")
pR = replace.pm.character.pm(p, "xy", toPoly.pm("2*x*y"))
pR = sortColumns.pm(pR)
pR
replace.pm(pR, c(xy=-1)) # 2*x*y = -1 => 8*b3 + 1;
replace.pm(pR, c(xy=-2)) # => - b3;


###
p = toPoly.pm("(x + y + b0)^3")
pR = replace.pm.character.pm(p, "xy2", toPoly.pm("x*y^2"))
pR = replace.pm.character.pm(pR, "x2y", toPoly.pm("x^2*y"))
pR


###
p = toPoly.pm("(x+y)^4 + b3*(x+y)^3")
pR = replace.pm.character.pm(p, "xy", toPoly.pm("x*y"))
pR


########################
########################

### Multiply List of Polynomials

cat("\n### Section: Poly List\n\n")

###
p = lapply(seq(5), function(n) toPoly.pm("x^n - 1"))
pR = mult.lpm(p)
print.pm(pR)


cat("\n### Section: Complex Numbers\n\n")

p[[6]] = 2 - 1i;
pR = mult.lpm(p)
print.pm(pR)

sapply(seq(1, 5), function(n) {
	m = unity(n, all=FALSE);
	err = eval.pm(pR, list(x=m));
	round0(err);
})


pR = mult.pm(pR, toPoly.pm("x^2 - (1 + 1i)*x - 2 - 1i"))
print.pm(pR)
print.pm(pR, brackets.complex=FALSE)
B0.pm(pR) # B0 = 5;


########################
########################

cat("\n### Section: Parse Characters as Vars\n\n")

### Symmetric Polynomials
n = 3
b = paste0("b", seq(n))
plst = lapply(b, function(b) toPoly.pm("x^2 + b[1]*x + 1"))
p = mult.lpm(plst)
p = sort.pm(p, xn="x")
# p
checkEmpty.pm(diff.pm(p, rev.pm(p, xn="x")))


########################
########################

cat("\n### Section: Differentiation\n\n")

### D
n = 3
p1 = toPoly.pm("x^n + b1*x + b0")
p2 = toPoly.pm("x^n + c2*x^2 + c1*x")
p = dp.exp.pm(list(Poly=p1, Exp=p2))
p$Poly = sort.pm(p$Poly, "x", xn2=c("c2","b1","c1"))
print.pm(p$Poly, leading="x")

# 3*x^5 + 2*c2*x^4 + c1*x^3 + 3*b1*x^3 + 3*x^2 + 3*b0*x^2 + 2*c2*b1*x^2 +
#	+ 2*c2*b0*x + c1*b1*x + c1*b0 + b1


###
cat("\n\n All Tests: Success!\n\n")

