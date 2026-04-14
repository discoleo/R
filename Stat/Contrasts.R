
### Contrasts

# - Polynomial contrasts are generated
#   using function contr.poly;
# - Hardcoded Limit: n > 95;

### Helper Functions

# Inner function in contr.poly;
make.poly <- function(n, scores = seq(n)) {
    y = scores - mean(scores);
    X = outer(y, seq_len(n) - 1, '^');
    QR = qr(X);
    z  = QR$qr;
    z  = z * (row(z) == col(z)); # z = diag(z)
    raw = qr.qy(QR, z);
    Z <- sweep(raw, 2L, apply(raw, 2L, function(x) sqrt(sum(x^2))), 
            '/', check.margin = FALSE)
    colnames(Z) <- paste0("^", 1L:n - 1L)
    Z
}


# Note:
# if(n %% 2 == 1) poly(n)[seq(2, n-1, by = 2), (n+1)/2]
# should be probably == 0;
# FAILS for n = 23;

sapply(seq(15, 41, by=2), \(n) {
	# n = 21: still (approx) 0;
	# BUT for n = 23: jumps to -0.122;
	n2 = (n+1)/2;
	make.poly(n)[n2, c(n-3, n-1)];
})

# Separation: Probably 0 vs Non-Zero
rg = sapply(seq(3, 40), \(id, tol = 1E-9) {
	x = make.poly(id, seq(id));
	x = abs(x[x != 0]);
	isE = x <= tol;
	if(sum(isE) > 0) {
		err = x[isE];
		x = x[! isE]
	} else return(c(0, min(x)));
	c(max(err), min(x));
})
cbind(c(NA,NA), c(NA,NA), rg)

# Note:
# - lowest abs values for n = 23 is probably
#   around 1.36E-6 / 2 = 7E-7 (and therefore larger than 4.7E-7);
# - for n = 24: 3.5E-7;


### Upper Limit

# BUT: Failure of P[22] already for n = 23;
sapply(seq(100), \(id, tol = 1E-10) {
	x = make.poly(id, seq(id));
	x = abs(x[x != 0]);
	# 1E-8 may still occur for n > 26;
	sum(x <= tol);
})


###############

### Case: n = 3
m = matrix(c(1, 1, 1, -1, 0, 1, 1, 0, 1), nrow = 3);

### Eigen:
eigen(m)
#
eVal = c(2, 1i, -1i);
# Vectors:
mV = matrix(c(
	- c(2,1,3) / sqrt(14),
	c(1,-1i,-1) / sqrt(3),
	c(1, 1i,-1) / sqrt(3)), nrow = 3);
m %*% mV - mV %*% diag(eVal);


### QR:
mQR = matrix(c(
	c(-3,1,1) / sqrt(3),
	0, - sqrt(2), (1 + sqrt(3)) / sqrt(8),
	c(-1,0,1) * sqrt(c(4,0,2)/3)), nrow = 3);

# Note:
sin(pi/12) - (sqrt(6) - sqrt(2)) / 4 # == 0;
qraux = c(1 + 1/sqrt(3), 1 + (sqrt(3) - 1) / sqrt(8), sqrt(2/3));

#
sgn = c(1,1,1,-1,0,1,1,-1,1);
mP = sgn * sqrt(c(1/3,1/3,1/3, 1/2,0,1/2, 1/6,2/3,1/6));
mP = matrix(mP, nrow = 3)
make.poly(3); mP;


### Components
# Note:
# - Signs NOT set;

# Component 1: 1 / sqrt(n);
make.poly(11)[1:5, 1]; 1 / sqrt(11);

# Component 2:
# Odd & Even (but misses 1 coefficient)
n = 35;
n2 = (n-1)/2;
make.poly(n)[seq(1, n2), 2]; seq(n2, 1) / sqrt((n^3 - n) / 12);

# For even:
n = 36;
n2 = (n-1)/2; 
make.poly(n)[seq(1, n2+1), 2]; c(seq(n2, 1/2, by=-1)) / sqrt((n^3 - n) / 12);

# Explicit:
make.poly( 5)[, 2]; seq(2, 1) / sqrt( 2* 5);
make.poly( 7)[, 2]; seq(3, 1) / sqrt( 4* 7);
make.poly( 9)[, 2]; seq(4, 1) / sqrt(20* 3); # ???
make.poly(11)[, 2]; seq(5, 1) / sqrt(10*11);
make.poly(13)[, 2]; seq(6, 1) / sqrt(14*13);
make.poly(15)[, 2]; seq(7, 1) / sqrt(56* 5); # ???
make.poly(17)[, 2]; seq(8, 1) / sqrt(24*17);
make.poly(19)[, 2]; seq(9, 1) / sqrt(30*19);
make.poly(21)[1:10, 2]; seq(10, 1) / sqrt(110* 7); # ???
make.poly(23)[1:11, 2]; seq(11, 1) / sqrt(4 * 11*23);
make.poly(25)[1:12, 2]; seq(12, 1) / sqrt(4 * 13*25);
make.poly(27)[1:13, 2]; seq(13, 1) / sqrt(2 * 91* 9);
make.poly(29)[1:14, 2]; seq(14, 1) / sqrt(2 * 35*29);
make.poly(31)[1:15, 2]; seq(15, 1) / sqrt(2 * 40*31);


# Component 3: Quadratic
# Note: different Structure;

# Gen: Component 3
poly.cp3 = function(n) {
	nm  = (n+1)/2; nc = seq(1, nm, by = 1);
	# div = sqrt(nm * (nm-1) * (8*nm^3 - 12*nm^2 - 2*nm + 3) / 5);
	div = sqrt((8*nm^5 - 20*nm^4 + 10*nm^3 + 5*nm^2 - 3*nm) / 5);
	(3*nc^2 - 3*(n+1)*nc + 2*nm^2 + nm) / div;
}
poly.cp3.mpfr = function(n) {
	v1  = mpfr(1, getPrec(n));
	nm  = (n+1)/2; nc = seq(v1, nm, by = 1);
	# div = sqrt(nm * (nm-1) * (8*nm^3 - 12*nm^2 - 2*nm + 3) / 5);
	div = sqrt((8*nm^5 - 20*nm^4 + 10*nm^3 + 5*nm^2 - 3*nm) / 5);
	(3*nc^2 - 3*(n+1)*nc + 2*nm^2 + nm) / div;
}

# Only with Sign:
n = 17; # n = 35; # n = 20;
make.poly(n)[1:((n+1)/2), 3]; poly.cp3(n);

# Note:
# - Algorithm based on Matrix factorization
#   is remarkably stable up to n = 162; [for Component 3]
n = 161
z  = make.poly(n)[1:((n+1)/2), 3]
n1 = mpfr(n, 240); ze = poly.cp3.mpfr(n1);
summary(abs(z - ze))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 2.900e-21 7.708e-18 1.531e-17 2.270e-17 2.361e-17 4.534e-16 

n = 162
z  = make.poly(n)[1:((n+1)/2), 3]
n1 = mpfr(n, 240); ze = poly.cp3.mpfr(n1);
summary(abs(z - ze))
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 5.888e-19 9.664e-18 2.523e-17 2.487e-17 3.503e-17 1.156e-16
min(abs(ze)) # is of Order 9.9E-4 (not very small);

# Matrix factorization fails for n >= 163;


### Specific Cases:
n = 9;
nm = (n+1)/2; nc = seq(1, nm); # NO DIV 3;
make.poly(n)[1:nm, 3]; (3*nc^2 - 3*(n+1)*nc + 2*nm^2 + nm) / sqrt(4 * 7*9*11);
n = 11;
nm = (n+1)/2; nc = seq(1, nm);
make.poly(n)[1:nm, 3]; (3*nc^2 - 3*(n+1)*nc + 2*nm^2 + nm) / 3 / sqrt(6 * 11*13);
n = 13;
nm = (n+1)/2; nc = seq(1, nm);
make.poly(n)[1:nm, 3]; (3*nc^2 - 3*(n+1)*nc + 2*nm^2 + nm) / 3 / sqrt(11*13*14);
n = 15;
nm = (n+1)/2; nc = seq(1, nm); # NO DIV 3;
make.poly(n)[1:nm, 3]; (3*nc^2 - 3*(n+1)*nc + 2*nm^2 + nm) / sqrt(8 * 21*13*17);
n = 17;
nm = (n+1)/2; nc = seq(1, nm);
make.poly(n)[1:nm, 3]; (3*nc^2 - 3*(n+1)*nc + 2*nm^2 + nm) / 3 / sqrt(24 * 17*19);
n = 19;
nm = (n+1)/2; nc = seq(1, nm);
make.poly(n)[1:nm, 3]; (3*nc^2 - 3*(n+1)*nc + 2*nm^2 + nm) / 3 / sqrt(42 * 17*19);
n = 21;
nm = (n+1)/2; nc = seq(1, nm); # NO DIV 3;
make.poly(n)[1:nm, 3]; (3*nc^2 - 3*(n+1)*nc + 2*nm^2 + nm) / sqrt(42 * 11*19*23);
n = 23;
nm = (n+1)/2; nc = seq(1, nm);
make.poly(n)[1:nm, 3]; (3*nc^2 - 3*(n+1)*nc + 2*nm^2 + nm) / 3 / sqrt(20 * 77*23);

# Explicit:
make.poly( 5)[1:3, 3]; c(2,1,2)   / sqrt(2 * 7);
make.poly( 7)[1:4, 3]; c(5,0,3,4) / sqrt(4 * 3 * 7); # has 0;
make.poly( 9)[1:5, 3]; c(28,7,8,17,20)  / sqrt(4 * 7*9*11);
make.poly(11)[1:6, 3]; c(15,6,1,6,9,10) / sqrt(6 * 11*13);
make.poly(13)[1:7, 3]; c(22,11,2,5,10,13,14) / sqrt(11*13*14);
make.poly(15)[1:8, 3]; c(91,52,19,8,29,44,53,56) / sqrt(8 * 21*13*17);
make.poly(17)[1:9,  3]; c(40,25,12,1,8,15,20,23,24) / sqrt(24 * 17*19);
make.poly(19)[1:10, 3]; c(51,34,19,6,5,14,21,26,29,30) / sqrt(42 * 17*19);


################

### Component 4: Cubic

make.poly( 5)[1:2, 4]; c(-1,2)   / sqrt(10);
make.poly( 7)[1:3, 4]; c(-1,1,1) / sqrt(6);
make.poly( 9)[1:4, 4]; c(-14,7,13,9)     / sqrt(9*110);
make.poly(11)[1:5, 4]; c(-30,6,22,23,14) / sqrt(6*715); # 5*6*11*13;
make.poly(13)[1:6, 4]; c(-11,0,6,8,7,4)  / sqrt(4*143);

