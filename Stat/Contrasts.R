
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
# if(n %% 2 == 2) poly[seq(2, n, by = 2), (n+1)/2]
# should be probably == 0;
# FAILS for n = 23;


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


# Note:
# Componenet 1: 1 / sqrt(n);
make.poly(11)[1:5, 1]; 1 / sqrt(11);

# Component 2:
make.poly( 9)[, 2]; seq(4, 1) / sqrt(20* 3); # ???
make.poly(11)[, 2]; seq(5, 1) / sqrt(10*11);
make.poly(13)[, 2]; seq(6, 1) / sqrt(14*13);
make.poly(15)[, 2]; seq(7, 1) / sqrt(56* 5); # ???
make.poly(17)[, 2]; seq(8, 1) / sqrt(24*17);
make.poly(19)[, 2]; seq(9, 1) / sqrt(30*19);
make.poly(21)[, 2]; seq(10, 1) / sqrt(110* 7); # ???
make.poly(23)[, 2]; seq(11, 1) / sqrt( 44*23);

