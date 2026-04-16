

####################

### Helper Functions

source("Contrasts.Helper.R")


########################
########################

### Basic Decompositions

### Eigenvectors, Eigenvalues & QR-Decomposition

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

########################
########################

### Polynomial Contrasts

### Component 2:
# Linear Component

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


########################

### Component 3:
# Quadratic Component

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


#######################

### Component 4: Cubic

# n =  9: (5*nc^3 -  75*nc^2 + 316*nc - 330) / 6;
# n = 11: (5*nc^3 -  90*nc^2 + 451*nc - 546) / 6;
# n = 13: (5*nc^3 - 105*nc^2 + 610*nc - 840) / 5;
# n = 15: (5*nc^3 - 120*nc^2 + 793*nc - 1224) / 6;
# n = 17; (5*nc^3 - 15*9*nc^2 + 1000*nc - 1710) / 30;
#
make.poly( 5)[1:2, 4]; poly.c4c( 5) / sqrt(10) / 6;
make.poly( 7)[1:3, 4]; poly.c4c( 7) / sqrt(6*25)  / 6;
make.poly( 9)[1:4, 4]; poly.c4c( 9) / sqrt(9*110) / 6;
make.poly(11)[1:5, 4]; poly.c4c(11) / sqrt(6*715) / 6;
make.poly(13)[1:6, 4]; poly.c4c(13) / sqrt(100*143) / 6;
make.poly(15)[1:7, 4]; poly.c4c(15) / sqrt(13*3060) / 6;
make.poly(17)[1:8, 4]; poly.c4c(17) / sqrt(17*5700) / 6;
make.poly(19)[1:9, 4]; poly.c4c(19) / sqrt(213180)  / 6;
make.poly(21)[1:10, 4]; poly.c4c(21) / sqrt(432630)  / 6;


# Explicit:
make.poly( 5)[1:2, 4]; c(-1,2)   / sqrt(10);
make.poly( 7)[1:3, 4]; c(-1,1,1) / sqrt(6);
make.poly( 9)[1:4, 4]; c(-14,7,13,9)     / sqrt(9*110);
make.poly(11)[1:5, 4]; c(-30,6,22,23,14) / sqrt(6*715); # 5*6*11*13;
make.poly(13)[1:6, 4]; c(-11,0,6,8,7,4)  / sqrt(4*143);
make.poly(15)[1:7, 4]; c(-91,-13,35,58,61,49,27) / sqrt(13*3060);
make.poly(17)[1:8, 4]; c(-28,-7,7,15,18,17,13,7) / sqrt(17*228);

#######################

### Component 5:
# P[4]

make.poly( 5)[1:3, 5]; c(1,-4,6) / sqrt(70);
make.poly( 6)[1:3, 5]; c(1,-3,2) / sqrt(28);
make.poly( 7)[1:4, 5]; c(3, -7, 1,6) / sqrt(154);
make.poly( 8)[1:4, 5]; c(7,-13,-3,9) / sqrt(616);
make.poly( 9)[1:5, 5]; c(14,-21,-11, 9,18) / sqrt(2002);
make.poly(10)[1:5, 5]; c(18,-22,-17, 3,18) / sqrt(2860);
make.poly(11)[1:6, 5]; c( 6, -6, -6, -1, 4, 6) / sqrt(286);
make.poly(12)[1:6, 5]; c(33,-27,-33,-13,12,28) / sqrt(8008);
make.poly(13)[1:7, 5]; c( 99,-66, -96,-54, 11,64, 84) / sqrt(68068);
make.poly(14)[1:7, 5]; c(143,-77,-132,-92,-13,63,108) / sqrt(136136);


###
nc = seq(6); # n = 12;
(7*nc^4 - 182*nc^3 + 1565*nc^2 - 4966*nc + 4368) / 24;


#######################

### Component 6:

make.poly(10)[1:6, 6]; c(-6,14, -1,-11,-6,6) / sqrt(780);

