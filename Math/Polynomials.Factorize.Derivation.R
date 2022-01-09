########################
###
### Leonard Mada
### [the one and only]
###
### Multi-Variable Polynomials
### Factorize: Derivations
###
### draft v.0.1b


### Factorize Multi-Variable Polynomials
### Various Derivations

### A. Univariate Polynomials


######################

### Helper Functions

source("Polynomials.Helper.R")
# - is automatically loaded in: Polynomials.Helper.R;
# source("Polynomials.Helper.Factorize.R")


#######################
#######################

### A. Univariate Polynomials
### P[3] * P[3]

# [NO results yet]

# (x^3 + b2*x^2 + b1*x + 1) * (x^3 + c2*x^2 + c1*x + 1)
p = toPoly.pm("(x^3 + b2*x^2 + b1*x + 1) * (x^3 + c2*x^2 + c1*x + 1)");

b = c(3, -4); c = c(2, 5);
b1 = b[1]; b2 = b[2]; c1 = c[1]; c2 = c[2];
P = function(x) eval.pm(p, list(x=x, b1=b1, b2=b2, c1=c1, c2=c2));

#########
### Mod 7

### P(3) (mod 7)
(9*(3*b2 + b1)*(3*c2 + c1) - P(3)) %% 7
((3*b2 + b1)*(3*c2 + c1) + 3*P(2)) %% 7

### P(5) (mod 7)
(25*(5*b2 + b1)*(5*c2 + c1) - P(5)) %% 7
((5*b2 + b1)*(5*c2 + c1) - 2*P(5)) %% 7

### P(6) (mod 7)
((b2 - b1)*(c2 - c1) - P(6)) %% 7

### Aux: P(1) (mod 7)
((b2 + b1 + 2)*(c2 + c1 + 2) - P(1)) %% 7
### Aux: P(2) (mod 7)
(4*(2*b2 + b1 + 1)*(2*c2 + c1 + 1) - P(2)) %% 7
((2*b2 + b1 + 1)*(2*c2 + c1 + 1) - 2*P(2)) %% 7

# System =>
(4*b1*c1 + 5*b2*c2 + 3*P(2) - 3*P(6)) %% 7
(6*b1*c1 + 2*b2*c2 - 2*P(5) - 4*P(6)) %% 7
# =>
(b1*c1 + 3*b2*c2 - P(2) + P(6)) %% 7
(b1*c1 - 2*b2*c2 + 2*P(5) + 4*P(6)) %% 7
# =>
(b1*c1 + P(2) + 4*P(5)) %% 7
(b2*c2 - 3*P(2) + P(5) - 2*P(6)) %% 7
((b1*c2 + b2*c1) - 2*P(2) - 2*P(5) - P(6)) %% 7
# Aux =>
(b1*c1 + b2*c2 + (b1*c2 + b2*c1) + 2*(b1+b2+c1+c2) + 4 - P(1)) %% 7
(2*(b1+b2+c1+c2) + 4 - P(1) + 4*P(2) - 3*P(5) + 3*P(6)) %% 7
((b1+b2) + (c1+c2) + 2 - 4*P(1) + 2*P(2) + 2*P(5) - 2*P(6)) %% 7
# =>
((b1+b2)*(c1+c2) - 4*P(2) + 3*P(5) - 3*P(6)) %% 7
# =>
x = (b1 + b2);
(x^2 + (2 - 4*P(1) + 2*P(2) + 2*P(5) - 2*P(6))*x - 4*P(2) + 3*P(5) - 3*P(6)) %% 7

b1.f = function(P) (2 - 4*P(1) + 2*P(2) + 2*P(5) - 2*P(6)) %% 7;
b0.f = function(P) (- 4*P(2) + 3*P(5) - 3*P(6)) %% 7;

### P1:
b1.f(P); b0.f(P);
# x^2 + x + 0 = 0 =>
# x =  0: c1 + c2 = 0 OR
# x = -1: b1 + b2 = -1;
# =>
b1 + b2 + 1 # = 0
c1 + c2 # = 0
b1*c1 + 1 # = 0
b2*c2 - 1 # = 0
# =>
(1 + b1)*(c1) - 1 # = 0 (Mod 7)
c(c1 - 2, c2 + 2) %% 7 # = 0
c(b1 + 4, b2 + 4) %% 7 # = 0


### P2: NO Solution
P2 = function(x, b1=2) x^6 + b1*x + 1;
b1.f(P2); b0.f(P2);
# x^2 + x + 5 = 0
sapply(seq(7), function(x) (x^2 + x + 5) %% 7)
x = c(1, 5)
((x - 1)*(x - 5)) %% 7
# =>
b1 + b2 - 1 # = 0
c1 + c2 - 5 # = 0
b1*c1 + 5 # = 0
b2*c2 + 1 # = 0
# =>
(1 - b1)*(5 - c1) + 1 # = 0
b1*c1 - 5*b1 - c1 + 6 # = 0
5*b1 + c1 - 1 # = 0
# =>
5*b1^2 - b1 - 5 # = 0
b1^2 - 3*b1 - 1 # = 0
sapply(seq(7), function(x) (x^2 - 3*x - 1) %% 7)
# NO Solution
# => NOT Factorisable!



#########
### Mod 9
# - offers some interesting reductions: 2^3 + 1, 8^3 + 1;
# - BUT NO pathway to solve;

### P(2) (mod 9)
(4*(2*b2 + b1)*(2*c2 + c1) - P(2)) %% 9
((2*b2 + b1)*(2*c2 + c1) + 2*P(2)) %% 9

### P(8) (mod 9)
((b2 - b1)*(c2 - c1) - P(8)) %% 9

### P(1/2) (mod 9)
((-2*b2 - 4*b1)*(-2*c2 - 4*c1) - P(5)) %% 9
((b2 + 2*b1)*(c2 + 2*c1) + 2*P(5)) %% 9

### Aux: P(1) (mod 9)
((b2 + b1 + 2)*(c2 + c1 + 2) - P(1)) %% 9
### Aux: P(4) (mod 9)
(4*(b2 - 2*b1 - 1)*(c2 - 2*c1 - 1) - P(4)) %% 9
((b2 - 2*b1 - 1)*(c2 - 2*c1 - 1) + 2*P(4)) %% 9
### Aux: P(3) (mod 9)
((3*b1 + 1)*(3*c1 + 1) - P(3)) %% 9

# System =>
(b1*c1 + 4*b2*c2 + 2*(b2*c1 + b1*c2) + 2*P(2)) %% 9
(4*b1*c1 + b2*c2 + 2*(b2*c1 + b1*c2) + 2*P(5)) %% 9
(b1*c1 + b2*c2 - (b2*c1 + b1*c2) - P(8)) %% 9
(b1*c1 + b2*c2 + (b2*c1 + b1*c2) + 2*(b1+b2+c1+c2) + 4 - P(1)) %% 9
(4*b1*c1 + b2*c2 - 2*(b2*c1 + b1*c2) + (2*b1 - b2 + 2*c1 - c2) + 1 + 2*P(4)) %% 9
(3*(b1+c1) + 1 - P(3)) %% 9
# Linear System =>
(P(2) + P(5) - 2*P(8)) %% 9 # == 0!
# Sum: Eq 3, Eq 4 =>
(2*(b1*c1 + b2*c2) + 2*(b1+b2+c1+c2) + 4 - P(1) - P(8)) %% 9
((b1*c1 + b2*c2) + (b1+b2+c1+c2) + 2 + 4*P(1) + 4*P(8)) %% 9
# Eq  5 - Eq 6 =>
(4*b1*c1 + b2*c2 - 2*(b2*c1 + b1*c2) - (b1 + b2 + c1 + c2) + 2*P(4) + P(3)) %% 9
# Eq 1 + (Eq  5 - Eq 6) =>
(5*b1*c1 + 5*b2*c2 - (b1 + b2 + c1 + c2) + 2*P(2) + 2*P(4) + P(3)) %% 9
(2*b1*c1 + 2*b2*c2 + 2*(b1+b2+c1+c2) + 4 - P(1) - P(8)) %% 9
# =>
(3*(b1+b2+c1+c2) + 8 - 2*P(1) - 2*P(8) + 2*P(2) + 2*P(4) + P(3)) %% 9
(3*(b2+c2) + 7 - 2*P(1) - 2*P(8) + 2*P(2) + 2*P(4) + 2*P(3)) %% 9

#
(5*(b1*c1 + b2*c2) + 4*(b2*c1 + b1*c2) + 2*P(2) + 2*P(5)) %% 9
((b1*c1 + b2*c2) - (b2*c1 + b1*c2) + 4*P(2) + 4*P(5)) %% 9 # redundant


# Always True:
P2 = function(x, b1=-2) x^6 - b1*x + 1
(P2(2) + P2(5) - 2*P2(8)) %% 9

p2 = toPoly.pm("(x^3 + b2*x^2 + b1*x + 1) * (x^3 + c2*x^2 + c1*x + 1) * (x^2 - 5*x + 1)");
P2 = function(x) eval.pm(p2, list(x=x, b1=b1, b2=b2, c1=c1, c2=c2));
p3 = toPoly.pm("(x^8 + x^2 - 2*x + 1)");
P3 = function(x) eval.pm(p3, list(x=x));
# Always True:
(P2(2) + P2(5) - 2*P2(8)) %% 9
(P3(2) + P3(5) - 2*P3(8)) %% 9

