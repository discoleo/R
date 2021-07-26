########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems:
### Heterogeneous Symmetric S4:
### Mixed Type with Resonances
###
### draft v.0.1a


### Heterogeneous Symmetric
### Polynomial Systems: 4 Variables
### Mixed: Hetero + Symmetric

### Example:
# x1^p*x2^n*x3^n + x2^p*x3^n*x4^n + x3^p*x4^n*x1^n + x4^p*x1^n*x2^n = R1
# x1^k + x2^k + x^3^k + x^4^k = R2
# ...


###############

###############
### History ###
###############


### draft v.0.1a:
# - Proof of concept:
#   Rotations Entangled with Roots of Unity;


####################
####################

### helper functions

library(polynom)
library(pracma)

# the functions are in the file:
# Polynomials.Helper.R

### other functions
# ...


#####################
#####################

##################
### Rotations: ###
### x1^n*x2*x3 ###
##################

### Order: 3+1+1

### Test
x1^3*x2*x3 + x2^3*x3*x4 + x3^3*x4*x1 + x4^3*x1*x2 # - R1
x1^5 + x2^5 + x3^5 + x4^5 # - R2
x1^10 + x2^10 + x3^10 + x4^10 # - R3
(x1*x2)^5 + (x2*x3)^5 + (x3*x4)^5 + (x4*x1)^5 # - R4


### Debug:
R = c(1, -2, -1, 3);
x1 =  0.0451568676 + 0.1440038706i;
x2 =  0.9916018740 + 0.5052387172i;
x3 = -0.1370521599 - 1.1076829770i;
x4 =  0.7725898004 + 0.1223347600i;
x = c(x1, x2, x3, x4);

m = unity(5, all=TRUE);
x = sapply(seq_along(m), function(id) x*m[id]);
x = t(x)
x1 = x[,1]; x2 = x[,2]; x3 = x[,3]; x4 = x[,4];

x0 = x[1,];
sol1 = x0 * m[2]^c(1,2,0,4);
sol3 = x0 * m[2]^c(3,2,4,0);
x = rbind(x, sol1, sol3);
x = rbind(x, x[,c(2,3,4,1)], x[,c(3,4,1,2)], x[,c(4,1,2,3)]);
x1 = x[,1]; x2 = x[,2]; x3 = x[,3]; x4 = x[,4];



##################
##################
##################

##################
### Resonances ###
##################

### Order: 3+1+1
### 4 Variables
# i = 4; p = c(3,1,1)
# p = Divisors(3^3*5^2);
### Trivial:
p = 5;
### Non-Trivial & Combinations:
p = c(15);
# - the computed ones do NOT work;
# - but there are quasi-non-trivial solutions for p = 5;
### Ex: p = 5 =>
# new "non-trivial" solution:
# (x1,x2,x3,x4) * (m^1, m^2, m^0, m^4);
# (1,2,0), (2,0,4), (0,4,1), (4,1,2)
# (3,2,4), (2,4,0), (4,0,3), (0,3,2)
### Ex: p = 15 =>
# new "non-trivial" solution:
# (x1,x2,x3,x4) * (m^1, m^11, m^1, m^11);
# (1,11,1), (11,1,11), (1,11,1), (11,1,11)
