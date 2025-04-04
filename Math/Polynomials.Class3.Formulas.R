
## Polynomials: Class 3
## Formulas

# "Nightmare on Elm Street"
# Any students interested in eternal fame?


###################

### Helper Functions

# source("Polynomials.Helper.R")

source("Polynomials.Helper.EP.R")


####################
####################

###############
### Order 5 ###
###############

r = 2*cos(2*pi*seq(5)/11)
x = r;

#
cc = c(1,-2,3,-4)

### Coeff b4
# Direct Formula:
sum(c(5, -1, 9, -4, 25) * c(0, cc))

###
# b4 = - sum(cc * S);
x  = r;
s1 = sum(x)
s2 = sum(x^2)
s3 = sum(x^3)
s4 = sum(x^4)
Sn = round(c(s1,s2,s3,s4))
S  = Sn[1];
b4 = - sum(cc * Sn)
print(Sn); print(b4)

# General:
E2 = (S^2 - s2)/2;
E3 = - (S^3 - s3 - 3*E2*S)/3;
E4 = (S^4 - s4 - 4*S^2*E2 + 2*E2^2 + 4*S*E3)/4;
E5 = prod(x);

### Coeff b3
sum(c(100,-9,147,-36,184) * cc[4] * c(0, cc),
	- c(16, 21, 20, 31) * cc[3] * c(0, cc[1:3]),
	  c(36, -5, 28) * cc[2] * c(0, cc[1:2]),
	- c(4, 4)*cc[1]*c(0, cc[1]), 10*0^2);


### Derivation:
# see:
# print.pm(Epoly.adv(pow, 5, 2));
# print.pm(Epoly.distinct(c(pow), 5));
e22 = E2^2 - 2*S*E3 + 2*E4;
e33 = E2^3 - 3*S*E2*E3 + 3*E3^2 + 3*S^2*E4 - 3*E2*E4 - 3*S*E5;
e44 = E2^4 - 4*S*E2^2*E3 + 2*S^2*E3^2 + 4*E2*E3^2 + 4*S^2*E2*E4 +
	- 4*E2^2*E4 - 8*S*E3*E4 + 6*E4^2 - 4*S^3*E5 + 8*S*E2*E5 - 4*E3*E5;
#
e21 = E2*S - 3*E3;
e31 = E2*S^2 - E3*S - 2*E2^2 + 4*E4;
e41 = E2*S^3 - E3*S^2 - 3*E2^2*S + E4*S + 5*E2*E3 - 5*E5;
e32 = - 2*E3*S^2 + E2^2*S + 5*E4*S - E2*E3 - 5*E5;
e42 = - 2*E3*S^3 + E2^2*S^2 + 2*E4*S^2 + 4*E2*E3*S - 6*E5*S +
	- 2*E2^3 - 3*E3^2 + 2*E2*E4;
e43 = 3*E4*S^3 - 3*E2*E3*S^2 - 7*E5*S^2 + E2^3*S + 5*E3^2*S +
	- 2*E2*E4*S - E2^2*E3 - 5*E3*E4 + 7*E2*E5;
#
b3 = sum(cc^2 * c(E2,e22,e33,e44), cc[1]*cc[-1]*c(e21,e31,e41),
	cc[2]*cc[c(3,4)]*c(e32,e42), cc[3]*cc[4]*e43);
print(b3)



### Coeff b2
# TODO


###
cc = c(1,-2,3,-4)
x = rbind(cc) %*% rbind(r,r^2,r^3,r^4);
poly.calc0(x)
x^5 + 131*x^4 + 4484*x^3 + 48403*x^2 + 56433*x + 17071


###
cc = c(5,-1,0,1)
x = rbind(cc) %*% rbind(r,r^2,r^3,r^4);
poly.calc0(x, digits = 3)
x^5 - 11*x^4 - 55*x^3 + 2*121*x^2 + 3*121*x - 121


### Coeff b0: The Monster
# see file: Polynomials.Class3.R
# Note: s1,...,s4 in that file represents cc;

