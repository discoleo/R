########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S5: Hetero-Symmetric
### Useful Formulas
###
### draft v.0.1a-test


### Formulas:

# - Formulas & derivations;
# - Useful for S5 Ht Systems;

# - Applicable for systems described in:
#   TODO


### Sections

### Basic:
# A.) E11a
# B.) Higher Powers


####################

### Helper Functions


source("Polynomials.Helper.R")


### Debug
x = sqrt(c(2,3,5,7,11));
x[1] = - x[1]; x[5] = - x[5];
x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4]; x5 = x[5];

### Notation:
S  = sum(x);
E5 = prod(x);

s1 = x1 + x3; s2 = x2 + x4;
p1 = x1 * x3; p2 = x2 * x4;
ps = s1 * s2; sp = p1 + p2;
# E2 = x1*(S - x1) + x2*(x3 + x4 + x5) + x3*(x4 + x5) + x4*x5;
E2 = sp + ps + x5*(S - x5);
# E3 = x1*x2*(x3 + x4 + x5) + x3*x4*(x1 + x2 + x5) + (x1 + x2)*(x3 + x4)*x5;
E3 = p1*s2 + p2*s1 + x5*(sp + ps);
# E4 = x1*x2*x3*x4 + x1*x2*x3*x5 + x1*x2*x4*x5 + x1*x3*x4*x5 + x2*x3*x4*x5;
E4 = p1*p2 + x5*(p1*s2 + p2*s1);

E11a = x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x1;
E11b = x1*x3 + x2*x4 + x3*x5 + x4*x1 + x5*x2;

E11a + E11b - E2 # = 0

### Note:
# - the remaining cyclic permutations equal E11b & E11a;
#   Perm(S5, by = 3) = rev(E11b) = E11b;
#   Perm(S5, by = 4) = rev(E11a) = E11a;

#######################

# TODO


#######################

### Disaster:
pE1 = toPoly.pm("x1 + x2 + x3 + x4 + x5 - 1"); # S = 1; !!!
pE2a = toPoly.pm("x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x1"); # E11a = 0; !!!
pE2b = toPoly.pm("x1*x3 + x2*x4 + x3*x5 + x4*x1 + x5*x2"); # E11b = 0; !!!
# E3 = 0; E4 = 0; !
pE3 = toPoly.pm("x1*x2*(x3 + x4 + x5) + x3*x4*(x1 + x2 + x5) + (x1 + x2)*x4*x5 + x3*x5*(x1 + x2)");
pE4 = toPoly.pm("x1*x2*x3*x4 + x1*x2*x3*x5 + x1*x2*x4*x5 + x1*x3*x4*x5 + x2*x3*x4*x5");

pR = solve.lpm(pE1, pE2a, pE2b, pE3, pE4, xn=c("x5", "x4", "x3", "x2"))


### Overflows!
pE1 = toPoly.pm("s1 + s2 + x5 - 1"); # S = 1; !!!
pE2 = toPoly.pm("p1 + p2 + s1*s2 + x5*(1 - x5)"); # E2 = 0; S = 1;
pE3 = toPoly.pm("p1*s2 + p2*s1 + x5*(s1*s2 + p1 + p2)");
pE4 = toPoly.pm("p1*p2 + x5*(p1*s2 + p2*s1)");
pE5 = toPoly.pm("p1*p2*x5 - 1"); # E5 = 1;

pR = solve.lpm(pE1, pE5, pE2, pE3, pE4, xn=c("x5", "s2", "s1", "p2"))
max(pR[[4]]$Rez$p1)


###
# spp = sp + ps; pp4 = p1*p2;
# - Solvable, but just ordinary P[5];
pE2 = toPoly.pm("spp + x5*(S - x5) - E2");
pE3 = toPoly.pm("A + x5*spp - E3");
pE4 = toPoly.pm("pp4 + x5*A - E4");
pE5 = toPoly.pm("pp4*x5 - E5");

pE4 = toPoly.pm("x5^2*(E3 - x5*spp) + E5 - E4*x5")

pR = solve.pm(pE2, pE4, xn=c("spp"))
pR = pR$Rez;
print.pm(pR, lead="x5") # trivial: ordinary P[5];

