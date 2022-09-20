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


### Formulas:

# - Formulas & derivations;
# - Useful for S5 Ht Systems;

# - Applicable for systems described in:
#   TODO


### Sections

### Basic:
# A.) E11a
# B.) Higher Poers


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

E2 = x1*(S - x1) + x2*(x3 + x4 + x5) + x3*(x4 + x5) + x4*x5;
E3 = x1*x2*(x3 + x4 + x5) + x3*x4*(x1 + x2 + x5) + x1*x4*x5 + x3*x5*(x1 + x2);
E4 = x1*x2*x3*x4 + x1*x2*x3*x5 + x1*x3*x4*x5 + x2*x3*x4*x5;

E11a = x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x1;
E11b = x1*x3 + x2*x4 + x3*x5 + x4*x1 + x5*x2;

E11a + E11b - E2 # = 0

### Note:
# - the remaining cyclic permutations equal E11b & E11a;
#   Perm(S5, by = 3) = rev(E11b) = E11b;
#   Perm(S5, by = 4) = rev(E11a) = E11a;

#######################

# TODO

