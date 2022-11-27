########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems
### S7: Hetero-Symmetric
### Useful Formulas
###
### draft v.0.1a


### Formulas:

# - Formulas & derivations;
# - Useful for S7 Ht Systems;

# - Applicable for systems described in:
#   TODO


### Sections

### Basic:
# A.) E11a
# B.) Higher Powers

### S7 Ht: Non-Oriented Ht
# - allow an inversion;
# - Example: E11a;

# If (x1,x2,x3,x4,x5,x6,x7) is a solution, then:
# - all 7 cyclic permutations are also solutions;
# - the triple permutation {(1,3), (4,7), (5,6)}
#   => (x3,x2,x1,x7,x6,x5,x4) is also a solution;
# - the 7 triple-permutations are cyclic permutations of each  other;


####################

### Helper Functions


source("Polynomials.Helper.R")


E2f = function(x) {
	len = length(x);
	if(len <= 1) return(0);
	if(len == 2) return(x[1]*x[2]);
	#
	S  = x[1] + x[2];
	E2 = x[1] * x[2];
	if(len == 2) {
		return(E2 + S*x[3]);
	}
	#
	for(id in seq(3, len)) {
		E2 = E2 + S*x[id];
		S  = S + x[id];
	}
	return(E2);
}

### Debug
x = sqrt(c(2,3,5,7,11,13,17));
x[1] = - x[1]; x[5] = - x[5]; x[7] = - x[7];
x1 = x[1]; x2 = x[2]; x3 = x[3]; x4 = x[4];
x5 = x[5]; x6 = x[6]; x7 = x[7];

### Notation:
S  = sum(x);
E7 = prod(x);

E2 = E2f(x);
# TODO
xs3 = x1 + x2 + x3; xs3c = S - xs3; xs4c = S - xs3 - x4;
E3 = x1*(x2*(S - x1 - x2) + x3*xs3c + x4*xs4c + x5*(x6 + x7));
E4 = 0; # TODO
E5 = 0; # TODO
E6 = E4 * sum(1/x);

E11a = x1*x2 + x2*x3 + x3*x4 + x4*x5 + x5*x6 + x6*x7 + x7*x1;
E11b = x1*x3 + x2*x4 + x3*x5 + x4*x6 + x5*x7 + x6*x1 + x7*x2;
E11c = x1*x4 + x2*x5 + x3*x6 + x4*x7 + x5*x1 + x6*x2 + x7*x3;

E11a + E11b + E11c - E2 # = 0

### Note:
# - the remaining cyclic permutations equal E11b & E11a;
#   Perm(S5, by = 4) = rev(E11c) = E11c;
#   Perm(S5, by = 5) = rev(E11b) = E11b;
#   Perm(S5, by = 6) = rev(E11a) = E11a;

#######################

# TODO

