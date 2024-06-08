

### Complex Symmetries
### S 3x3 Ht


####################

### Helper Functions

as.variables = function(x, n=4) {
	tmp = which(x != 0);
	id1 = (tmp - 1) %/% n + 1;
	id2 = (tmp - 1) %% n + 1;
	ss = letters[seq(n)];
	sn = ss[id1];
	sn = paste0(sn, id2);
	return(sn);
}


#####################

### System 3x3 Ht

# x1^n + x2^n + x3^n = R1
# y1^n + y2^n + y3^n = R1
# z1^n + z2^n + z3^n = R1
#
# x1^2 + y1^2 + z1^2 = R2
# x2^2 + y2^2 + z2^2 = R2
# x3^2 + y3^2 + z3^2 = R2
#
# x1*y2 + y2*z3 + z3*x1 = R3
# x2*y3 + y3*z1 + z1*x2 = R3
# x3*y1 + y1*z2 + z2*x3 = R3

### Note:
# If ((x1,x2,x3), (y1,y2,y3), (z1,z2,z3)) is a solution:
# - then any simultaneous cyclic permutation of
#   the rows or the columns or the diagonals are also a solution;
# => 9 permutations;
# - this particular type & n = even: (-1) * solutions is also solution;

# If n = 2: => R1 = R2, but system is indeterminate;
# Aux eq: x1*y1*z1 + x2*y2*z2 + x3*y3*z3 = R4

# TODO: solve;


#####################
#####################

### Matrix
### Complex Symmetries

### Size = 3
m = matrix(0L, nrow=9, ncol=9)
for(i in seq(3)) m[i, c(1,2,3) + (i-1)*3] = 1L;
for(i in seq(3)) m[i + 3, c(0,3,6) + i] = 1L;
for(i in seq(3)) m[i + 6, c(i, i %% 3 + 4, (i+1) %% 3 + 7)] = 1L;

det(m)


### Size = 4
m = matrix(0L, nrow=16, ncol=16)
for(i in seq(4)) m[i, c(1,2,3,4) + (i-1)*4] = 1L;
for(i in seq(4)) m[i + 4, c(0,4,8,12) + i] = 1L;
for(i in seq(4)) m[i + 8, c(i, i %% 4 + 5, (i+1) %% 4 + 9, (i+2) %% 4 + 13)] = 1L;
for(i in seq(4)) m[i + 12, c((i+2) %% 4 + 1, (i+1) %% 4 + 5, i %% 4 + 9, (i-1) %% 4 + 13)] = 1L;

det(m)

tmp = apply(m, 1, function(x) paste0(as.variables(x), collapse = " + "))

