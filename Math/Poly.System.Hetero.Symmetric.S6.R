########################
###
### Leonard Mada
### [the one and only]
###
### Polynomial Systems: S6
### Heterogeneous Symmetric
###
### draft v.0.1a


##############
### Simple ###
##############

### x[i]^n + b*x[i+1] = R

# x1^n + b*x2 = R
# x2^n + b*x3 = R
# x3^n + b*x4 = R
# x4^n + b*x5 = R
# x5^n + b*x6 = R
# x6^n + b*x1 = R

### Solution:

### Case 1: x1=x2=x3=x4=x5=x6
# P[n]: n solutions;

### Case 2: x1=x3=x5, x2=x4=x6, distinct
# P[2] o P[(n^2 - n)/2];

### Case 3: x1=x4, x2=x5, x3=x6, distinct
# P[3] o P[(n^3 - n)/3];

### Case 4: all distinct;
# P[6] o P[(n^6 - n^3 - n^2 + n)/6];

# - System is decomposable into 3 subsystems:
#   P[n] * (P[2] o P[(n^2 - n)/2]) *
#   (P[3] o P[(n^3 - n)/3]) *
#   (P[6] o P[(n^6 - n^3 - n^2 + n)/6]);


###############
### Order 2 ###
###############

### x[i]^2 + b*x[i+1] = R

### Solution:

### Case 1: all equal
# - 2 solutions;

### Case 2: x1=x3=x5, x2=x4=x6, distinct
# - P[2] o P[1]: 2 solutions;

### Case 3: x1=x4, x2=x5, x3=x6, distinct
# - P[3] o P[2]: 6 solutions;

### Case 4: all distinct;
# P[6] o P[9]: 54 solutions;
# S: P[9];


# TODO;

