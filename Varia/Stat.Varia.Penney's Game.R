
### Statistics: Varia
### Penney's Game
###
### Leonard Mada
### 2020-01-16

# draft v0.2:
# - added Infinite Series;

#####################
### Penney's Game ###

# A coin is given with:
# p(H) = p, p(T) = 1 - p;

# Calculate:
# A.) E(HHT), E(HTH);
# B.) p, so that E(HHT) = E(HTH);
# C.) Infinite Series: sum( (i+3)/x^i );


### A.) Computing the Expectations
### A.a.) Mathematical approach
### A.b.) TODO: Simulation approach

### A.a.) Mathematical approach

### HHT
# let E(HHT) = x; X = "HHT";
#
# x = p * E(X|H) + (1-p) * E(X|T);
#  E(X|T) = x + 1;
#  E(X|H) = p * E(X|HH) + (1-p) * E(X|HT);
#    E(X|HT) = x + 2;
#    E(X|HH) = 3 * (1-p) + p * (E(X|HH) + 1);
#    => E(X|HH) = (3 - 2*p) / (1-p);
# after some calculations:
# x = 1 / ( p^2 * (1-p) )

### HTH
# let E(HTH) = x; X = "HTH";
#
# x = p * E(X|H) + (1-p) * E(X|T);
#  E(X|T) = x + 1;
#  E(X|H) = p * E(X|HH) + (1-p) * E(X|HT);
#    E(X|HH) = E(X|H) + 1; # needs back-substitution;
#    E(X|HT) = 3 * p + (1-p) * (x + 3);
#    => E(X|HT) = 3*p + (1-p) * (x + 3);
#  => E(X|H) = (1-p)*(x+3) + 3*p + p/(1-p);
# after some calculations:
# x = (-p^2 + p + 1) / ( p^2 * (1-p) )

### Test
p = 1/2
pHHT = 1 / ( p^2 * (1-p) )
pHTH = (-p^2 + p + 1) * pHHT
pHHT
pHTH

p = 0.4
pHHT = 1 / ( p^2 * (1-p) )
pHTH = (-p^2 + p + 1) * pHHT
pHHT
pHTH

p = 1/4
pHHT = 1 / ( p^2 * (1-p) )
pHTH = (-p^2 + p + 1) * pHHT
pHHT
pHTH

### A.b.) Simulation approach

### TODO


### B.) p, so that E(HHT) = E(HTH)

# pHHT = (-p^2 + p + 1) * pHHT;
# => -p^2 + p + 1 = 1
# => p^2 - p = 0
# p = 0 or p = 1;

### both are extreme probability values;

####################

### C.) Infinite Series

# Series:
# sum( (i+3)/x^i ), from i=1 to Infinite;
# where x = 1/p;

# Ex.) 4/3 + 5/3^2 + 6/3^3 + 7/3^4 + ...;

### Estimated sum
sum.series = function(p) {
	i = 1: 25
	sum((i+3) * p^i)
}

### Exact sum
sum.poly = function(p) {
	# ((3-2*p)/(1-p) + 1) * p/(1-p)
	p * (4-3*p) / (1-p)^2
}

p = 1/2
sum.series(p)
# 4.9999999
sum.poly(p)
# 5

p = 2/9
sum.series(p)
# 1.22449
sum.poly(p)
# 1.22449
