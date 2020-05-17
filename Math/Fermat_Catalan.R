
### Leonard Mada
###
### Fermat's Theorem & Catalan's Conjecture
###
### draft 0.1

### Topics:
### A.) Decomposition of sums of Powers into Products
### Applications:
### B.) Fermat's Theorem
### C.) Catalan's Conjecture

################
###
###  Theory  ###
###
################

### A.) Decomposition of sums of Powers into Products

# a^n - b^n = prod(a - b*m^j)
# where m^n = 1, m = root of unity;
# j = 0 to n-1;


### B.) Decomposition of Fermat's Theorem

# a^n + b^n = c^n
# =>
# prod(a+b - c*m^j) == sum(choose(n, k) * a^k * b^(n-k)),
# where the product goes from j=0 to j=n-1;
# m^n = 1, m = root of unity;
# the sum goes from *1* to *(n-1)*!

# if SOLUTION, then: prod() == sum()

# Notes:
# - if a,b > 1 and n = prime => sum is divisible by n*a*b;

### TODO:
# 1.) Prove Fermat's Theorem using this approach;


### C.) Catalan's Conjecture

# a^m - b^n = 1
# =>
# prod(a^(1/n) - b^(1/m)*w^j) = 1
# where w^(m*n) = 1, w = root of unity;

# this is a generalization of the decomposition formula from [A];

# Notes:
# - the decomposition is a product:
#   it should be easier to work with products;
#   [high hopes to solve the problem using this approach]

####################
####################

### helper functions

unity = function(n, all=TRUE) {
	m = complex(re=cos(2*pi/n), im=sin(2*pi/n))
	if(all) {
		m = m^(0:(n-1))
	}
	return(m)
}

# Test Fermat
fermat = function(n, a, b, c) {
	id = 1:(n-1)
	# only if: a^n + b^n = c^n
	s = sum(choose(n,id) * a^id * b^(n-id))
	# but that is not the case (for n > 2)
	r.corect = s + a^n + b^n - c^n
	return(list("sum"=s, "actual"=r.corect))
}

########################


#############################
###
### A.) Basic Decompositions:
### Sums of Powers into Products
###
#############################

### n = 5
n = 5
m = unity(n)

# a=3, b=2
a = 3
b = 2
a^n - b^n
prod(a - b*m)
# 211


# a=5, b=3
a = 5
b = 3
a^n - b^n
prod(a - b*m)
# 2882

### Generalization
n = 15
m = unity(n)

# a^5 - b^3
a = 5
b = 3
a^5 - b^3
prod(a^(5/15) - b^(3/15)*m)
# 3098

# a^5 - b^3
a = 3
b = 4
a^5 - b^3
prod(a^(5/15) - b^(3/15)*m)
# 179


########################
###
### B.) Fermat's Theorem
###
########################

### some Tests

### n = 5
n = 5
m = unity(n)

# a=2, b=3, c=4
prod(2 + 3 - 4*m)
# 2101
f.s = fermat(n, 2, 3, 4)
f.s
# sum = 2850: NO solution to Fermat's problem;
# corrected: 2101 # the calculations are correct ;-)
f.s$actual - f.s$sum
# - 749
2^n + 3^n - 4^n
# - 749


# a=5, b=6, c=7
prod(5 + 6 - 7*m)
# 144244
fermat(n, 5, 6, 7)
# sum = 150150: NO solution to Fermat's problem;
# corrected: 144244

### TODO:
# find a suitable methodology to prove Fermat's Theorem;
# Note:
# - sum is divisible by n*a*b;
# - Q: Can this fact be used?


########################
########################


########################
###
### Catalan's Conjecture
###
########################


### Working solution
n = 6
m = unity(n)

# a^3 - b^2
a = 3
b = 2
a^2 - b^3
prod(a^(2/6) - b^(3/6)*m)
# 1 ### THIS IS THE SOLUTION !


### Other Solutions?
# lets find more solutions ;-)

# a^3 - b^2
a = 5
b = 3
a^2 - b^3
prod(a^(2/6) - b^(3/6)*m)
# -2 # close, but NO solution;


