########################
###
### Leonard Mada
### [the one and only]
###
### k-SAT
###
### v.0.1a


to.matrix.SAT = function(m) {
	len = max(m)
	apply(m, 2, function(x) {
		r = rep(0, len)
		r[abs(x)] = sign(x)
		return(r)
	})
}
summary.SAT = function(m) {
	# TODO
}


n = 30 # Clauses
v = 26 # Variables
k = 3 # Variables per clause
### generate Matrix
m = sapply(seq(n), function(id) sample(v, k))
ms = to.matrix.SAT(m)
### Negated
# - TODO: non-sparsely negated;
nn = sample(seq(n*k), 1)
id.nn = sample(seq(n*k), nn)
ms[id.nn] = - ms[id.nn]

ms
table(ms)

####################
####################

### SAT Algorithms

### Step 1:
# - count number of occurrences of each variable:
#   positive & negated;


### Step 2:
# - select variable with highest total count;
# - ties: select variable with smallest difference
#   between positive & negated counts;


### Step 3:
# - select value for this variable:
#   TRUE if positive >= negated, otherwise FALSE;


### Step 4:
# - evaluate reduced k-SAT;
# - resolve remaining variables;
# - repeat Step 1 if NO further resolution possible;

# TODO: more work;

