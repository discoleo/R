########################
###
### Leonard Mada
### [the one and only]
###
### k-SAT
###
### v.0.1c-step1


to.matrix.SAT = function(m) {
	len = max(m)
	apply(m, 2, function(x) {
		r = rep(0, len)
		r[abs(x)] = sign(x)
		return(r)
	})
}
sat.gen = function(n, v, k=3, p=1/3) {
	### generate Matrix
	m = sapply(seq(n), function(id) sample(v, k))
	ms = to.matrix.SAT(m)
	### Negated
	if(p == 0) return(ms);
	# - non-sparsely negated;
	isVar = (ms == 1);
	len = n*k; # sum(isVar)
	if(p < 0) {
		nn = sample(seq(len), 1)
	} else {
		nn = round(len * p)
	}
	id.nn = sample(seq(n*v)[isVar], nn)
	ms[id.nn] = - ms[id.nn]
	return(ms)
}
summary.sat = function(m) {
	# TODO: more statistics
	st.total = apply(m, 1, function(x) sum(abs(x)))
	st.diff  = apply(m, 1, sum)
	data.frame("id"=seq(nrow(m)), "total"=st.total, "diff"=st.diff, "dabs"=abs(st.diff))
}

################
################

n = 30 # Clauses
v = 26 # Variables
k = 3 # Variables per clause

### generate Matrix
ms = sat.gen(n, v, k)
ms[1:6, 1:6]


### basic statistics
table(ms)
# number of negations per Clause
table(apply(ms, 2, function(x) sum(x < 0)))

st = summary.sat(ms)
# st
st[st$total > 3 & abs(st$dabs) <= 1, ]


id = order(st$total, st$dabs, decreasing=c(TRUE, FALSE))
st = st[id,]
st


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

### Init
st$val = NA

id.all = seq(nrow(st))
for(id in id.all) {
	### 1st non-specified value
	na.id = match(NA, st)
	st$val[na.id] = if(st$diff[na.id] >= 0) TRUE else FALSE;
	# TODO
}
