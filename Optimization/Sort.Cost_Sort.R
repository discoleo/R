##################
###
### Cost-Sort
###
### Leonard Mada
###
### draft v.0.1b


### Problem: Optimize the Cost-Sort sort

### Data
# S = set of items X[i];
### Properties of X:
# X[i]: can be ranked/ordered;
# X[i]: each x[i] has an associated cost;

### Cost:
# Cost(Swap(x[i], x[j])) = mean(Cost(x[i]), Cost(x[j]));
# Cost(Sort(S)) = sum(Cost(Swaps(sorting of S)));

### Task:
# Minimize(Cost(Sort(S)));

### Ties:
# Variant 1: NO ties;
# Variant 2: all equivaent ties have the same cost;
# Variant 3: S may contain ties, i.e. x[i] = x[j], with different cost;

###############

###############
### History ###

### draft v.0.1b:
# - added costsort.gen() function:
#   generates data sets;


################
################

### Analysis

# TODO

# Can we compute the optimal value?
# ;-)

# minimal Cost = minimum(S => Sorted(S));


####################

### helper Functions

### generate data-sets
costsort.gen = function(n, unique=TRUE, scale=1/4, doRank=TRUE) {
	id = 1:n
	x = sample(id, n, replace = ! unique)
	### possible options
	# x.max = n
	# x = runif(n, 0, x.max)
	### Cost
	cost = sample(id, n) * scale
	x.df = data.frame(id=id, x=x, cost=cost)
	if(doRank) {
		x.df$rank = rank(cost)
	}
	return(x.df)
}

##############

x = costsort.gen(10)
x

