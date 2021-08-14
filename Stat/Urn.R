
### Various Urn Processes


###################

shift = function(x, by=1) {
	x = c(tail(x, by), head(x, -by));
	return(x);
}

############
### Urns ###
############

### Non-random Urns
urn.simple = function(n, balls.cat) {
	# Urns: each containing balls of same type;
	# balls.cat = types of balls;
	len = length(n);
	if(len == 1) {len = len + 1; n = c(n,n);}
	if(missing(balls.cat)) {balls.cat = seq(0, len-1);}
	else if(length(balls.cat) > len) warning("More categories than urns!");
	b.unqcat = sort(unique(balls.cat)); # unique categories;
	zero = rep(0, length(b.unqcat));
	#
	urn = array(numeric(0), c(length(b.unqcat), len));
	for(i in seq(len)) {
		balls.count = zero;
		balls.count[which(balls.cat[i] == b.unqcat)] = n[i];
		urn[,i] = balls.count;
	}
	urn = t(urn);
	class(urn) = c("urn", class(urn));
	attr(urn, "cat") = b.unqcat;
	attr(urn, "plen") = length(b.unqcat);
	return(urn)
}

#################
### Explicit Urns

# - specific realization with balls; 

### Random Urns

### TODO: p = matrix;
rurn = function(n, p=1/2) {
	# Equal Urns
	# length(p) = categories of balls;
	len = length(n);
	if(len == 1) {len = 2; n = rep(n, len);}
	plen = length(p);
	add.cat = FALSE;
	if(plen == 1) {add.cat = TRUE;}
	else if(round(sum(p), 5) < 1) {add.cat = TRUE;}
	### n Balls
	ngr = list();
	for(i in seq_along(n)) {
		if(add.cat) {
			tmp = round(n[i] * p);
		} else {
			tmp = round(n[i] * p[-plen]);
		}
		ngr[[i]] = c(tmp, n[i] - sum(tmp));
	}
	if(add.cat) plen = plen + 1;
	# Balls by Category
	urn = list();
	for(i in seq(length(n))) {
		balls = rep(seq(0, plen-1), ngr[[i]]);
		urn[[i]] = sample(balls, n[i]);
	}
	class(urn) = c("urnb", class(urn));
	attr(urn, "plen") = plen;
	return(urn)
}
### Non-random Urns
rurn.simple = function(n, nbc) {
	# Urns: each containing balls of same type;
	# nbc = count with the types of balls;
	len = length(n);
	if(len == 1) {len = len + 1; n = c(n,n);}
	if(missing(nbc)) nbc = seq(0, len-1);
	urn = list();
	for(i in seq(len)) {
		urn[[i]] = rep(nbc[i], n[i]);
	}
	class(urn) = c("urnb", class(urn));
	attr(urn, "plen") = len;
	return(urn)
}

### Urn Processes

swap.urnb = function(urn, iter=10) {
	nall = length(urn);
	n = sapply(urn, length);
	#
	ids = sapply(seq(nall), function(idu) sample(seq(n[idu]), iter, replace=TRUE));
	ids = t(ids);
	for(i in seq(iter)) {
		balls = sapply(seq(nall), function(idu) urn[[idu]][ids[idu, i]]);
		# swap the balls
		balls = shift(balls, by=1);
		# insert balls back
		for(idu in seq(nall)) {
			urn[[idu]][ids[idu, i]] = balls[idu];
		}
	}
	return(urn);
}


################
################

### Count-Urns:
urn = urn.simple(c(20,10,10), c(3,3,1))
urn

#################

### [old]
# - explicit urns;

### Ex: Eq urns
urn = rurn(10)
sapply(urn, mean)

urn2 = swap.urnb(urn)
sapply(urn2, mean)


### Ex: non-Eq urns
urn = rurn(c(20, 10))
sapply(urn, mean)

urn2 = swap.urnb(urn)
sapply(urn2, mean)


### Ex: non-Eq urns
urn = rurn(c(30, 20, 10))
sapply(urn, mean)

urn2 = swap.urnb(urn)
sapply(urn2, mean)


### 1-Coloured Urns
urn = rurn.simple(c(30, 20, 10))
sapply(urn, mean)

urn2 = swap.urnb(urn, iter=100)
sapply(urn2, mean)

### TODO:
# - Mixing times;

