
### 2.) Primes

all.primes.f = function(lastN) {
	# Init
	p.i = 3:lastN
	p = c(2)
	
	MAX_LEN = sqrt(p.i[[length(p.i)]])
	p.next = p[[1]]
	
	# the sieve
	for(i in 1:MAX_LEN) {
		is.p = ! (p.i %% p.next == 0)
		p.i = p.i[is.p]

		if(length(p.i) <= 0) {
			break;
		}
		p.next = p.i[[1]]
		p = c(p, p.next)
	}

	# all prime numbers
	p = c(head(p, -1), p.i)
	return(p)
}

### Generate ###

p = all.primes.f(1000)

###################

### Twin Primes ###

p.prev = c(0, head(p, -1))
is.tw = (p - p.prev == 2)
is.tw = is.tw | c(tail(is.tw, -1), FALSE)
# ALL primes <= 2 apart
# (includes also 2)
p[is.tw]

# df consumes too much screen (only head())
is.tw = (p - p.prev == 2) & (p != 2)
p.df = data.frame("p"=p[c(tail(is.tw, -1), FALSE)], p.tw=p[is.tw])
head(p.df, n=20)

#######################



p.card = function(p) {
	pc = data.frame(x.p1=0, x.p2=0, p2=0, p3=0, p5=0)
	
	for(i1 in 2:(length(p)-1)) {
		a = 1:(p[i1] - 1)
		a2 = a^2 + 1
		a3 = a^3 + 2*a + 1
		for(i2 in 1:(i1 - 1)) {
			# p2
			tbl2 = table( (a2 - p[i2]) %% p[i1] == 0)
			# p3
			tbl3 = table( (a3 - p[i2]) %% p[i1] == 0)
			if(length(tbl2) == 1) {
				p2 = 0
			} else {
				p2 = tbl2[[2]]
			}
			if(length(tbl3) == 1) {
				p3 = 0
			} else {
				p3 = tbl3[[2]]
			}
			pc = rbind(pc, c(p[i1], p[i2], p2, p3, 0))
		}
	}
	return(pc)
}

pc = p.card(p)

head(pc, 20)


table(pc[,3] > 0)
table(pc[,4] > 0)

is.mult = pc[,3] > 0 & pc[,4] > 0
tbl.pc = table(is.mult)
tbl.pc


head(pc[ is.mult , ], 20)

is.any = pc[,3] > 0 | pc[,4] > 0
head(pc[ ! is.any , ], 20)


# FALSE  TRUE
# p2
#  7226  6636
# p3
#  4500  9362
# mult
#  9418  4444


### a^i + a
# FALSE  TRUE
# p2
#  7060  6802
# p3
#  4532  9330
# mult
#  9270  4592

