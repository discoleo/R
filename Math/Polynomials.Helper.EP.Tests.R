########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### Elementary Polynomials
### Tests


### Tests for Elementary Polynomials


#######################


### Helper Functions

### Polynomial Tools
source("Polynomials.Helper.EP.R")


### fast load
# source("Polynomials.Helper.EP.Tests.R")


########################
########################

#############
### Tests ###
#############

### TODO:
# - convert lists to data.frames where applicable;


### E Polynomials:
### TODO:
# - thoroughly check!
Epoly.gen(5, 4)

###
x1 = sqrt(2); x2 = sqrt(3); x3 = sqrt(5); x4 = sqrt(7); x5 = sqrt(11);
xx = c(x1, x2, x3, x4, x5);
E2 = eval.pm(perm.poly(5, c(1,1)), xx);
E3 = x1*x2*(x3+x4+x5) + (x1+x2)*x3*(x4+x5) + (x1+x2+x3)*x4*x5;
E4 = prod(xx) * sum(1/xx); E5 = prod(xx); S = sum(xx);
#
eval.pm(perm.poly(5, c(3,3,3)), xx)
eval.pm(Epoly.gen(3, 5, 3), c(S,E2,E3,E4,E5))

### TODO:
# check for larger n:
Epoly.gen(4, 5, 4)
#
eval.pm(perm.poly(5, rep(4, 4)), xx)
eval.pm(Epoly.gen(4, 5, 4), c(S,E2,E3,E4,E5))


###
test.epoly = function(pseq.all) {
	pseq = pseq.all[pseq.all != 0]
	s = sapply(1:10000, function(id) sample(pseq.all, N))
	s = t(s)
	s = unique(s)
	print(nrow(s));
	r.s = sum(sapply(seq(nrow(s)), function(id) prod(xx^s[id,])))
	r.p = eval.pm(Epoly.distinct(pseq, N), c(S, E2, E3, E4, E5))
	data.frame(r=c(r.s, r.p));
}
N = 5;
pseq.all = 0:4;
# pseq.all = c(0, 1,3,4,6); # pseq.all = c(0, 3,3,3,5);
# pseq.all = c(0, 2,2,3,3); # pseq.all = c(0, 2,2,4,4);
# pseq.all = c(0, 3,2,1,1); # pseq.all = c(0, 4,1,3,3);
test.epoly(pseq.all)


################
### Permutations
n = 4
p = prod.perm.poly(n)
p = sort.pm(p, paste0("x", 1:4));
p

print.pm(p)

apply(perm3(4, p=c(3,2,1)), 1, sum)
table(duplicated(perm3(4, p=c(3,2,1))))


###################
###################

### Test Sum & Diff
# Decomposition of:
# x^n + y^n, x^n - y^n

x = sqrt(2:3)

n = 8
diff(x^n)
eval.pm(diff.E2.pm(n), c(sum(x), prod(x), diff(x)))

n = 7
sum(x^n)
eval.pm(sum.E2.pm(n), c(sum(x), prod(x)))

