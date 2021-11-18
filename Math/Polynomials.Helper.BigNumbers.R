########################
###
### Leonard Mada
### [the one and only]
###
### Helper Functions
### BigNumbers


#######################

### fast load:
# source("Polynomials.Helper.BigNumbers.R")

### Helper Functions

# already loaded!
# library(polynom)
# library(pracma)

library(gmp) # BigNumbers

### helper Functions

# - hack to work with BigNumbers;
# - workaround for: aggregate();

aggregate0.pm = function(p) {
	# hack for bigz numbers
	coeff = p$coeff;
	p$coeff = seq(nrow(p));
	aggr.bignum = function(id) {
		x.df = data.frame(coeff=0);
		s = sum(coeff[id]);
		x.df$coeff = s;
	}
	p.r = aggregate(coeff~., p, aggr.bignum);
	return(p.r);
}

#############

