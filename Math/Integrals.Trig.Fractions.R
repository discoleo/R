


### Integrals: Trig Fractions

### Note:
# - Fractions of type: 1 / (sin(x)^n + cos(x)^m)
#   are in file: Integrals.Fractions.Trig.R;


####################

### Helper Functions

### I( 1 / (sin(x) + a) )
int.sinfr = function(x, a) {
	if(a == 1) {
		r = 2*sin(x/2) / (sin(x/2) + cos(x/2));
		return(diff(r));
	}
	a.sq = sqrt(as.complex(1 - a^2));
	# with fractions; alternative: atan;
	r1 = log((a.sq - cos(x))/(a.sq + cos(x)));
	# use result from 1/(sin(x)^2 - a^2);
	r2 = log((a.sq * tan(x) - a)/(a.sq * tan(x) + a));
	r = diff(r1 - r2);
	if(Im(a) != 0 ) {
		sg = if(Im(a) < 0) 1 else -1;
		r = r + 2i*sg*pi;
	}
	r = r / (2*a.sq);
	return(r)
}


########################

### I( 1 / (sin(x)^n + 1) )

### n = 3
n = 3;
cs = cos(2*pi/n); sn = sin(2*pi/n);
integrate(\(x) 1 / (sin(x)^n + 1), 0, pi/2)
(int.sinfr(c(0, pi/2), 1) +
	+ (cs + 1i*sn) * int.sinfr(c(0, pi/2), cs + 1i*sn) +
	+ (cs - 1i*sn) * int.sinfr(c(0, pi/2), cs - 1i*sn)) / n;
# TODO: compute expression;


# Fraction Decomposition:
x = pi/7; # Test
n = 3;
cs = cos(2*pi/n); sn = sin(2*pi/n);
#
1 / (sin(x)^n + 1)
1/n / (sin(x) + 1) + 2/n * (cs*sin(x) + 1) / ((sin(x) + cs)^2 + sn^2)
#
1/n / (sin(x) + 1) +
	+ 1/n * (cs - 1i*sn) / (sin(x) + cs - 1i*sn) +
	+ 1/n * (cs + 1i*sn) / (sin(x) + cs + 1i*sn)


#########
### n = 5
n = 5;
integrate(\(x) 1 / (sin(x)^n + 1), 0, pi/2)
# TODO


# Fraction Decomposition:
x = pi/7; # Test
n = 5; id = c(1,2);
cs = cos(2*id*pi/n); sn = sin(2*id*pi/n);
#
1 / (sin(x)^n + 1)
1/n / (sin(x) + 1) + 2/n * sum( (cs*sin(x) + 1) / ((sin(x) + cs)^2 + sn^2) )

