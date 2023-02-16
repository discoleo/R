########################
###
### Leonard Mada
### [the one and only]
###
### Data Tools: Date-Time
###
### draft v.0.1c


### Tools to Process/Transform Dates

# - Basic Tools;

######################

### Format Time

# Note: may be useful in:
# - research: biochemistry, cellular biology;
# - clinical trials: pharmacology, hormone assays;

# s.drop = drop seconds;
# collapse = convenience argument;
format.time = function(x, sep=":", h.width=2, s.drop=FALSE, collapse=NULL) {
	h  = x %/% 3600;
	ms = x %%  3600;
	m  = formatC(ms %/% 60, width = 2, flag = 0);
	if(h.width != 1) h = formatC(h, width = h.width, flag = 0);
	if(s.drop) {
		tt = paste(h, m, sep=sep);
	} else {
		s = formatC(ms %% 60, width = 2, flag = 0);
		tt = paste(h, m, s, sep=sep);
	}
	if( ! is.null(collapse)) tt = paste(tt, collapse=collapse);
	return(tt);
}


### Extract Year/Month

extract.YMDate = function(x) {
	xlt = as.POSIXlt(x);
	xmonth = xlt$mon + 1;
	xyear  = xlt$year + 1900;
	rdt = data.frame(year=xyear, month=xmonth);
	return(rdt);
}

is.LeapYear = function(year) {
	return(year %% 4 == 0 && year %% 400 != 0);
}

# Note:
# - benchmarking experiment vs Internal R code;
# - NOT PRODUCTION CODE!
toDate.R = function(s) {
	sf = function(s) {
		bdt = as.integer(charToRaw(s)[-c(5,8)]);
		bdt = bdt - 48L;
		bdt = bdt * c(1000L, 100L, 10L, 1L, 10L, 1L, 10L, 1L);
		year = sum(bdt[1:4]);
		month = bdt[5] + bdt[6];
		day   = bdt[7] + bdt[8] + days[month];
		return(c(year, day));
	}
	dt0 = s[1]; year0 = extract.YMDate(dt0)$year;
	# assuming same year (or only years of same type: Leap *OR* non-Leap)!
	days = as.integer(start.month(year0));
	r = lapply(s, sf);
	r = do.call(rbind, r);
	scale.day = 86400L;
	scale.year = 31556926L;
	r[,2] = r[,2] * scale.day;
	r[,1] = (r[,1] - 1970L) * scale.year + r[,2];
	r = as.POSIXct(r[,1], origin="1970-01-01");
	return(r);
}


##############
### Months ###
##############

### Start of each Month
# - returns also the start of the new year: the 13th month;
start.month = function(year=NULL, day0=0) {
	days = c(day0,31,28, 31,30,31,30,31,31,30,31,30,31);
	if( ! is.null(year) && is.LeapYear(year)) {
		days[3] = 29;
	}
	days = cumsum(days);
	return(days);
}
### Start of each Month
# - after start.dt;
shift.days1 = function(start.dt, day0=0) {
	xlt = as.POSIXlt(start.dt);
	year = xlt$year + 1900;
	diff.days = as.Date(start.dt) - as.Date(trunc.POSIXt(xlt, "years"));
	days = start.month(year=year, day0=day0);
	offset.days = days - diff.days;
	isBefore = (offset.days < 0); # move to next year;
	ydays = if(is.LeapYear(year)) 366 else 365;
	days = c(offset.days[ ! isBefore], ydays + offset.days[isBefore]);
	return(days);
}
shift.days = function(start.dt, day0=0) {
	sapply(start.dt, shift.days1, day0=day0);
}


###################

### Test

shift.days(as.Date(c("2019-01-10", "2020-01-10", "2019-03-24", "2021-06-19")))


# TODO: repeat benchmarks;
toDate.R(c("2019-03-24", "2021-06-19"))


