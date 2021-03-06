#########################
###
### Population Statistics
###
### Leonard Mada


### Populations
# Generate random populations
# - based on the code examples from package MicSim;

### Github:
# https://github.com/discoleo/R/tree/master/Stat/Population.R


#######################

###############
### History ###
###############

### draft v.0.1c - v.0.1c-fix:
# - calculate exact age; [fixed]
### draft v.0.1b:
# - new parameter: proportion males;
### draft v.0.1a:
# - basic functionality;


is.LeapYear = function(y) {
	(y %% 4 == 0) & ((y %% 100 != 0) | (y %% 400 == 0));
}
age.years = function(d, at.date, asFraction=FALSE) {
	# TODO: check thoroughly !!!
	y0   = as.numeric(format(d, format="%Y"));
	y.cr = as.numeric(format(at.date, format="%Y"));
	age  = y.cr - y0 - 1;
	age[age < 0] = 0;
	# Days
	days0 = (as.Date(paste0(y0 + 1, "-01-01")) - d); # days in 1st life-year;
	overhead.days = (d - as.Date(paste0(y0, "-01-01")));
	days.cr = (at.date - as.Date(paste0(y.cr, "-01-01")));
	# print(overhead.days); print(days.cr);
	isLeap0 = is.LeapYear(y0);
	isLeap.cr = is.LeapYear(y.cr);
	dLeap = isLeap0 - isLeap.cr;
	isBirthday = overhead.days <= (days.cr + dLeap);
	age[isBirthday] = age[isBirthday] + 1;
	if(asFraction) {
		# days.startYear = days.cr + dLeap - overhead.days;
		days.startYear = ifelse(isBirthday,
			days.cr + dLeap - overhead.days,
			days0 + days.cr);
		# TODO: Leap Years;
		age = age + days.startYear / (365 + isLeap.cr); # TODO: proper Leap Year;
	}
	return(age);
}

### Test
# TODO: needs a lot of debugging!
age.years(as.Date(c("2020-04-20", "2000-04-20", "2000-04-21", "2000-04-19")),
	as.Date("2021-04-20"), asFraction=T)


rpopulation.gen = function(n, startDate, endDate, simHorizon, format="%d/%m/%Y",
	fertility.lambda=1, collapse=NULL, sex=c("M", "F"), sex.p=0.5,
	marital = c("NM","M","D","W"), edu = c("no","low","med","high")) {
	# Definition of an initial population
	# (for illustration purposes, create a random population)
	if(missing(endDate)) {
		if(length(startDate) == 1) stop("End Date needed!")
		else endDate = startDate[2];
	}
	if(is.character(startDate)) startDate = as.Date(startDate, format=format);
	if(is.character(endDate)) endDate = as.Date(endDate, format=format);
	birthDatesRange = c(startDate[1], endDate[1]);
	
	birthDates <- birthDatesRange[1] +
		runif(N, min=0, max=diff(birthDatesRange))
	
	age <- trunc(as.numeric(simHorizon[1] - birthDates)/365.25)
	s1 = if(length(sex.p) == 1 && sex.p[1] == 1/2) sample(sex, N, replace=TRUE)
		else if(length(sex) - length(sex.p) == 1) sample(sex, N, replace=TRUE, prob=c(sex.p, 1 - sum(sex.p)))
		else sample(sex, N, replace=TRUE, prob=sex.p)
	s2 = rpois(N, ifelse(age <= 18, 0, fertility.lambda)) # Fertility
	# Marital
	s3 = ifelse(age <= 18, marital[1],
			ifelse(age <= 22,
				sample(marital[1:3], N, replace=TRUE),
				sample(marital, N, replace=TRUE)))
	# Education
	s4 = ifelse(age <= 7, edu[1],
			ifelse(age <= 18, edu[2],
			ifelse(age <= 23,
				sample(edu[2:3], N, replace=TRUE),
				sample(edu[-1], N, replace=TRUE))))
	
	if( ! is.null(collapse)) {
		# TODO: can implement option with upper bound on fertility;
		s2ch = as.character(s2);
		s2ch[s2 >= 3] = "3+";
		s = paste(s1,s2ch,s3,s4, sep=collapse)
		initState = data.frame(ID=1:N, birthDate=birthDates, age=age, s);
		return(initState);
	}
	initState = data.frame(ID=1:N, birthDate=birthDates, age=age, s1, s2, s3, s4)
	return(initState)
}

#############

### Examples:

N = 100
dates = c("31/12/1950","01/01/2014")
simHorizon = as.Date("01/01/2021", format="%d/%m/%Y")

initPop = rpopulation.gen(N, dates, simHorizon=simHorizon)
head(initPop)
table(initPop$s3)
mean(initPop$age)
boxplot(initPop$age[initPop$s3 == "W"])


initPop = rpopulation.gen(N, dates, simHorizon=simHorizon, sex.p=0.4)
table(initPop$s1)

initPop = rpopulation.gen(N, dates, simHorizon=simHorizon,
	sex=c("M", "F", "FF"), sex.p=c(0.4, 0.4))
table(initPop$s1)


### with Collapse
initPop = rpopulation.gen(N, dates, simHorizon=simHorizon, collapse="/")
head(initPop)

