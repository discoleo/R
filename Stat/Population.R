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

### draft v.0.1b:
# - new parameter: proportion males;
### draft v.0.1a:
# - basic functionality;


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
		else if(length(sex.p) == 1) sample(sex, N, replace=TRUE, prob=c(sex.p, 1 - sex.p))
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


### with Collapse
initPop = rpopulation.gen(N, dates, simHorizon=simHorizon, collapse="/")
head(initPop)

