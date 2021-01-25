###################
###
### Big Data
### Process Excel files
###
### Leonard Mada
### v 0.2a


# install.packages("readxl")
# install.packages("xts")

library(readxl)



setwd(".../Data/historical-files-up-to-oct-1-included")


####################

### helper functions

get_sheet = function(file.name, sheet.pattern="^[^ _]+(?=[ _])", DEBUG=FALSE) {
	sheets = excel_sheets(file.name)
	if(length(sheets) == 1) {
		return(sheets)
	} else {
		name.reg = regexec(sheet.pattern, file.name, perl=T)
		sName = regmatches(file.name, name.reg)[[1]]
		isSheet = grepl(paste0("^", sName), sheets)
		if(DEBUG) {
			len = length(isSheet[isSheet])
			if(len > 1) {
				print(sheets[isSheet])
			} else if(len == 0) {
				print(paste0("NO sheets: ", sheets))
			}
		}
		sheet = sheets[isSheet][[1]]
		return(sheet)
	}
}
get_sheet.old = function(id, files, sheet.pattern="^[^ _]+(?=[ _])", DEBUG=FALSE) {
	sheets = excel_sheets(files.n[id])
	if(length(sheets) == 1) {
		return(sheets)
	} else {
		name.reg = regexec(sheet.pattern, files.n[id], perl=T)
		sName = regmatches(files.n[id], name.reg)[[1]]
		isSheet = grepl(paste0("^", sName), sheets)
		if(DEBUG && length(isSheet[isSheet]) > 1) {
			print(sheets[isSheet])
		}
		sheet = sheets[isSheet][[1]]
		return(sheet)
	}
}
find_last = function(chars, str) {
	# last occurance
	has = grepl(chars, str)
	len = length(has[has])
	if(len == 1) {
		return(str[has])
	} else {
		return(str[has][len])
	}
}
readNA.excel = function(files, col.name="KWh", sheet.pattern="^[^ _]+(?=[ _.])", PRINT_TABLE = TRUE) {
	# PRINT_TABLE: for debug purposes;
	x.df = data.frame(
		id=seq(length(files.n)),
		na=NA, Zero=NA, name=files.n, sheet=NA)
	col.name.pattern = paste0("^", col.name)
	last.n = nrow(x.df) + 1
	#
	for(first.pass in 1:2) {
		if(first.pass == 1) {
			all.ids = x.df$id
		} else if(nrow(x.df) >= last.n) {
			all.ids = x.df$id[last.n:nrow(x.df)]
		} else {
			break;
		}
	for(id in all.ids) {
		# Sheet name
		sheets = excel_sheets(files.n[id])
		if(length(sheets) == 1) {
			x = read_excel(files.n[id])
		} else {
			# sheet name extracted from filename;
			if( ! is.na(x.df$sheet[id])) {
				sheet = x.df$sheet[id]
				print("Processing multiple sheets!")
			} else {
				name.reg = regexec(sheet.pattern, files.n[id], perl=T)
				sName = regmatches(files.n[id], name.reg)[[1]]
				isSheet = grepl(paste0("^", sName), sheets)
				# if > 1 sheets
				sheets.all = sheets[isSheet]
				sheet = sheets.all[[1]]
				x.df$sheet[id] = sheet;
				len = length(sheets.all) - 1
				if(len >= 1) {
					x2.df = data.frame(
						id = nrow(x.df) + seq(len),
						na=NA, Zero=NA, name=files.n[id], sheet=sheets.all[-1])
					x.df = rbind(x.df, x2.df)
					cat("\n"); print(x2.df)
					files.n = c(files.n, rep(files.n[id], len))
				}
			}
			print(paste0("Opening: ID = ", id, "; Sheet = ", sheet))
			x = read_excel(files.n[id], sheet)
		}
		# Col names
		if( ! (col.name %in% names(x))) {
			print(paste0("Missing column: ID = ", id, "; Sheet = ", sheet))
			col.names = names(x)
			kwh.col = find_last(col.name.pattern, col.names)
			name.pos = match(kwh.col, col.names)
			names(x)[name.pos] = col.name
			if(PRINT_TABLE) {
				# Debug!
				print(table(x$KWh > 0))
			}
		}
		# x[ x[ , col.name] == 0 , col.name] = NA
		### NA
		isNA = is.na(x[,col.name])
		count = if(length(isNA) == 1) 0 else sum(isNA)
		x.df$na[id] = count
		### Zero
		isZero = x[ , col.name] == 0
		x.df$Zero[id] = length(isZero[isZero])
	}}
	return(x.df)
}

#################

#################
### All files ###

flt = "\\.xlsx$"

# All files in directory
files.n = list.files(pattern = flt)


### Read from ALL files

# some files contain 2 sheets with data;
# sapply(files.n, get_sheet, DEBUG=T)
# TODO: not yet extracted!

# How many NAs & Zero values:
x.df = readNA.excel(files.n)

# x.df
x.df = x.df[order(x.df$name),]
x.df

#################

#################
### Read data ###

# read a specific data set;

# file.n = "ALM KWh data.xlsx"

id = 3
file.n = files.n[id]

x = read_excel(file.n, get_sheet(files.n[id]))
# x = read_excel(file.n)
# x = x[,1:4]
head(x)

### NA: 1st value in series
table(is.na(x$KWh))
x = x[ ! is.na(x$KWh), ]

################
################

################
### TS

library(xts)

### whole Days
day = 4*24
nrow(x) %% day
# missing time points
day - (nrow(x) %% day)
### TODO: imputation


x.ts <- xts(x$KWh, order.by=x$timeLogged)
attr(x.ts, 'frequency') = 4*24
names(x.ts) = "KWh"
x.ts = x.ts[ ! is.na(x.ts) ]
head(x.ts)

plot(x.ts)


simple.ret <- function(x, col.name, k=1){
  x[,col.name] / lag(x[,col.name], k=k) - 1
}
plot(x.ts, FUN=simple.ret, col.name="KWh", k=4*24) # NO effect


### Filter
filter.val = function(x, limit=c(0, 200), frequency=4*24) {
	# daily = 4*24
	# 4 times per day = 4*24 / 4
	# weekly = 4*24*7
	attr(x, 'frequency') = frequency;
	x$KWh[x$KWh > limit[2]] = 0
	x$KWh[x$KWh < limit[1]] = 0
	return(x)
}
table(x.ts$KWh > 200)
table(x.ts$KWh > 100)
table(x.ts$KWh < 0)


### Decompose
plot.decomposed.ts = function (x, start=1, len=end(x.decomp$x), ...) {
	end = start + len - 1
	flt = function(x) {
		window(x, start=start, end=end)
	}
    xx <- flt(x$x)
    if (is.null(xx)) 
        xx <- with(x, if (type == "additive") 
            flt(random) + flt(trend) + flt(seasonal)
        else flt(random) * flt(trend) * flt(seasonal))
    plot(cbind(observed = xx, trend = flt(x$trend), seasonal = flt(x$seasonal), random = flt(x$random)), 
         main = paste("Decomposition of", x$type, "time series"), ...)
}

### Decompose
# - but only 1 seasonality!
freq = 4*24;
freq = 4*24*7;
x2.ts = filter.val(x.ts, frequency=freq)
x.decomp <- decompose(as.ts(x2.ts), type="additive")

# daily
plot(x.decomp, start=50, len=90)
# weekly
plot(x.decomp, start=3, len=7)


### Multiple Seasonalities
library(forecast)

x.mts = msts(x2.ts, seasonal.periods=c(4*24, 4*24*7))

x.fit <- tbats(x.mts)
plot(forecast(x.fit))


hw.fit = HoltWinters(x.mts)
hw.fcst = forecast(hw.fit, h=8)

# analyze the autocorrelation of the residuals
# (a good model should lead to no autocorrelation in residuals)
acf(hw.fcst, lag.max=20)

# a good model should also have the forecast errors normally distributed
# with mean zero and constant variance
plot.ts(hw.fcst$residuals)


### Day of week: simple Mean
ts.m = matrix(x2.ts, nrow=4*24*7)
ts.mm = apply(ts.m, 1, mean)
plot(as.ts(ts.mm))
plot(ts.mm) # low consumption during certain hours


### TODO:
# fasster: Forecasting multiple seasonality with state switching
# https://www.youtube.com/watch?v=6YlboftSalY


################
### Analyse data

table(is.na(x$KWh))
x[is.na(x$KWh),]

# extreme outliers
boxplot(x$KWh)

isOutlier = x$KWh > 10
table(isOutlier)
x[isOutlier, ]

isOutlier = x$KWh <= 0
table(isOutlier)
x[isOutlier, ]


##################
### Partition Data

### TODO
med = quantile(x$KWh, 0.5, na.rm=T)
x.norm = (x$KWh - med)
x.norm[is.na(x.norm)] = -1000

x.cl = kmeans(x.norm, 4)
table(x.cl$cluster)

