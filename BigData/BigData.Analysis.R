###################
###
### Big Data
### Process Excel files
###
### Leonard Mada
### v 0.10


# install.packages("readxl")

library(readxl)



setwd("...\\BigData\\Data\\historical-files-up-to-oct-1-included")


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
		na=NA, Zero=NA, name=files.n)
	col.name.pattern = paste0("^", col.name)
	#
	for(id in x.df$id) {
		# Sheet name
		sheets = excel_sheets(files.n[id])
		if(length(sheets) == 1) {
			x = read_excel(files.n[id])
		} else {
			# sheet name extracted from filename;
			name.reg = regexec(sheet.pattern, files.n[id], perl=T)
			sName = regmatches(files.n[id], name.reg)[[1]]
			isSheet = grepl(paste0("^", sName), sheets)
			sheet = sheets[isSheet][[1]]
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
		# NA
		isNA = is.na(x[,col.name])
		count = if(length(isNA) == 1) 0 else sum(isNA)
		x.df$na[id] = count
		# Zero
		isZero = x[ , col.name] == 0
		x.df$Zero[id] = length(isZero[isZero])
	}
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
sapply(files.n, get_sheet, DEBUG=T)
# TODO: not yet extracted!

# How many NAs & Zero values:
x.df = readNA.excel(files)

x.df


#############
### Read data

# file.n = "ALM KWh data.xlsx"

id = 3
file.n = files.n[id]

x = read_excel(file.n, get_sheet(files.n[id]))
# x = read_excel(file.n)


################
################

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

