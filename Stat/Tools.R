
### Various Tools

library(magrittr)
library(ggplot2)

# some of the tools are based on plyr
library(plyr)


### String Tools

seq.character = function(from, to) {
	if(nchar(from[1]) > 1) {
		# convenience functionality;
		if(is.numeric(to)) return(seq.nchar(from, to));
	}
	from = as.integer(charToRaw(from));
	isUpper = if(from <= 90) TRUE else FALSE;
	off = if(isUpper) 64 else 96;
	from = from - off;
	# Length:
	if(is.numeric(to)) {
		if(isUpper) return(LETTERS[seq(from, from + to - 1)]);
		return(letters[seq(from, from + to)]);
	}
	to = as.integer(charToRaw(to)) - off;
	if(isUpper) return(LETTERS[seq(from, to)]);
	return(letters[seq(from, to)]);
}
# simple replacement for paste(from, seq(n), sep="")
seq.nchar = function(from, length.out, start.at=1, sep="") {
	paste(from, seq(start.at, length.out=length.out), sep=sep)
}


########################
### Data Transformations

### Aggregation
# - using native R;
# - using plyr::ddply;


### using native R
countUnique = function(x, cols) {
	# using native R
	x2 = x[, cols]
	x.agg = aggregate(formula(paste0(cols[1], "~.")), x2, FUN=function(x) length(unique(x)))
	names(x.agg)[length(cols)] = "Count"
	return(x.agg)
}


### using plyr::ddply
ddcountUnique = function(x, cols) {
	x.agg = x %>% select(cols)
	x.agg = eval(substitute(
		ddply(x.agg, y, summarise,
		Count = length(unique(get(z)))), list(y=cols[-1], z=cols[1]) ),
		envir=x.agg, enclos=parent.frame())
	return(x.agg)
}


### Test

# x = some data set
id.name = c("ID1", "ID2")
col.names = c(id.name, "luna")

### native R:
x2 = countUnique(x, col.names)

x2 %>% 
  ggplot(aes(x=luna, y=Count))+
  geom_line()+
  facet_wrap(~ID_Judet)
###


### plyr version: using ddply;
x %>% ddcountUnique(col.names) %>%
  ggplot(aes(x=luna, y=Count))+
  geom_line()+
  facet_wrap(~ID_Judet)


