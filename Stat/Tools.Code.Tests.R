

source("Tools.Code.R")


###########

### Example

pkg = "partitions"
ls(getNamespace(pkg))

###
extract.str.pkg(pkg)

###
extract.str.pkg("sp") # relatively fast
# 1st run: ages! afterwards: 0.2 s;
extract.str.pkg("NetLogoR")

### All strings in "base"
extract.str.pkg("base")

###
pkg = "partitions"
fn = "restrictedparts"
extract.str.fun(fn, pkg)

### C-Code:
# Type == 50
s = deparse.fun("C_acf", "stats");
npos = parse.RC(s)
# npos = cut.code(npos)
# head(npos)
# preserve NL;
cat(extract.str(s, npos[npos$Type == 50, ]), sep="\n")

###
pkg = "stats"
nms = ls(getNamespace(pkg));
# takes ~ 10-20 s;
ids = sapply(nms, function(FUN.nm) {
	s = deparse.fun(FUN.nm, pkg);
	isC = is.code.RC(s);
	return(isC);
})
print(ids[ids])


###################

### All code tokens

### Ex 1:
extract.str.fun("compositions", "partitions", type=99)

### Ex 2:
extract.str.fun("summary.Spatial", "sp", type=99)

### Ex 3:
extract.str.fun("StrSpell", "DescTools", type=99)
#
cat(extract.str.fun("StrTrim", "DescTools", type=99, trim.regex=regex.trim()), sep="\n")
cat(extract.str.fun("StrTrim", "DescTools", type=24, trim.regex=regex.trim()), sep="\n\n")

### Ex 4:
s = deparse.fun("PlotMosaic", "DescTools")
npos = parse.fun("PlotMosaic", "DescTools")
head(npos)
npos = cut.code(npos)
head(npos)
# preserve NL;
cat(extract.str(s, npos, trim.regex=regex.trim()), sep="\n")


### TODO:
# - functions to filter the tokens:
#   e.g. simple calls "()", etc;
# - functions to extract the names of the called functions;


####################

# setwd(r'(...)')

file = "mixregEngine.R"
s = read.code(file=file)
s = format.code(s)
#
file.out = file(sub("\\.R$", ".tmp.R", file))
writeLines(paste0(s, collapse=""), file.out)
close(file.out)

####################

### Other
paste0(deparse(partitions::restrictedparts), collapse=" ");


####################

#############
### Tests ###
#############

s = "# Comment"
nchar(s)
parse.simple(s)


s = "# Comment\n"
nchar(s)
parse.simple(s)


s = r"(# Comment
# Comment 2
s = "String")"
nchar(s)
npos = parse.simple(s)
npos
cat(extract.str(s, npos), sep="\n")


s = "\"this is \\\" still a String\"
# Comment with \"pseudo-string
s = \"String \\n \\\"2\\\"\""
nchar(s)
npos = parse.simple(s)
npos
cat(extract.str(s, npos), sep="\n")


gsub("[\\s]++", "", "s  \t\nsbc\nbb\uA0UU", perl=T)


######################
######################

### Function Arguments

e = formals(stats:::plot.lm)
a = summary.args(e)
aggregate(rep(1, nrow(a)) ~ type, a, FUN=length)


f = ls(getNamespace("partitions"))
e = parse(text=paste0("formals(partitions:::", f[39], ")"))
e = eval(e)
a = summary.args(e)
aggregate(rep(1, nrow(a)) ~ type, a, FUN=length)


summary.all.args("partitions")

