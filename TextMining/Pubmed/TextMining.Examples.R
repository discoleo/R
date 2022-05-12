

library(xml2)
library(corporaexplorer)

source("Pubmed.XML.R")
source("TextMining.R")

# source("Tools.CRAN.R")
# - needed for formatting;


#####################

### Examples

# - using the previously saved xml search results;


### Data Set
# - see Pubmed.Examples.R;
xn = "_Pubmed_Covid_Review_part.xml"
x = read_xml(xn)


################
### Keywords ###
################

keys = extractKeywords(x)
scroll.txt(keys, len=10)

w = tableWords(keys$Content)
length(w)
head(w, 20)

w2 = splitDF(as.data.frame(w, stringsAsFactors = FALSE), ncol=4)
scroll.txt(w2, w=rep(c(18,8), 4), len=20)


#################
### Abstracts ###
#################

abstracts = extractAbstract(x)

scroll.txt(abstracts, len=10)

scroll.txt(abstracts, start=30, len=10)

### Abstracts: Parenthesis

### Examples:
# - for 1 particular abstract;
id = 3
parenth = parseParenth(abstracts$Abstract[id])
parenth

# - for 1 particular abstract;
id = 37
extractParenth(abstracts$Abstract[id])

# Extract strings:
# - may take some time!
allParenth = lapply(seq(nrow(abstracts)), function(id) {
	extractParenth(abstracts$Abstract[id], warn=id);
})
allParenth[(1:10) + 40]

### Errors
# - particular Error in abstract with id = ...;
id = 87
allParenth[id]
scroll.txt(abstracts, start=id, len=1)


### All Parenthesis & Errors
# - output: list with npos;
allParenth = parseParenth(abstracts$Abstract);
isErr = isErrParenth(allParenth)
sum(isErr)
countErrParenth(allParenth, isErr)

# only abstracts with errors:
abstractsErr = abstracts[isErr,]
scroll.txt(abstractsErr, len=10)

# count(errors in abstract)
# same as countErrParenth(...)
sapply(allParenth[isErr], function(x) {
	sum(x$Err != 0);
})

# hasP = hasParenth(allParenth);
summary.AbstractParenth(abstracts, allParenth)

# Section Methods: contains an error;
scroll.txt(abstracts[abstracts$PMID == 34558742, ], len=10)

absErr = cbind(abstractsErr, countErrParenth(allParenth, isErr))
scroll.txt(absErr, len=10, w=c(8,76,4), start=20)

# code to count the 22 true mismatches
# (based on the example above)
tmp = read.csv("Pubmed.Abstracts.Corrections.Parenthesis.csv")
length(unique(tmp$PMID))


### Nested
strNested = extractNested(abstracts$Abstract, pos=allParenth)
length(strNested)
strNested[1:20]


##################

### Visualization
### using corporaexlorer

id = grep("^Abstract", names(abstracts));
if(length(id) == 1) {
	names(abstracts)[id] = "Text"
}

xc = prepare_data(abstracts,
	grouping_variable="PMID",
	date_based_corpus=FALSE)

explore(xc)
