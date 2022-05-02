

library(xml2)

source("Pubmed.XML.R")
source("TextMining.R")

# source("Tools.CRAN.R")
# - needed for formatting;


#####################

### Examples

# - using the previously saved xml search results;


###
xn = "_Pubmed_Covid_Review_part.xml"
x = read_xml(xn)

###
abstracts = extractAbstract(x)
scroll.txt(abstracts, len=10)
scroll.txt(abstracts, start=30, len=10)

parenth = parseParenth(abstracts$Abstract[3])
parenth

id = 37
extractParenth(abstracts$Abstract[id])

# may take some time!
allParenth = lapply(seq(nrow(abstracts)), function(id) {
	extractParenth(abstracts$Abstract[id], warn=id);
})
allParenth[(1:10) + 40]

# Errors
id = 87
allParenth[id]
scroll.txt(abstracts, start=id, len=1)

### Errors
allParenth = parseParenth(abstracts$Abstract);
isErr = isErrParenth(allParenth)
sum(isErr)
abstractsErr = abstracts[isErr,]
scroll.txt(abstractsErr, len=10)


### Nested
strNested = extractNested(abstracts$Abstract, pos=allParenth)
length(strNested)
strNested[1:20]

