

library(xml2)

source("Pubmed.XML.R")

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
	extractParenth(abstracts$Abstract[id]);
})
allParenth[(1:10) + 40]

