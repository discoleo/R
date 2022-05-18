

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
##################

### Tokenization & Dependency Parsing

library(udpipe)
library(textplot)

# udmodel = udpipe_download_model(language = "english-ewt")
# udmodel = udpipe_load_model(file = udmodel$file_model)
model.file = paste0(getwd(), "/english-ewt-ud-2.5-191206.udpipe")
udmodel = udpipe_load_model(file = model.file)

# text = uses the abstracts object;
# abstracts = extractAbstract(x)

txtSection = "(?i)^(?:Objectives|Results|Conclusions)\\:"

id = 100
tmp = sub(txtSection, "", abstracts$Abstract[[id]])
ann = udpipe_annotate(udmodel, x = tmp)
ann = as.data.frame(ann, detailed = TRUE)
str(ann)

# Sentence
idS  = 2
txtS = ann[ann$sentence_id == idS,]
textplot_dependencyparser(txtS)
scroll.txt(abstracts, start=id, len=1)


layoutSplit = function(x, y=NULL, nr=2, scale=c(1,20), dx=3) {
	len = nchar(x);
	if( ! is.null(y)) len = pmax(len, nchar(y));
	len = len + dx;
	if(nr > 1) {
		nm = length(len) %% nr;
		if(nm != 0) len = c(len, rep(0, nr - nm));
	} else nm = 0;
	nc  = length(len) %/% nr;
	len = matrix(len, ncol=nr);
	# x-Start = 1;
	# len = rbind(1, len[ - nc, ]);
	len = apply(len, 2, cumsum);
	dim(len) = NULL;
	len = len * scale[1];
	xdf = data.frame(
		x = len,
		y = rep(seq(nr, 1, by=-1) * scale[2], each=nc));
	if(nm > 0) {
		xdf = xdf[ seq(nr - nm) - 1 - nrow(xdf), ];
	}
	return(xdf)
}

# using a hacked version of textplot_dependencyparser.default
textplot_dependencyparser(txtS, layout = layoutSplit(txtS$token, txtS$upos, nr=2), nudge_y = -5)

# => layout = layout;
# ggraph::ggraph(g, layout = "linear") +
# ...


txtSection = extract.regex(abstracts$Abstract, "(?i)^([^:]{1,25}+)\\:", gr=1)
table(txtSection)


##################

### Visualization
### using corporaexlorer

# text = uses the abstracts object;
# abstracts = extractAbstract(x)

# convert col-name to "Text"
id = grep("^Abstract", names(abstracts));
if(length(id) == 1) {
	names(abstracts)[id] = "Text"
}

xc = prepare_data(abstracts,
	grouping_variable="PMID",
	date_based_corpus=FALSE)

explore(xc)


### CSV-File
x = read.csv("Example_Abstracts_Title_Pubmed.csv")

x = x[ -3]; # remove Year column
names(x)[2] = "Text"
xc = prepare_data(x,
	grouping_variable="PMID",
	date_based_corpus=FALSE)

explore(xc)
