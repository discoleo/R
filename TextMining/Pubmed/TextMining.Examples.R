

library(xml2)
library(corporaexplorer)

source("Pubmed.XML.R")
source("TextMining.R")

# source("Tools.CRAN.R")
# - needed for formatting;
# - see files in the "/R/Stat" folder;


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

# Visualize
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


###################

### All Parenthesis & Errors
# - output: list with npos;
# - data.frame(nS, nE, ..., Err = 0);

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


### Strip Parenthesis
id = c(3,4)
txt = stripParenth(abstracts$Abstract[id], pos=allParenth[id])

scroll.txt(cbind(abstracts$PMID[id], txt), len=10)


### Contents of Parenthesis

allParenth = parseParenth(abstracts$Abstract);
isErr = isErrParenth(allParenth)
tmp = extractParenthV(abstracts$Abstract[! isErr], allParenth[ ! isErr])
head(tmp)

tmp = unlist(tmp);

# Numeric
isNum = is.string.numExt(tmp)
table(isNum)

head(tmp[isNum], n=20)
# unique(tmp[isNum])
# table(gsub("\\d", "n", tmp[isNum]))
table(encode.num(tmp[isNum]))

# Specific Types
tmp = unlist(tmp);
tmp = unique(tmp);
hasNum = grepl("\\d", tmp);
tmp = tmp[hasNum];
length(tmp)
# Numeric
isNum = is.string.numExt(tmp)
table(isNum)
tmp = tmp[ ! isNum]
# P-Values
isPVal = is.string.pVal(tmp)
table(isPVal)
tmp = tmp[ ! isPVal]
# N-Values
isNVal = is.string.nVal(tmp)
table(isNVal)
tmp = tmp[ ! isNVal]

tmp[1:100]

# also "95 % OR"
regCI = paste0(
	"^95\\%[ ]*+CI", "[\\:, =]++",
	"\\d++(?:[.\uB7]\\d++)?+", "\\%?+", "[- ,;\u2009]++",
	"\\d++(?:[.\uB7]\\d++)?+", "\\%?+$")
isCI = grepl(regCI, tmp, perl=TRUE)
table(isCI)

tmp = tmp[ ! isCI]


##################
##################

### Tokenization & Dependency Parsing

library(udpipe)
library(textplot)

# hacked function:
# source("textplot.Hack.R")

layoutSplit = function(x, y=NULL, nr=2, scale=c(1, 30), dx=3) {
	len = nchar(x);
	if( ! is.null(y)) len = pmax(len, nchar(y));
	len = len + dx;
	if(nr > 1) {
		nm = length(len) %% nr;
		if(nm != 0) len = c(len, rep(0, nr - nm));
	} else nm = 0;
	nc  = length(len) %/% nr;
	len = matrix(len, ncol=nr);
	# x-Coord = cumsum() + mid;
	pos = apply(len, 2, cumsum);
	pos = rbind(0, pos[ - nrow(pos), ]);
	pos = pos + len / 2; # midpoint
	dim(pos) = NULL;
	pos = pos * scale[1];
	xdf = data.frame(
		x = pos,
		y = rep(seq(nr, 1, by=-1) * scale[2], each=nc));
	if(nm > 0) {
		xdf = xdf[ seq(nr - nm) - 1 - nrow(xdf), ];
	}
	return(xdf)
}

# udmodel = udpipe_download_model(language = "english-ewt")
# udmodel = udpipe_load_model(file = udmodel$file_model)
model.file = paste0(getwd(), "/english-ewt-ud-2.5-191206.udpipe")
udmodel = udpipe_load_model(file = model.file)

txtSection = "(?i)^(?:Objectives?|Results|Conclusions?|Introduction)\\:"

### Examples:
# text = uses the abstracts object;
# abstracts = extractAbstract(x)

id = 100
tmp = sub(txtSection, "", abstracts$Abstract[[id]])
ann = udpipe_annotate(udmodel, x = tmp)
ann = as.data.frame(ann, detailed = TRUE)
str(ann)

# Sentence
idS = 1
txtAnn = ann[ann$sentence_id == idS,]
textplot_dependencyparser(txtAnn)
scroll.txt(abstracts, start=id, len=1)


# using a hacked version of textplot_dependencyparser.default
nrows = 3;
textplot_dependencyparser(txtAnn, layout = layoutSplit(txtAnn$token, txtAnn$upos, nr=nrows), nudge_y = -5)

# Hack => layout = layout;
# ggraph::ggraph(g, layout = "linear") +
# ...


### Section Titles
txtSection = extract.regex(abstracts$Abstract, "(?i)^([^:]{1,25}+)\\:", gr=1)
table(txtSection)


# Strip Section title:
tmp = stripSection(abstracts$Abstract)
scroll.txt(cbind(abstracts$PMID, tmp), len=10, start=300)


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
