########################
###
### Leonard Mada
### [the one and only]
###
### Tools: Packages & CRAN
### Examples
###
### draft v.0.1a



source("Tools.CRAN.R");


#####################
#####################

################
### Examples ###
################

####################
### Package Size ###
####################

# - Analyse Size of Installed Packages
# Note: takes ages!
if(FALSE) {
	# !! setwd(...); !!
	x = size.pkg();
}
if(FALSE) {
	system.time({
		x = size.pkg(file=NULL);
	})
	# elapsed time: 509 s !!!
	# 512 Packages; 1.64 GB;
}

# using previously saved data:
# - previously saved data;
xsz = read.csv("Packages.Size.csv")

### Size
head(xsz, 20)


###############
###############

###############
### Imports ###
###############

# - much faster: but NO size;
p = info.pkg();
f = imports.pkg();


#####################
### Data Analysis ###


### Description
# - packages which do NOT import any other package;
format.lines(p[is.na(p$Imports), ][1:20, -6])

# - pretty print:
cat.mlines(format.lines(p[is.na(p$Imports), c(1,5,2,3,4)][1:20, ]))

# - pretty print:
scroll.pkg(p[is.na(p$Imports), c(1,5,2,3,4)], start=30)


### Exploratory analysis

# - some are NOT Bioconductor packages;
# - TODO: filter by biocViews?
p[is.na(p$Repository), 1:4]


# No imports
table(is.na(p$Imports))
# Most imported
head(f, 20)

# imported only once:
f$Name[f$Freq == 1]


match.imports("hunspell", p)
match.imports("labeling", p)
match.imports("rpart.plot", p)

# Concept?
match.imports(c("pROC", "ROCR"), p)


###################

### Find in Package

p = info.pkg();

nrow(find.pkg("(?i)matrix", pkg=p))
scroll.pkg(find.pkg("(?i)matrix", pkg=p), start=1)

scroll.pkg(find.pkg("(?i)colou?+r", pkg=p), start=1)

scroll.pkg(find.pkg("(?i)dendro|phylo|tree", pkg=p), start=1)
scroll.pkg(find.pkg("(?i)dendro|phylo", pkg=p), start=1)


#####################
#####################

###################
### Search CRAN ###
###################

### Examples:

### Text-Processing

x = searchCran("text", from=60, sep.h="-")

scroll.pkg(x, start=20, len=21)


### Other packages
# TODO: explore;
# - sources: pubmed.mineR, rplos, rbhl (biodiversity), rcoreoa, biorxivr, pubchunks, jaod,
#   synthesisr, miRetrieve, inpdfr;
# - tools: diffr, cheatR, similr (?), stringdist, sourcetools (src C++), asciiruler;
# - NLP: LDAShiny, corporaexplorer, tokenizers.bpe, wordpiece, text,
#   phm, textmineR, ruimtehol, RKEA, oolong, CSeqpat (?);
# - output: ..., grobblR, REPLesentR, rdoc, GIFTr, formattable, textutils, sassy;
# - other: quanteda.textplots (wordcloud);
# - other search words:
#   mining, language, NLP, LDA, phrase, content, corpora,
#   wordcloud, bibliometric, intrusion;


###############

### Dendrograms
x = searchCran("dendro")

scroll.pkg(x, start=20, len=21)


### PDB
x = searchCran("pdb")

scroll.pkg(x, start=20, len=21)


### Img Processing
x = searchCran("texture")

scroll.pkg(x, start=20, len=21)


### ...
# invasion/invasive, intruder, speciation
x = searchCran("intruder")

scroll.pkg(x, start=20, len=21)


### Percolation
# percol, pore, poros/porou, adsorb
x = searchCran("pore")

scroll.pkg(x, start=20, len=21)

