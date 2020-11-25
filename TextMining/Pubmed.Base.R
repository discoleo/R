##################
###
### Pubmed: Text Mining in R
### Extract Data from XML files
###
### Leonard Mada
###  2020-11-25
### [draft v.0.2]


### History

### draft v.0.2:
# - shortcut functions to extract abstract or titles;
# - correct handling of multiple AbstractText nodes:
#   concatenation of these nodes;
### draft v.0.1:
# - initial cloned version;
# - separated xml-structure code (separate R script);

### Note:
# - modified & improved version;
# - initial R script:
#   https://github.com/discoleo/R/blob/master/XML/R.xml.R


###################

###################
###  Libraries  ###
###################

### 1.) xml2
### 2.) xml: prior to xml2; (not maintained)

### Other packages
### flatxml: may be used by other packages;

### Other uses of xml
### A.) GIS data:
### B.) Sensor data:
### B.1.) Sensor Observation Services: sos4R
###       https://www.opengeospatial.org/standards/sos
### C.) Ontologies: ontologyIndex, rols, ontoCAT
### - 243 ontologies available in the OLS;
### [OWL and OBO formats]

#########################

# install.packages("xml2")
# install.packages("magrittr")
# install.packages("ggraph")

library(xml2)
# uses internally the 'libxml2' C library

# helper
library(magrittr)

###################

### set the Working Directory
setwd("...\\DB")

##################

############
### Data ###

### XML datasets:
# - all datasets are downloaded from Pubmed;

### 200 MB xml
# file.name = "Allergies_Pubmed.part-6.xml"
### 22 MB xml
file.name = "term=(malaria[TIAB])+AND+2017[pdat].xml"

x = read_xml(file.name)

####################

### XPATH Base-Paths

base_path_articles = "/PubmedArticleSet/PubmedArticle"
base_path_medline  = paste(base_path_articles, "/MedlineCitation", sep="")
base_path_article  = paste(base_path_medline, "/Article", sep="")
base_path_journal  = paste(base_path_article, "/Journal", sep="")
base_path_author   = paste(base_path_article, "/AuthorList/Author", sep="")
base_path_text     = paste(base_path_article, "/Abstract/AbstractText", sep="")
base_path_title    = paste(base_path_article, "/ArticleTitle", sep="")

path_text = "/Abstract/AbstractText"

####################

### helper functions

### Filter by PMID
filter.id.xpath = function(pmid, subpath="/Article/Abstract/AbstractText") {
	return(paste0(base_path_medline, "[./PMID/text() =\"", pmid, "\"]", subpath))
}
### Extract Abstract
extract = function(x, path, type=c("Abstract", "Title"), debug=TRUE) {
	### Extract ALL abstracts
	if(missing(path)) {
		type = match.arg(type)
		print(type)
		path = switch(type,
			Abstract = base_path_text,
			Title = base_path_title,
			NA);
		path = paste(path, "/text()", sep="")
	}
	r = xml_find_all(x, path)
	if(debug) {
		print(head(r))
		print(length(r))
	}
	return(r)
}
### Count <AbstractText> nodes
count = function(x) {
	path = paste0(base_path_article, "/Abstract")
	r.all = xml_find_all(x, path)
	print(length(r.all))
	
	### Count <AbstractText> nodes
	# includes also Other Nodes!
	# len = sapply(r.all, function(n) xml_length(n))
	### Alternative:
	# - takes some time!
	len = sapply(r.all, function(x) xml_find_first(x, "count(./AbstractText)"))

	return(len)
}
c.xml = function(x, rep) {
	id = seq(length(rep))
	id.rep = unlist(sapply(id, function(id) rep(id, rep[id])))
	print(length(id.rep))
	tapply(x, id.rep, function(txt) paste(txt, collapse=" "))
}

#################
#################


### "Shortcut"

########################
### Extract Abstract ###
r = extract(x)
length(r)

len = count(x)
length(len)
sum(len)
head(len, n=20)

# nicely concatenated
r = c.xml(r, len)
substr(r.c[1:5], 1, 100)
length(r)


#####################
### Extract Title ###
r = extract(x, type="T")



#################
### Structure ###

src.base = "C:/Users/Leo Mada/Desktop/Practica"
src.path = paste0(src.base, "/MSc/TextMining/Pubmed.Structure.R")
src2.path = paste0(src.base, "/MSc/TextMining/")

### external R script:
# - explore structure of XML file;
# - can be skipped;
SKIP_STRUCTURE = TRUE
if( ! SKIP_STRUCTURE) {
	source(src.path, echo=TRUE, keep.source=TRUE)
}

###########################
###########################

### Extract Information ###

### XPATH: important to select specific nodes;
# xml_find_all(...)

### 1.) Extract Year
### 2.) Extract Article ID: PMID
### 3.) Extract Abstracts
### 4.) Extract Authors
### 5.) Extract Affiliations

########################
########################

########################
### 1.) Extract Year ###
########################

path = paste(base_path_journal, "/JournalIssue/PubDate/Year/text()", sep="")
path
r = xml_find_all(x, path)

head(r)
length(r)
# Result: is a collection/set of xml nodes!
# table(r) # Error

### Extract: as_list()
date = as_list(r)
head(date)

# conversion to linear Vector
date = unlist(date)
date = as.numeric(date)
head(date)
table(date, useNA="ifany")


### Year Errors:
# - Missing Years

### Year Element, but NO Year text
r = xml_find_all(x, paste(base_path_journal, "/JournalIssue/PubDate/Year[not(text())]", sep=""))
head(r)

### NO Year Element
r = xml_find_all(x, paste(base_path_journal, "/JournalIssue[not(./PubDate/Year)]", sep=""))
head(r)
length(r)

### Alternatie Date: MedlineDate

### Only MedlineDate
# /PubDate[not(./Year)]
# /MedlineDate/text()
r = xml_find_all(x, paste(base_path_journal,
		"/JournalIssue/PubDate[not(./Year)]/MedlineDate/text()", sep=""))
head(r)
length(r)

# Extract Year
r.date = regmatches(r, regexpr("^[0-9]+", r))
head(r.date)
r.date = as.numeric(r.date)
table(r.date)


### NO Year
r = xml_find_all(x, paste(base_path_journal,
		"/JournalIssue/PubDate[not(./Year) and not(./MedlineDate)]", sep=""))
head(r)
length(r)
# NO records without Year!

####################
####################

####################
### ALL Articles ###

base_path_article
r = xml_find_all(x, base_path_article)
length(r)


########################
### 2.) Article ID: PMID

path = paste(base_path_medline, "/PMID/text()", sep="")
path
r = xml_find_all(x, path)
head(r)
length(r)


################

################
### 3.) Abstract

### Extract 1 article by PMID
pmid = 28584181
path = filter.id.xpath(pmid)
r = xml_find_all(x, path)
head(r)
sapply(0:10, function(start) substr(r, 80*start, 80*(start+1)))


### Extract ALL abstracts
path = paste(base_path_text, "/text()", sep="")
path
r = xml_find_all(x, path)
head(r)
length(r)

### TODO:
# - clarify differences in length:
#   3330 abstracts vs 1.647 PMIDs;
# - Text Mining on Abstracts;

### Multiple <AbstractText> nodes
path = paste0(base_path_article, "[count(./Abstract/AbstractText) > 1]")
r = xml_find_all(x, path)
head(r)
length(r)


### Count <AbstractText> nodes
path = paste0(base_path_article, "/Abstract")
r.all = xml_find_all(x, path)
length(r.all)

# takes some time!
len = sapply(r.all, function(x) xml_find_first(x, "count(./AbstractText)"))

### old
# includes also other Nodes (see Other Nodes);
len = sapply(r.all, function(n) xml_length(n))
head(len)

# Siblings:
path = paste0(base_path_article, "/Abstract/AbstractText")
r.all = xml_find_all(x, path)
length(r.all)
sibl = sapply(r.all, xml_siblings)
length(sibl)

len = sapply(sibl, length)
len = len + 1
head(len, n=20)


### Other Nodes:
path = paste0(base_path_article, "/Abstract/*[not(self::AbstractText)]")
r = xml_find_all(x, path)
head(r)
length(r)


###############
###############

###############
### 4.) Authors

### Articles with only 1 Author
path = paste(base_path_article,
		"[count(./AuthorList/Author) = 1]/ArticleTitle/text()", sep="")
path
r = xml_find_all(x, path)
head(r)
length(r)

# Number of Authors per Article
path = paste(base_path_article,
		"[count(./AuthorList/Author) > 0]/ArticleTitle/text()", sep="")
path
r.Title = xml_find_all(x, path)

path = paste("count(", base_path_article,
		"[count(./ArticleTitle) > 0]/AuthorList/Author)", sep="")
path # does NOT work


### Authors per Articles
path_article = paste(base_path_article, "[count(./ArticleTitle) > 0]", sep="")
path_article
r = xml_find_all(x, path_article)
length(r)

# Count Authors: does NOT work
# various techniques explored
a.nr = sapply(r, function(x) as_list(xml_find_all(x, "count(./AuthorList/Author)")))

### xml_length()
# Are there any other types of AuthorList Nodes?
path_not_authors = paste(path_article, "/AuthorList[count(./*) > count(./Author)]", sep="")
length(xml_find_all(x, path_not_authors))
#
# path_authors = paste(path_article, "/AuthorList/Author", sep="")
path_authors = paste(path_article, "/AuthorList/Author", sep="")
path_authors
a.r = xml_find_all(x, path_authors)
length(a.r)
a.r = xml_parent(a.r)
length(a.r)
#
a.nr = xml_length(a.r)
head(a.nr)
length(a.nr)

### TODO:
author.count = data.frame("Nr"=a.nr, "title"=unlist(as_list(r.Title)))
head(author.count)


#################
###  Authors  ###
#################

############################
### 4.) Extract Affiliations

# XPATH 2: tokenize() does NOT work
path = paste(base_path_author,
#		"/AffiliationInfo/Affiliation/tokenize(text(), '[\\s,;]+')", sep=""))
		"/AffiliationInfo/Affiliation/text()", sep="")
path
r = xml_find_all(x, path)

head(r)
length(r)

affil.r = unlist(as_list(r))
words.r = strsplit(affil.r, "[\\s,;]+", perl=TRUE)
head(words.r)

### TODO:
