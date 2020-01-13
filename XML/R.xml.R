
### XML Data in R
###
### Leonard Mada
### [draft 0.1]
###  2020-01-13

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
# visualization of trees/graphs
library(ggraph)
library(tidygraph)
library(igraph)

###################

### set the Working Directory
setwd("...") # !!!


### 20 - 200 MB xml
# 22 MB xml
# ArticleSet downloaded from Pubmed:
x = read_xml("term=(malaria[TIAB])+AND+2017[pdat].xml")

### XPATH Base-Paths

base_path_articles = "/PubmedArticleSet/PubmedArticle"
base_path_article  = paste(base_path_articles, "/MedlineCitation/Article", sep="")
base_path_journal  = paste(base_path_article,  "/Journal", sep="")
base_path_author   = paste(base_path_article,  "/AuthorList/Author", sep="")


################

### Structure
# xml_structure(...)

# select only 1 Article: position() <= 1
q.xpath = paste(base_path_articles, "[position() <= 1]", sep="")
q.xpath
x.small = xml_find_all(x, q.xpath)

# display the structure:
# unfortunately, it is only printed!
xml_structure(x.small, indent = 2)

### a more detailed Structure:
# xml_children(): does NOT extract all descendants!
# xml_children(x.small) # only the 2 children;
### extract all nodes: "/descendant::*";
q.xpath = paste(base_path_articles, "[position() <= 1]/descendant::*", sep="")
q.xpath
x.nodes = xml_find_all(x, q.xpath)
x.nodes

# XPATH: "name()": does not work directly;
xml_name(x.nodes)

# xml_path(): extract the full path;
# shortcut: "//*" == "/descendant::*"
q.xpath = paste(base_path_articles, "[position() <= 1]//*", sep="")
q.xpath
strip.len = nchar(base_path_articles) + 1
# extract the full path:
x.small %>% xml_find_all( q.xpath ) %>%
	xml_path() %>% substr(strip.len, 0xFFFF)

### Create Graph:
# - extract individual Nodes from Path;
# - group nodes to create the edges;
gr.split = function(s) {
	s.split = unlist(strsplit(s, "/"))
	len = length(s.split)
	# Edges
	pair = sapply(2:len, function(npos) c(s.split[npos-1], s.split[npos]))
	return(t(pair))
}
gr.vsplit = function(s) {
	sapply(s, function(s) gr.split(unlist(s)) )
}
nodes.rez = x.small %>% xml_find_all( q.xpath) %>%
	xml_path() %>%
	substr(strip.len, 0xFFFF) %>%
	gr.vsplit()
head(nodes.rez)

# add all edges to 1 Matrix:
# [code may not be particularly efficient!]
m = matrix(character(0), ncol=2)
tmp = lapply(nodes.rez, function(l) {m <<- rbind(m, l); return();})
# probably the correct alternative:
# m = matrix(unlist(nodes.rez), ncol=2) # TODO: check if it works correct!
colnames(m) = c("from", "to")
head(m)

### group certain Vertices together
f.group_by = function(v, edges, group.tag, gr=1, group) {
	group[match(group.tag, v$name)] = gr
	MAX_DEPTH = 10
	for(i in 1:MAX_DEPTH) {
		# TODO: process only newly selected nodes;
		L2 = unique(edges$to[edges$from %in% v$id[group == gr]])
		group[v$id %in% L2] = gr
	}
	return(group)
}

# ALL vertices
v = unique(matrix(m, ncol=1))
v.names = data.frame(name = v, stringsAsFactors=FALSE, id=factor(v))
# ALL edges
v.edges = data.frame(from=factor(m[,1], levels=v), to=factor(m[,2], levels=v))
head(v.edges)

# Group vertices together
v.names$group = 0
group.tag = "AuthorList" # AuthorList
v.names$group = f.group_by(v.names, v.edges, group.tag, 1, v.names$group)
group.tag = "Journal" # Journal
v.names$group = f.group_by(v.names, v.edges, group.tag, 2, v.names$group)
v.names[ v.names$group > 0 , ]

# there is a BUG in tbl_graph:
# geom_node_label() will associate the labels wrongly!
# gr = tbl_graph(nodes=v.names, edges=v.edges, directed = TRUE)
gr <- graph_from_data_frame( v.edges, vertices=v.names )
head(gr)

### plot ggraph
ggraph(gr, layout = 'lgl', circular = FALSE) +
	# layout = 'dendrogram': some nodes are NOT unique!
	### Edge type
	# geom_edge_diagonal() +
	geom_edge_link() +
	### Labels
	geom_node_label(aes(x = x, y=y, label=name, color=group)) +
	# geom_node_point() +
	theme_void()


### igraph: too much overlap!
g.igr = graph_from_data_frame(v.edges, directed = FALSE, vertices =v.names$v)
plot(g.igr)

### interactive plot!
tkplot(gr, canvas.width = 550, vertex.color="lightblue")

###########################
###########################

### Extract Information ###

### XPATH: important to select specific nodes;
# xml_find_all(...)

### 1.) Extract Year
### 2.) Extract Article ID: PMID
### 3.) Extract Authors
### 4.) Extract Affiliations

####################
### 1.) Extract Year
path = paste(base_path_journal, "/JournalIssue/PubDate/Year/text()", sep="")
path
r = xml_find_all(x, path)

head(r)
length(r)
# result: are still nodes in the xml
table(r) # Error

### Extract: as_list()
date = as_list(r)
head(date)

# conversion to linear Vector
date = unlist(date)
date = as.numeric(date)
head(date)
table(date, useNA="ifany")


### Year Errors: Missing Years

### Year Element, but NO Year text
r = xml_find_all(x, paste(base_path_journal, "/JournalIssue/PubDate/Year[not(text())]", sep=""))
head(r)

### MedlineDate
r = xml_find_all(x, paste(base_path_journal, "/JournalIssue[not(./PubDate/Year)]", sep=""))
head(r)
length(r)

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


####################
### ALL Articles ###

base_path_article
r = xml_find_all(x, base_path_article)
length(r)


################################
### 2.) Extract Article ID: PMID

path = paste(base_path_articles,
		"/MedlineCitation/PMID/text()", sep="")
path
r = xml_find_all(x, path)
head(r)


#######################
### 3.) Extract Authors

# Articles with only 1 Author
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
