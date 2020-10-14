
##################
### XML Data in R:
### XML Structure
###
### Leonard Mada
###  2020-10-15
### [draft v.0.1]

### Note:
# - modified & improved version;
# - initial R script:
#   https://github.com/discoleo/R/blob/master/XML/R.xml.R


###################
###  Libraries  ###
###################

### 1.) xml2
### 2.) xml: prior to xml2; (not maintained)

### Other packages
### flatxml: may be used by other packages;


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

### Data
# - the data is imported outside this "package";
# - variables:
#   x = xml data;

###################

### XPATH Base-Paths
# - imported outside this "package":
#   base_path_... = the various paths;
# - see file: Pubmed.Base.R;

###################
###################

### helper Functions
limit.xpath = function(base_path, n=1) {
	paste0(base_path, "[position() <= ", n, "]//*")
}

#################


#################
### Structure ###

### very simple structure:
# Function: xml_structure(...)

# select only 1 Article: position() <= 1
q.xpath = paste(base_path_articles, "[position() <= 1]", sep="")
q.xpath
x.small = xml_find_all(x, q.xpath)

# display the structure:
# unfortunately, it is only printed!
xml_structure(x.small, indent = 2)


### Detailed Structure:
# - code needed to extract a more detailed Structure:

### XML Functions:
# xml_children(): does NOT extract all descendants!
# xml_children(x.small) # only the 2 children;
### extract all nodes: "/descendant::*";
q.xpath = paste(base_path_articles, "[position() <= 1]/descendant::*", sep="")
q.xpath
x.nodes = xml_find_all(x, q.xpath)
x.nodes # Nodes + Text

# XPATH: "name()": does not work directly;
xml_name(x.nodes)


# xml_path(): extract the full path;
# shortcut: "//*" == "/descendant::*"
q.xpath = paste(base_path_articles, "[position() <= 1]//*", sep="")
q.xpath
# strip "root"
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


