#########################
###
### Text Mining
### Project & Lab 5:
### LDA, LSA, Text Clustering
###
### Leonard Mada
### draft v.0.3d


### History

### draft v.0.3d:
# - improved stemming of nouns;
### draft v.0.3c:
# - Clustering: top terms;
### draft v.0.3b:
# - Clustering: Hierarchical;
### draft v.0.3a:
# - Clustering: k-Means;
# - more text-processing helper functions;
### draft v.0.2 - v.0.2b:
# - first draft on Github;
# - added minimalistic visualization of the LDA topic model; [v.0.2b]


# original on:
# https://github.com/discoleo/R/blob/master/TextMining/Lab_5.LDA_LSA.R


#########################
#########################

# install.packages("lsa")
# install.packages("proxy")

library(ggplot2)

# processing
library(dplyr)
library(tidyr)
library(tidytext)
# text mining
library(topicmodels) 
library(tm)
library(lsa)
# Dist Function: Clustering
library(proxy)


setwd(...)


####################

### helper functions

print.art = function(r, tlen=81, max=0) {
	len = nchar(r)
	tlen = 81;
	part = len %/% tlen
	id = (0:part)*tlen;
	if(max > 0) {
		id = id[1:max]
	}
	s = sapply(id, function(id) substr(r, id + 1, id + tlen))
	sapply(s, function(s) {cat(s); cat("\n");})
	cat("\n")
	invisible(s)
}
jitter.txt = function(x, len=4, NOP=FALSE, rnd=FALSE, seed=1234) {
	if(length(x) == 1) {
		x = id = seq(x);
	} else {
		id = seq_along(x);
	}
	if(NOP) return(x);
	if(rnd) {
		set.seed(seed)
		len = sample(1:len, length(id) %/% 2, replace=TRUE)
		jt.txt = sapply(len, function(len) paste(rep(" ", len), sep="", collapse=""))
	} else {
		jt.txt = paste(rep(" ", len), sep="", collapse="")
	}
	x[id %% 2 == 0] = paste0(x[id %% 2 == 0], jt.txt)
	return(x)
}
corpus.f = function(corpus, ...) {
	# TODO: implement as options;
	x.corpus = corpus;
	x.corpus <- tm_map(x.corpus, tolower)
	x.corpus <- tm_map(x.corpus, removePunctuation)
	x.corpus <- tm_map(x.corpus, removeNumbers) 
	x.corpus <- tm_map(x.corpus, function(x) removeWords(x, stopwords("english"))) 
	x.corpus <- tm_map(x.corpus, stemDocument, language = "english") 
	x.corpus <- tm_map(x.corpus, PlainTextDocument)
	return(x.corpus)
}
dtm.f = function(x) {
	corpus <- VCorpus(VectorSource(x),
	readerControl = list(language = 'english'))
	# DTM
	dtm <- DocumentTermMatrix(corpus, control = list(
		weighting = function(x) weightTfIdf(x, normalize = FALSE),
		stopwords = c(stopwords("en"), stopWords)))
	dtm
}
### Post-Processing
sparse.dtm = function(dtm, sparse=0.92, print=TRUE, print.n=100) {
	# Remove sparse terms
	dtm <- removeSparseTerms(dtm, sparse=sparse)
	if(print) {
		print(dtm)
		print(head(dtm$dimnames$Terms, print.n))
	}
	dtm;
}
normalize.dtm = function(dtm) {
	# Normalize
	dtm.wgt <- weightTfIdf(dtm)
	dtm.m <- as.matrix(dtm.wgt)
	rownames(dtm.m) <- 1:nrow(dtm.m)
	# Term Normalization
	normalise_dtm <- function(y) {
		y / apply(y, MARGIN=1, FUN=function(k) sum(k^2)^.5)
	}
	dtm_norm <- normalise_dtm(dtm.m)
}
### Text
clean.txt = function(x, rmBrackets=TRUE, rmPunctuation=TRUE) {
	x = gsub("&amp;", "&", x)
	x = gsub("&gt;", " > ", x)
	x = gsub("&lt;", " < ", x)
	x = gsub("&#039;", "'", x)
	x = gsub("[\r\n\t \uA0]+", " ", x)
	x = gsub("P[.]?(?= falc| vivax)", "Plasmodium", x, perl=TRUE)
	x = gsub("'s(?=[ .)])", "", x, perl=TRUE)
	x = gsub("'(?=[ .)])", "", x, perl=TRUE)
	x = gsub("(?<=[ (])'", "", x, perl=TRUE)
	x = gsub("(?<=[0-9])-(?=[0-9])", " ", x, perl=TRUE)
	x = gsub("\"", "", x)
	x = gsub("\\: ", " : ", x)
	x = gsub("\\$(?=[0-9])", " $NN $", x, perl=TRUE)
	if(rmBrackets) x = gsub("[()]", "", x)
	if(rmPunctuation) x = gsub("[,.;](?= |$)", "", x, perl=TRUE)
	return(x)
}
stem.txt = function(x) {
	x = gsub("(?<=[rnltg]|[lt]e)s(?= |$)|(?<=to)es(?= |$)", "", x, perl=TRUE)
	x = gsub("(?<=[a-z]{3}ce(?<!ice))s(?= |$)", "", x, perl=TRUE)
	x = gsub("(?<=[a-z]{3}de(?<!ede|oide))s(?= |$)", "", x, perl=TRUE)
	x = gsub("(?<=[a-z]ge)s(?= |$)", "", x, perl=TRUE) # fails: Georges
	x = gsub("(?<=[f])ied(?= |$)", "y", x, perl=TRUE)
	x = gsub("(?<=[tgdpfhlmns])ies(?= |$)", "y", x, perl=TRUE)
	x = gsub("(?<!fa|[sS]pe)cies(?= |$)", "cy", x, perl=TRUE)
	x = gsub("(?<!se|ca)ries(?= |$)", "ry", x, perl=TRUE)
	x = gsub("(?<=rm)s(?= |$)", "", x, perl=TRUE) # {farms, terms}; could fail: HRMS
	x = gsub("diseases", "disease", x)
	x = gsub("(?<=[Vv])accines", "accine", x, perl=TRUE)
	x = gsub("scores", "score", x)
	x = gsub("biotics", "biotic", x)
	x = gsub("(?<=[a-z]{3}d(?<!aid))s(?= |$)", "", x, perl=TRUE) # exclude: (?i)AIDS
	x = gsub("olds(?= |$)", "old", x, perl=TRUE)
	x = gsub("areas(?= |$)", "area", x, perl=TRUE)
	x = gsub("zones(?= |$)", "zone", x, perl=TRUE)
	x = gsub("assays(?= |$)", "assay", x, perl=TRUE)
	x = gsub("(?<=nd)s(?= |$)", "", x, perl=TRUE) # fails: NLs, keep "seconds";
	x = gsub("RDTs(?= |$)", "RDT", x, perl=TRUE)
	x = gsub("genes(?= |$)", "gene genes", x, perl=TRUE) # keep plural;
	x = gsub("enzymes(?= |$)", "enzyme enzymes", x, perl=TRUE) # keep plural;
	return(x)
}
add.txt = function(x) {
	x = gsub("IL(?=-[0-9])", "IntL IL", x, perl=TRUE)
	x = gsub("sero-?(?=preval)", "prevalence sero", x, perl=TRUE)
}
stopWords = c("also", "although", "can", "may", "among", "respectively", "95%",
	"due", "one", "two", "three", # TODO: DM
	"showed", "found", "suggest", "whilst", "whereas", "yet",
	"however", "thereby", "thereof", "therefore", "though")

find.nextWord = function(x, s, next.reg=" ([^ ]+)", sort=TRUE, perl=TRUE, case.ignore=TRUE) {
	case.txt = if(case.ignore) "(?i)" else "";
	m = regexec(paste0(case.txt, s, next.reg, collapse=""), x, perl=perl)
	m.txt = regmatches(x, m)
	isM = sapply(m.txt, function(l) length(l) > 0)
	m.txt = sapply(m.txt[isM], function(txt) txt[2])
	if(sort) m.txt = sort(m.txt)
	return(m.txt)
}
find.nouns = function(x, noun.reg, base="[a-z]", sort=TRUE, case.ignore=FALSE) {
	# plural of potential nouns (+ verbs ending in 's');
	if(missing(noun.reg)) noun.reg = paste0(".(", base, ")s$")
	case.txt = if(case.ignore) "(?i)" else "";
	m = regexec(paste0(case.txt, noun.reg, collapse=""), x, perl=TRUE)
	m.txt = regmatches(x, m)
	isM = sapply(m.txt, function(l) length(l) > 0)
	m.txt = sapply(m.txt[isM], function(txt) txt[2])
	txt = cbind(m.txt, x[isM])
	if(sort) txt = txt[order(txt[,1], txt[,2]) , ]
	return(txt)
}
### Clusters
top.terms = function(cl, top=8) {
	ord = apply(cl$centers, 1, function(w) tail(order(w), top))
	ord = apply(ord, 2, rev) # top: most important term;
	m = attr(cl$centers, "dimnames")[[2]][ord]
	matrix(m, nrow=top)
}

#################

#################
### Section 1 ###
#################


############
### Data ###

### r = Abstracts from Pubmed.Base.R

size = 100; # 100; 500
s.id = sample(seq_along(r), size=size, replace=FALSE)

textdata = r[s.id]
textdata = clean.txt(textdata)

# Abbreviations
table(grepl(" A[.]", r))
table(grepl(" [B-OQ-Z][.]", r))
table(find.nextWord(x, "", "([a-z]{1,20}ges) "))

print.art(textdata[grepl("[AB-OQ-Z][.]", textdata)][1])


#######################

#######################
### Text Processing ###

### Corpus
corpus <- VCorpus(VectorSource(as.character(textdata)),
	readerControl = list(language = 'english'))

# the simple Corpus does NOT work with the LSA code!
# corpus <- Corpus(VectorSource(as.character(textdata)), readerControl = list(language = 'english'))


###########
### TDM ###

wordStemmer = function(x, language = 'english') {
	# stemming = _this_function();
	print(x[1:2]); # TODO: NO processing!
	stemDocument(x, language = meta(x, language))
}

control <- list(bounds = list(local = c(1, Inf)),
	language = 'english', tolower = TRUE,
	removeNumbers = TRUE, removePunctuation = TRUE, stripWhitespace = TRUE,
	stopwords = TRUE,
	stemming = TRUE, wordLengths = c(3,20), weighting = weightTf)

# DocumentTermMatrix != TermDocumentMatrix
tdm <- DocumentTermMatrix(corpus, control = control)

tdm$dimnames$Terms[1:100]

###########

###########
### LDA ###

numberOfTopics = 5
lda <- LDA(tdm, numberOfTopics) 
terms(lda)

lda@terms[1:10]

### TODO:
# - extract additional information from the topic models;

### Analysis
w.topics <- tidy(lda)
w.topics

top_terms <- w.topics %>%
  group_by(topic) %>%
  top_n(10, beta) %>%
  ungroup() %>%
  arrange(topic, -beta)

top_terms

print(top_terms, n=60)

### Visualization

theme_set(theme_bw())

top_terms %>%
  mutate(term = reorder_within(term, beta, topic)) %>%
  ggplot(aes(term, beta)) +
  geom_bar(stat = "identity") +
  scale_x_reordered() +
  facet_wrap(~ topic, scales = "free_x")


#################
### LDA + VEM ###

# https://www.rdocumentation.org/packages/topicmodels/versions/0.2-11/topics/TopicModelcontrol-class

k = numberOfTopics;
control_LDA_VEM = list( estimate.alpha = TRUE, alpha = 50/k, estimate.beta = TRUE,
	verbose = 0, prefix = tempfile(), save = 0, keep = 0, seed = as.integer(Sys.time()),
	nstart = 1, best = TRUE, var = list(iter.max = 500, tol = 10^-6),
	em = list(iter.max = 1000, tol = 10^-4), initialize = "random")

lda = LDA(tdm, k, method = "VEM", control = control_LDA_VEM)
terms(lda)

######################
######################

###########
### LSA ###

### TDM

x.corpus = corpus.f(corpus)
doc = inspect(x.corpus[1:10])

### TDM
# tdm.m = as.matrix(tdm)
tdm = TermDocumentMatrix(x.corpus)
tdm.m <- as.matrix(tdm)
tdm.lsa <- lw_bintf(tdm.m) * gw_idf(tdm.m)


view <- as.factor(paste0("A", sort(s.id)))
head(view)


### LSA
lsaSpace <- lsa(tdm.lsa)


dist_lsa <- dist(t(as.textmatrix(lsaSpace)))
# dist_lsa


### Plot
fit <- cmdscale(dist_lsa, eig = TRUE, k = 2) # !!!
points.df <- data.frame(x = fit$points[, 1], y = fit$points[, 2])

# png(file="Malaria.Articles.Dist.png")
ggplot(points.df, aes(x = x, y = y)) +
	geom_point(data = points.df, aes(x = x, y = y, size=5, color=view), show.legend = FALSE) +
	geom_text(data = points.df, aes(x = x, y = y - 1.2, label = view))

# dev.off()

print.art(r[570], max=10)
print.art(r[906], max=10)

print.art(r[949], max=10)

# with initial non-document TDM:
# rownames(fit$points)[40]


########################
########################

#################
### Section 2 ###
#################


##################
### Clustering ###
##################

### Data

x = read.csv("Malaria.Abstracts.csv", stringsAsFactor=F)
x = x[,1]

x = clean.txt(x)
x = stem.txt(x)
x = add.txt(x)
print.art(x[1], max=20)

# TODO: lowercase;

### Exclude Corrections
isCorrect = grepl("^\\[This correct", x)
sum(isCorrect)
x = x[ ! isCorrect]


###################

### Section 2.A.)

### Corpus
corpus <- VCorpus(VectorSource(x),
	readerControl = list(language = 'english'))

### DTM: StopWords
dtm <- DocumentTermMatrix(corpus, control = list(
	weighting = function(x) weightTfIdf(x, normalize = FALSE),
	stopwords = c(stopwords("en"), stopWords)))
dtm
dtm$dimnames$Terms[1:100]
# dtm$dimnames$Terms[grepl("ies$", dtm$dimnames$Terms)]


### Post-Processing

# Remove sparse terms
dtm <- sparse.dtm(dtm, 0.92)

# Normalize
dtm.n = normalize.dtm(dtm)

##################
### Clustering ###

### K-Means
cl <- kmeans(dtm.n, 8)
str(cl)

# Check the number of objects in each cluster
table(cl$cluster)

# Check the cluster assigned to each object
head(cl$cluster)
# check the center of each clusters
cl$centers[, 1:20]
# Within cluster sum of squares by cluster
head(cl$withinss)


table(find.nextWord(x, "based"))
table(find.nextWord(x, "gene"))
table(find.nextWord(x, "(?= [^ ]+ 86Y)"))


table(grepl("(?i)Anophele|Cul(?:ex|icidae)", x))

table(grepl("(?i)urban|city|rural|region", x))
table(grepl("(?i)river|estuary", x))

table(grepl("(?i)reservoir", x))

table(grepl("(?i)18S", x))
table("R"=grepl("(?i)ribosom", x), "S"=grepl("(?i)18S", x))
# print.art(x[grepl("(?i)S18", x)][1])
table(grepl("(?i)Hsp", x))
table("Ch"=grepl("(?i)chaperon", x), "H"=grepl("(?i)Hsp", x))
table(grepl("(?i)target", x))


#################

### Text Analysis

# StopWords
match(stopWords, stopwords("en"))

head(find.nouns(dtm$dimnames$Terms), 100)

head(find.nouns(dtm$dimnames$Terms, base="[a-z]d(?<!aid)"), 100)

#################
#################

### Section 2.B.)

dtm = dtm.f(x[grepl("(?i)target", x)])
dtm

# len = length(dtm$dimnames$Terms)
tail(dtm$dimnames$Terms, 200)

### Sparse
dtm = sparse.dtm(dtm, 0.93)

### Normalize
dtm.n = normalize.dtm(dtm)

##################
### Clustering ###

### K-Means
cl <- kmeans(dtm.n, 8)
str(cl)

# Check the number of objects in each cluster
table(cl$cluster)

# check the center of each clusters
cl$centers[, 1:20]

### Topics
# top terms defining the topics/clusters
top.terms(cl)


################
### Hierarchical

# library(proxy)

### Complexity (O(n^2))
distance <- dist(dtm_norm, method="cosine")
hc <- hclust(distance, method="average")

plot(hc, labels=jitter.txt(length(x), rnd=T, len=10, seed=1))

# draw dendogram with red borders around the 5 clusters 
rect.hclust(hc, k=5, border="red") # 3 major clusters;

group_clust <- cutree(hc, k=5) # cut tree into 5 clusters

