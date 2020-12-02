#########################
###
### Text Mining
### Lab 5: LDA, LSA, Text Clustering
###
### Leonard Mada
### draft v.0.3a


### History

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

library(ggplot2)

# processing
library(dplyr)
library(tidyr)
library(tidytext)
# text mining
library(topicmodels) 
library(tm)
library(lsa)


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
clean.txt = function(x, rmBrackets=TRUE, rmPunctuation=TRUE) {
	x = gsub("&gt;", " > ", x)
	x = gsub("&lt;", " < ", x)
	x = gsub("[\r\n\t \uA0]+", " ", x)
	x = gsub("P[.]?(?= falc| vivax)", "Plasmodium", x, perl=TRUE)
	x = gsub("'s(?=[ .)])", "", x, perl=TRUE)
	x = gsub("'(?=[ .)])", "", x, perl=TRUE)
	x = gsub("(?<=[ (])'", "", x, perl=TRUE)
	x = gsub("(?<=[0-9])-(?=[0-9])", " ", x, perl=TRUE)
	x = gsub("\"", "", x)
	x = gsub("\\$(?=[0-9])", " $NN $", x, perl=TRUE)
	if(rmBrackets) x = gsub("[()]", "", x)
	if(rmPunctuation) x = gsub("[,.;](?= |$)", "", x, perl=TRUE)
	return(x)
}
stem.txt = function(x) {
	x = gsub("(?<=[rnltg]|[lt]e)s(?= |$)|(?<=to)es(?= |$)", "", x, perl=TRUE)
	x = gsub("(?<=[f])ied(?= |$)", "y", x, perl=TRUE)
	x = gsub("(?<=[tgdpfhlmns])ies(?= |$)", "y", x, perl=TRUE)
	x = gsub("(?<!fa|spe)cies(?= |$)", "cy", x, perl=TRUE)
	x = gsub("(?<!se|ca)ries(?= |$)", "ry", x, perl=TRUE)
	x = gsub("diseases", "disease", x)
	x = gsub("vaccines", "vaccine", x)
	x = gsub("areas(?= |$)", "area", x, perl=TRUE)
	x = gsub("assays(?= |$)", "assay", x, perl=TRUE)
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
	"showed", "found", "suggest")

find.nextWord = function(x, s, sort=TRUE) {
	m = regexec(paste0(s, " ([^ ]+)", collapse=""), x)
	m.txt = regmatches(x, m)
	isM = sapply(regmatches(x, m), function(l) length(l) > 0)
	m.txt = sapply(m.txt[isM], function(txt) txt[2])
	if(sort) m.txt = sort(m.txt)
	return(m.txt)
}

#################

############
### Data ###

### r = Abstracts from Pubmed.Base.R

size = 100; # 100; 500
s.id = sample(seq_along(r), size=size, replace=FALSE)

textdata = r[s.id]
textdata = clean.text(textdata)

# Abbreviations
table(grepl(" A[.]", r))
table(grepl(" [B-OQ-Z][.]", r))

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
	removeNumbers = TRUE, removePunctuation = TRUE, stopwords = TRUE, stripWhitespace = TRUE,
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

##################
### Clustering ###
##################

x = read.csv("Malaria.Abstracts.csv", stringsAsFactor=F)
x = x[,1]

x = clean.txt(x)
x = stem.txt(x)
x = add.text(x)
print.art(x[1], max=20)

# TODO: lowercase;

### Exclude Corrections
isCorrect = grepl("^\\[This correct", x)
sum(isCorrect)
x = x[ ! isCorrect]

### Corpus
corpus <- VCorpus(VectorSource(x),
	readerControl = list(language = 'english'))

# StopWords
match(stopWords, stopwords("en"))

dtm <- DocumentTermMatrix(corpus, control = list(
	weighting = function(x) weightTfIdf(x, normalize = FALSE),
	stopwords = c(stopwords("en"), stopWords)))
dtm
dtm$dimnames$Terms[1:100]
# dtm$dimnames$Terms[grepl("ies$", dtm$dimnames$Terms)]

# Remove sparse terms
dtm <- removeSparseTerms(dtm, sparse=0.92)
dtm
dtm$dimnames$Terms[1:100]

dtm_tx <- weightTfIdf(dtm)
mat_dtm <- as.matrix(dtm_tx)
rownames(mat_dtm) <- 1:nrow(mat_dtm)

# Term Normalization
normalise_dtm <- function(y) y/apply(y, MARGIN=1,
	FUN=function(k) sum(k^2)^.5)
dtm_norm <- normalise_dtm(mat_dtm)


### Clustering
cl <- kmeans(dtm_norm, 8)
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


