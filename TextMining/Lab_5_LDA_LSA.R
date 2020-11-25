#########################
###
### Text Mining
### Lab 5: LDA, LSA, Text Clustering
###
### Leonard Mada
### draft v.0.2


# install.packages("lsa")

library(ggplot2)

library(topicmodels) 
library(tm)
library(lsa)


setwd(...)


#################

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

#################

############
### Data ###

### r = Abstracts from Pubmed.Base.R

size = 500; # 100; 500
s.id = sample(seq_along(r), size=size, replace=FALSE)

textdata = r[s.id]

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
	stemDocument(x, language = meta(x, language))
}

control <- list(bounds = list(local = c(1, Inf)),
	language = 'english', tolower = TRUE,
	removeNumbers = TRUE, removePunctuation = TRUE, stopwords = TRUE, stripWhitespace = TRUE,
	stemWords=wordStemmer, wordLengths = c(3,20), weighting = weightTf)

tdm <- DocumentTermMatrix(corpus, control = control)

###########

###########
### LDA ###

numberOfTopics = 5
lda <- LDA(tdm, numberOfTopics) 
terms(lda)

lda@terms[1:10]

### TODO:
# - extract additional information from the topic models;


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

