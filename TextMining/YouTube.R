
# Global Immunotalks
# https://www.youtube.com/@GLOBALIMMUNOTALKS/videos

# setwd(...)

library(xml2)

x = read_html("YouTube.html");
# x = read_html("GLOBAL IMMUNOTALKS - YouTube.html")
# x = read_html("Brandeis Math Bio - YouTube.html")
# Note: html content needs to be manually copied from the web page to the file;

### All Nodes
nodes = xml_find_all(x, ".//a[@id=\"video-title-link\"]")

### Title
title = xml_text(nodes);
# cat(xml_text(nodes), sep = "\n")


### Link
link = xml_text(xml_find_all(nodes, ".//@href"));
link = paste0("https://www.youtube.com", link)
# cat(link, sep="\n")

x.df = data.frame(ID = seq(length(title)), Title = title, href = link);
head(x.df)

write.csv(x.df, file = "ImmunoTalks.csv", row.names = FALSE)


#####################

### Brandeis Playlist

### All Nodes
nodes = xml_find_all(x, ".//h3[contains(@class, \"ytd-playlist-video-renderer\")]/a")

### Title
title = gsub("^[ \t\r\n]+|[ \t\r\n]+$", "", xml_text(nodes));


### Link
link = xml_text(xml_find_all(nodes, ".//@href"));
link = paste0("https://www.youtube.com", link)
# cat(link, sep="\n")

x.df = data.frame(ID = seq(length(title)), Title = title, href = link);
head(x.df)

write.csv(x.df, file = "ImmunoTalks.csv", row.names = FALSE)

