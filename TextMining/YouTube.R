
# Global Immunotalks
# https://www.youtube.com/@GLOBALIMMUNOTALKS/videos

# setwd(...)

x = read_html("YouTube.html");


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

