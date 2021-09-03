

##################

### Package igraph

library(igraph)


### Connected components

bind.components = function(net, mode="strong") {
    cnet = components(net, mode = mode);
    lapply(seq(cnet$no),function(id)  attr(cnet$membership, "names")[cnet$membership == id]);
}

print.components = function(net, mode="strong", print=TRUE) {
    comp = bind.components(net, mode=mode);
    str = sapply(comp, function(s) paste0(s, collapse=", "));
    str = paste0("Component ", seq(length(str)), ": ", str, "\n");
    if(print) {cat(str, sep=""); cat("\n");}
    invisible(str);
}

### Examples:

# Cycle: 1 -> 2 -> 4 -> 1
net = graph_from_data_frame(cbind(c(1,1,2,3,4,6,7), c(2,3,4,5,1,7,6)), directed=TRUE)
### Weak components
bind.components(net, mode = "weak");
# Pretty-print:
print.components(net, mode = "weak");
### Strong components
bind.components(net);
# Pretty-print:
print.components(net);


