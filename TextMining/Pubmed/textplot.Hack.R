
textplot_dependencyparser.default = function (x, title = "Dependency Parser",
	subtitle = "tokenisation, parts of speech tagging & dependency relations", 
    vertex_color = "darkgreen", edge_color = "red", size = 3, 
    base_family = "", layout = "linear", nudge_y = -0.35, ...) 
{
    stopifnot(is.data.frame(x) & all(c("sentence_id", "token_id", 
        "head_token_id", "dep_rel", "token_id", "token", "upos") %in% 
        colnames(x)))
    dep_rel <- token <- upos <- NULL
    requireNamespace("ggraph")
    requireNamespace("ggplot2")
    requireNamespace("igraph")
    x <- x[!is.na(x$head_token_id), ]
    x <- x[x$sentence_id %in% min(x$sentence_id), ]
    edges <- x[x$head_token_id != 0, c("token_id", "head_token_id", 
        "dep_rel")]
    edges$label <- edges$dep_rel
    g <- igraph::graph_from_data_frame(edges, vertices = x[, 
        c("token_id", "token", "upos")], directed = TRUE)
    ggraph::ggraph(g, layout = layout) +
		ggraph::geom_edge_arc(ggplot2::aes(label = dep_rel, 
        vjust = -0.2), arrow = grid::arrow(length = ggplot2::unit(4, 
        "mm"), ends = "last", type = "closed"), end_cap = ggraph::label_rect("wordswordswords"), 
        label_colour = edge_color, check_overlap = TRUE, label_size = size) + 
        ggraph::geom_node_label(ggplot2::aes(label = token), 
            col = vertex_color, size = size, fontface = "bold") + 
        ggraph::geom_node_text(ggplot2::aes(label = upos), nudge_y = nudge_y, 
            size = size) + ggraph::theme_graph(base_family = base_family) + 
        ggplot2::labs(title = title, subtitle = subtitle)
}

