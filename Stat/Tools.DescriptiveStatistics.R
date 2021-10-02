########################
###
### Leonard Mada
### [the one and only]
###
### Tools: Descriptive Statistics
###
### draft v.0.1e


################
### Packages ###
################

### Review started:

# 1.) atable
# 2.) gtsummary # Formula-based
# 3.) qwraps2 # Formula-based
# 4.) vtable
# 5.) DescTools

### TODO:
# tables, tangram, rtables, dashPivottable, carpenter, ...


##################
##################

#############
### Tools ###
#############

# source("Tools.DescriptiveStatistics.R")

### Hack to split long statistics:
# - in gt: e.g. used by gtsummary;
# - functions which require overriding:
#   as.tags.gt_tbl & render_as_html;
# - the problem is in function: gt:::process_text;

### Usage:
# as.tags.gt_tbl(as_gt(tbl));


view.gtsummary = function(x, len=10, view=TRUE) {
	view = if(view) interactive() else NULL;
	as.tags.gt_tbl(as_gt(tbl), len=len, view = view);
}
split.stat = function(x, len=10, sep="<br/>", reg = "\\([0-9 ,]++\\)", BLOCK="_body") {
	# BLOCK = "_body"; # "table_body"
	nms = names(x[[BLOCK]]);
	id  = grepl("^stat_", nms);
	nms = nms[id];
	for(nm in nms) {
		stat = x[[BLOCK]][[nm]];
		npos = regexpr(reg, stat, perl=TRUE);
		isMatch = (! is.na(stat)) & (npos >= 0);
		LEN = attr(npos, "match.length");
		isLong  = (LEN >= len);
		isMatch = isMatch & isLong;
		sMatch  = stat[isMatch];
		s1 = substr(sMatch, 1, npos[isMatch] - 1);
		s2 = substr(sMatch, npos[isMatch], npos[isMatch] + LEN[isMatch]);
		stat[isMatch] = paste0(s1, sep, s2);
		x[[BLOCK]][nm] = stat;
	}
	return(x);
}
render_as_html = function (data, ...) 
{
	# various results may contain characters that need escaping;
	data <- gt:::build_data(data = data, context = "html")
	# splitting the lines after the escaping;
	data = split.stat(data, ...);
	data <- gt:::add_css_styles(data = data)
    caption_component <- gt:::create_caption_component_h(data = data)
    heading_component <- gt:::create_heading_component_h(data = data)
    columns_component <- gt:::create_columns_component_h(data = data)
    body_component <- gt:::create_body_component_h(data = data)
    source_notes_component <- gt:::create_source_notes_component_h(data = data)
    footnotes_component <- gt:::create_footnotes_component_h(data = data)
    table_defs <- gt:::get_table_defs(data = data);
    html_tbl = as.character(htmltools::tags$table(class = "gt_table", 
        style = table_defs$table_style, caption_component, table_defs$table_colgroups, 
        heading_component, columns_component, body_component, 
        source_notes_component, footnotes_component));
	invisible(html_tbl);
}
as.tags.gt_tbl = function (x, len=10, ..., view = interactive()) 
{
    table_id <- gt:::dt_options_get_value(x, option = "table_id")
    if (is.na(table_id)) {
        id <- gt:::random_id()
    }
    else {
        id <- table_id
    }
    html_table = render_as_html(data = x, len=len); # override original function;
    css <- gt:::compile_scss(data = x, id = id)
    container_overflow_x <- gt:::dt_options_get_value(x, option = "container_overflow_x")
    container_overflow_y <- gt:::dt_options_get_value(x, option = "container_overflow_y")
    container_width <- gt:::dt_options_get_value(x, option = "container_width")
    container_height <- gt:::dt_options_get_value(x, option = "container_height")
    html_tbl <- htmltools::tags$div(id = id, htmltools::tags$style(htmltools::HTML(css)), 
        style = htmltools::css(`overflow-x` = container_overflow_x, 
            `overflow-y` = container_overflow_y, width = container_width, 
            height = container_height), htmltools::HTML(html_table));
	if( ! is.null(view))
		print(html_tbl, browse = view, ...);
    invisible(html_tbl);
}

### Add Abbreviations:

add.abbrev = function(x, abbr, label, view=TRUE) {
	if(length(abbr) != length(label))
		stop("Error: Abbreviations must match labels!")
	if(inherits(x, "shiny.tag")) {
		html = x;
	} else if(inherits(x, "gtsummary")) {
		html = view.gtsummary(tbl, view=FALSE);
	}
	# XML
	h2 = read_html(html$children[[2]])
	# Footer
	nFoot = as.integer(xml_text(xml_find_all(h2, "//tfoot/tr/td/p/sup"), trim=TRUE));
	nFoot0 = max(nFoot); nFoot = nFoot0 + 1;
	# Add new Footnote:
	xml.foot = xml_find_first(h2, "//tfoot/tr/td");
	for(id in seq(length(abbr))) {
		foot.html = read_xml(paste0(
			"<p class=\"gt_footnote\"><sup class=\"gt_footnote_marks\"><em>",
			nFoot, "</em></sup>", abbr[id], " = ", label[id], "</p>"));
		xml.foot %>%
			xml_add_child(foot.html);
		nFoot = nFoot + 1;
	}
	# Add Ref in table:
	xml.ref = xml_find_all(h2, "//tbody/tr/td[contains(@class,'gt_left')]");
	for(node in xml.ref) {
		id = which(xml_text(node, trim=TRUE) == abbr);
		if(length(id) > 0) {
			foot.html = read_xml(paste0(
				"<sup class=\"gt_footnote_marks\"><em>", nFoot0 + id, "</em></sup>"));
			node %>%
				xml_add_child(foot.html);
		}
	}
	# write temp xml:
	out.html = tempfile("_out.tmp.html");
	write_xml(h2, out.html, options = c("no_declaration", "format"));
	# read new xml:
	h2 = readLines(out.html);
	# TODO: robust method
	h2 = h2[-1]; h2[length(h2)] = "</table>";
	s.tmp = h2[1]; h2[1] = substr(s.tmp, 13, nchar(s.tmp));
	# shiny html:
	h2 = paste0(h2, collapse="");
	attr(h2, "html") = TRUE; class(h2) = c("html", "character");
	html$children[[2]] = h2;
	unlink(out.html); rm(out.html);
	if(view) print(html, browse = interactive());
	invisible(html);
}


### Split Header
# - in package: gtsummary;

header = list(
	all_stat_cols(FALSE) ~ "**{level}**<br/> N = {style_number(n)}",
	label ~ paste0("**", "Characteristic", "**")
)
header0 = list(
	stat_0 ~ "**Overall**<br/> N = {style_number(N)}"
)


### Usage:
if(FALSE) {
# NOT run
some.data %>%
	rename("Full name" = "abbreviated.name") %>%
	tbl_summary(by = Dx) %>%
	modify_header(update = header) %>%
	add_p() %>%
	add_overall() %>%
	modify_header(update = header0);
}
