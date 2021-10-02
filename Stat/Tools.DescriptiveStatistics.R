########################
###
### Leonard Mada
### [the one and only]
###
### Tools: Descriptive Statistics
###
### draft v.0.1c


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


### Split Header
# - in gtsummary;

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
