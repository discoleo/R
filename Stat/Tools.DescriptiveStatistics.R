########################
###
### Leonard Mada
### [the one and only]
###
### Tools: Descriptive Statistics
###
### draft v.0.1k-fix


###############
### History ###
###############


### draft v.0.1k - v.0.1k-fix:
# - apply.html() function;
# - [improvement] xml_text vs xml_contents; [v.0.1k-impr]
# - what.types: better concept for text vs content; [v.0.1k-ref]
# - [fix] preserve tree structure; [v.0.1k-fix]
### draft v.0.1j - v.0.1j-ref4:
# - preparation for redesign:
#   moved *Hack* to bottom of file;
# - [refactor] split.stat() independent of data container;
# - [refactor] format.html.table() replaces initial hack;
# - [refactor] re-organized sections of code; [ref3 & ref4]
### draft v.0.1g - v.0.1i:
# - some fixes & better example;
# - improved example; [v.0.1g-ex]
# - convert xml directly to string; [v.0.1h]
# - as.html: preparation for redesign/refactoring;
### draft v.0.1f - v.0.1f-ref:
# - support differences in abbreviations in
#   table body and table footer;
# - [refactor] split into separate functions; [v.0.1f-ref]
### draft v.0.1e:
# - insert Reference in table body;


################
### Packages ###
################

### Review started:

### Formula based:
# 1.) atable
# 2.) gtsummary
# 3.) qwraps2

### TODO
# 4.) vtable
# 5.) DescTools

### TODO:
# stargazer, tables, tangram, rtables, dashPivottable, carpenter,
# huxtable, jtools, ...

### Output engines:
# Hmisc::latex
# officer::body_add_table
# flextable::regulartable
# rmarkdown
# ...

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
# - see section Examples;


### Basic helper functions:
view.gtsummary = function(x, len=10, view=TRUE) {
	view = if(view) interactive() else NULL;
	### Old Hack
	# as.tags.gt_tbl(as_gt(x), len=len, view = view);
	format.html.table(x, len=len, view=view);
}
as.html = function(x) {
	if(inherits(x, "shiny.tag")) {
		html = x;
	} else if(inherits(x, "gtsummary")) {
		html = gt:::as.tags.gt_tbl(as_gt(x));
	} else {
		stop("Other HTML formats: Not yet implemented!")
	}
	return(html);
}

### Split Long Results
split.stat.node = function(node, len=10, sep="<br/>", reg = "\\([-0-9 ,.]++\\)") {
	x = xml_text(node, trim=TRUE);
	npos = regexpr(reg, x, perl=TRUE);
	isMatch = (! is.na(x)) & (npos >= 0);
	LEN = attr(npos, "match.length");
	isLong  = (LEN >= len);
	isMatch = isMatch & isLong;
	if( ! isMatch) return(node);
	# Match:
	# TODO: only 1 Node vs Nodeset ?
	sMatch  = x[isMatch];
	s1 = substr(sMatch, 1, npos[isMatch] - 1);
	s2 = substr(sMatch, npos[isMatch], npos[isMatch] + LEN[isMatch]);
	# x[isMatch] = paste0(s1, sep, s2);
	# <r> = a root is needed;
	# new.line = read_xml(paste0("<r>", s1, sep, s2, "</r>"));
	# - with some effects ;-)
	new.line = read_xml(paste0("<r>", s1, sep, "<span style=\"color:#D03232\">", s2, "</span></r>"));
	# delete previous text
	xml_text(node) = "";
	for(nn in xml_contents(new.line)) {
		xml_add_child(node, nn);
	}
	return(node);
}
# [old] modifies the raw data in the gt data-structure;
split.stat.gt = function(x, len=10, sep="<br/>", reg = "\\([-0-9 ,.]++\\)") {
	# Regex: cover also negative numbers;
	nms = names(x);
	id  = grepl("^stat_", nms);
	nms = nms[id];
	for(nm in nms) {
		stat = x[[nm]];
		npos = regexpr(reg, stat, perl=TRUE);
		isMatch = (! is.na(stat)) & (npos >= 0);
		LEN = attr(npos, "match.length");
		isLong  = (LEN >= len);
		isMatch = isMatch & isLong;
		sMatch  = stat[isMatch];
		s1 = substr(sMatch, 1, npos[isMatch] - 1);
		s2 = substr(sMatch, npos[isMatch], npos[isMatch] + LEN[isMatch]);
		stat[isMatch] = paste0(s1, sep, s2);
		x[nm] = stat;
	}
	return(x);
}

### Note:
### [OLD HACK]
# - functions render_as_html() & as.tags.gt_tbl():
#   moved to the bottom of this file;


format.html.table = function(x, len=10, sep="<br/>", ..., view=FALSE) {
	html = as.html(x);
	# XML
	h2 = read_html(html$children[[2]])
	### Format table
	cells = xml_find_all(h2, "//tbody/tr/td[contains(@class,'gt_center')]");
	# update table:
	for(node in cells) {
		split.stat.node(node, len=len, sep=sep, ...);
	}
	# convert directly to text
	html$children[[2]] = read.shiny(text=as.character(h2), strip=TRUE);
	if(view) print(html, browse = interactive());
	invisible(html);
}

######################

### Add Abbreviations:

add.abbrev = function(x, abbr, label, view=TRUE, sep.eq = " = ") {
	if(length(abbr) != length(label))
		stop("Error: Abbreviations must match labels!");
	html = as.html(x);
	# XML
	h2 = read_html(html$children[[2]])
	### Footer
	nFoot = as.integer(xml_text(xml_find_all(h2, "//tfoot/tr/td/p/sup"), trim=TRUE));
	nFoot0 = max(nFoot); nFoot = nFoot0 + 1;
	# Add new Footnote:
	# - assumes only 1 table!
	xml.foot = xml_find_first(h2, "//tfoot/tr/td");
	for(id in seq(length(abbr))) {
		# Case: Abbreviation = NA
		if(is.na(abbr[id])) { abbr.txt = ""; eq.txt = ""; }
		else { abbr.txt = abbr[id]; eq.txt = sep.eq; }
		foot.html = read_xml(paste0(
			"<p class=\"gt_footnote\"><sup class=\"gt_footnote_marks\"><em>",
			nFoot, "</em></sup>", abbr.txt, eq.txt, label[id], "</p>"));
		xml.foot %>%
			xml_add_child(foot.html);
		nFoot = nFoot + 1;
	}
	### Add Ref in table:
	xml.ref = xml_find_all(h2, "//tbody/tr/td[contains(@class,'gt_left')]");
	# Names that are abbreviated:
	abbr = abbr[ ! is.na(abbr)];
	nms = if(is.null(names(abbr))) abbr else names(abbr);
	isEmpty = ! nzchar(nms); nms[isEmpty] = abbr[isEmpty];
	# update table:
	for(node in xml.ref) {
		id = which(xml_text(node, trim=TRUE) == nms);
		if(length(id) > 0) {
			foot.html = read_xml(paste0(
				"<sup class=\"gt_footnote_marks\"><em>", nFoot0 + id, "</em></sup>"));
			node %>%
				xml_add_child(foot.html);
		}
	}
	# write temp xml:
	if(FALSE) {
		out.html = tempfile("_out.tmp.html");
		# "no_declaration": seems to have no effect?
		write_xml(h2, out.html, options = c("no_declaration", "format"));
		# read new xml:
		html$children[[2]] = read.shiny(out.html, strip=TRUE);
		unlink(out.html); rm(out.html);
	} else {
		# convert directly to text
		html$children[[2]] = read.shiny(text=as.character(h2), strip=TRUE);
	}
	if(view) print(html, browse = interactive());
	invisible(html);
}

#####################

### HTML/XML
### Helper functions

read.shiny = function(file.html, text=NULL, strip=TRUE) {
	isFile = is.null(text);
	h2 = if(isFile) readLines(file.html) else text;
	# TODO: robust method
	if(strip && isFile) {
		h2 = h2[-1]; h2[length(h2)] = "</table>";
		s.tmp = h2[1]; h2[1] = substr(s.tmp, 13, nchar(s.tmp));
		h2 = paste0(h2, collapse="");
	} else if(strip) {
		h2 = sub("(?i)\\<\\!DOC[^<>]++\\>[ \n\r\t]*+", "", h2, perl=TRUE);
		h2 = sub("(?i)\\<html[^<>]*+\\>[ \n\r\t]*+", "", h2, perl=TRUE);
		h2 = sub("(?i)\\<body[^<>]*+\\>[ \n\r\t]*+", "", h2, perl=TRUE);
		# END Tag
		h2 = sub("(?i)\\</html[^<>]*+\\>[ \n\r\t]*+$", "", h2, perl=TRUE);
		h2 = sub("(?i)\\</body[^<>]*+\\>[ \n\r\t]*+$", "", h2, perl=TRUE);
	}
	# shiny html:
	attr(h2, "html") = TRUE; class(h2) = c("html", "character");
	return(h2);
}


#####################

# Note:
# - does NOT function properly if the selected nodes contain children;
apply.html = function(x, XPATH="//tbody/tr/td", FUN, ..., view=FALSE, what.type="Text") {
	html = as.html(x);
	# XML
	h2 = read_html(html$children[[2]])
	### Format table
	cells = xml_find_all(h2, XPATH);
	if(length(cells) <= 0) {
		print("No cells were updated!")
		invisible(html);
	}
	# "Nodes" = "Elements"
	what.type = pmatch(what.type, c("Text", "Nodes", "AsText"));
	if(is.na(what.type)) stop("Error: Processing type not supported!")
	txt = if(what.type == 1) xml_text(cells)
		else if(what.type == 2) lapply(cells, xml_contents) # preserve tree structure
		else lapply(cells, function(cell) as.character(xml_contents(cell)));
	txt = FUN(txt, ..., hasTags = (what.type > 1));
	### Update table:
	if(is.list(txt)) {
		# nodes which were updated
		cells = cells[txt$isUpdated];
		txt = txt$txt;
	}
	# proper update:
	# Note:
	# if(what.type == 2) assumes the nodes were updated by FUN;
	if(what.type != 2) {
	for(id in seq_along(cells)) {
		node = cells[id];
		new.line = read_xml(paste0("<r>", txt[[id]], "</r>", collapse=""));
		# delete previous text
		xml_text(node) = "";
		if(what.type > 1) {
			# remove all children;
			for(child in xml_children(node)) {
				xml_remove(child);
			}
		}
		for(nn in xml_contents(new.line)) {
			xml_add_child(node, nn);
		}
	}
	}
	# convert directly to text
	html$children[[2]] = read.shiny(text=as.character(h2), strip=TRUE);
	if(view) print(html, browse = interactive());
	invisible(html);
}

llammas.FUN = function(s, hasTags=FALSE, escape=TRUE) {
	### Note:
	# - incomplete: still misses the full visual effects
	#   & some relevant information is not yet printed:
	#   Was implemented at great expense and at the last minute using
	#   40 specially trained Ecuadorian mountain llamas,
	#   6 Venezuelan red llamas,
	#   142 Mexican whooping llamas,
	#   14 north Chilean guanacos (closely related to the llama);
	len = length(s);
	rgb = sample(seq(0, 255), 3*len, replace=TRUE);
	rgb = rgb / 255; len = len - 1;
	rgb = sapply(seq(0, len), function(id) {
		id0 = 3*id + 1;
		rgb(rgb[id0], rgb[id0 + 1], rgb[id0 + 2]);
	});
	# TODO: skip html tags;
	# Note: htmltools::htmlEscape may be required;
	s = if(escape) htmltools::htmlEscape(s) else s;
	s = paste0("<span style=\"color:", rgb, "\">", s, "</span>");
	return(s);
}


#####################
#####################

### Package gtsummary

### Split Header
# - in package: gtsummary;

header = list(
	all_stat_cols(FALSE) ~ "**{level}**<br/> N = {style_number(n)}",
	label ~ paste0("**", "Characteristic", "**")
)
header0 = list(
	stat_0 ~ "**Overall**<br/> N = {style_number(N)}"
)

################
################

################
### Examples ###
################

### Usage:
if(FALSE) {
# NOT run
some.data %>%
	# rename2(): see file Tools.Data.R;
	# - behaves almost the same to dplyr::rename();
	rename2("Full name" = "abbreviated.name") %>%
	tbl_summary(by = Dx) %>%
	modify_header(update = header) %>%
	add_p() %>%
	add_overall() %>%
	modify_header(update = header0) %>%
	# Split long results
	format.html.table(len=10) %>%
	# Add abbreviations
	add.abbrev(c("Abbr"), c("Full Name"));
}


if(FALSE) {
# NOT run / can be run
library(gtsummary)
library(xml2)
# Statistic
stats = list(all_continuous() ~ "{median} ({p25}, {p75})");
# variant:
stats = list(all_continuous() ~ "{median} ({min} - {max})");
# Data
mtcars %>%
	# rename2():
	# - see file Tools.Data.R;
	# - behaves in most cases the same as dplyr::rename();
	rename2("HP" = "hp", "Displ" = disp, "Wt (klb)" = "wt", "Rar" = drat) %>%
	# as.factor.df():
	# - see file Tools.Data.R;
	# - encode as (ordered) factor;
	as.factor.df("cyl", "Cyl ") %>%
	# the Descriptive Statistics:
	tbl_summary(by = cyl, statistic = stats) %>%
	modify_header(update = header) %>%
	add_p() %>%
	add_overall() %>%
	modify_header(update = header0) %>%
	# Hack: split long statistics
	format.html.table(len=8) %>%
	add.abbrev(
		c("Displ", "HP", "Rar", "Wt (klb)" = "Wt"),
		c("Displacement (in^3)", "Gross horsepower", "Rear axle ratio",
		"Weight (1000 lbs)"));

	# experimental: add above before format.html.table:
	# apply.html(FUN=llammas.FUN) %>%
}

### Note: XPATH
# - any table: relative path => "//"
# - only tbody => "//tbody"
# - table cells, not rows => "//tbody/tr/td"
# - exclude left column: in package gt class is "gt_left"
#   => "//tbody/tr/td[not(contains(@class, 'gt_left'))]"
#   => or "//tbody/tr/td[contains(@class, 'gt_center')]"


#########################
#########################


### [OLD]
# Code from the old Hack

# - a lot of code from the package gt is duplicated;
render_as_html.hack = function (data, ...) 
{
	# various results may contain characters that need escaping;
	data <- gt:::build_data(data = data, context = "html")
	# --- HACK ---
	# splitting the lines after the escaping;
	BLOCK = "_body";
	data[[BLOCK]] = split.stat.gt(data[[BLOCK]], ...);
	# --- END HACK ---
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
as.tags.gt_tbl.hack = function (x, len=10, ..., view = interactive()) 
{
    table_id <- gt:::dt_options_get_value(x, option = "table_id")
    if (is.na(table_id)) {
        id <- gt:::random_id()
    }
    else {
        id <- table_id
    }
    html_table = render_as_html.hack(data = x, len=len); # override original function;
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

