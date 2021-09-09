
format.style = function(s, style=NULL) {
	# TODO: style;
	# current default: color;
	reg = "(?<![\\\\])\\%([^%]*+)\\%";
	gr  = gregexpr(reg, s, perl=TRUE);
	grc = lapply(gr, function(gr) {
		if(length(gr) == 1) gr else {
			gs = as.array(attr(gr, "capture.start"));
			gl = as.array(attr(gr, "capture.length"));
			attr(gs, "match.length") = gl;
			attr(gs, "useBytes") = TRUE;
			return(gs);
		}
	})
	sgr = regmatches(s, grc);
	sgr = lapply(sgr, function(s) if(length(s) == 0) s else {
		isEndTag = (nchar(s) == 0);
		s[isEndTag] = "</span>";
		s[ ! isEndTag] = paste0("<span style=\"color:", s[ ! isEndTag], ";\">");
		return(s);
	})
	regmatches(s, gr) = sgr;
	hasStyle = sapply(sgr, function(s) (length(s) > 0));
	attr(s, "html") = hasStyle;
	return(s);
}

### HTML Tags

p = function(s, style=NULL) {
	# uses library("htmltools")
	s = format.style(s);
	hasStyle = attr(s, "html");
	lapply(seq(length(s)), function(id) {
		s = s[[id]]; attr(s, "html") = hasStyle[[id]];
		s = htmltools::p(s);
		if( ! is.null(style)) {
			s$attribs$style = style;
		}
		return(s);
	})
}


#################

### Test


p(c("A %red%red tag%%", "blue"))

