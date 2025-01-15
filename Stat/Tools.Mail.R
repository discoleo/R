
### Tools: eMail

### Ref:
# https://en.wikipedia.org/wiki/MIME
# https://en.wikipedia.org/wiki/Base64

decode.ro = function(x, rm.breaks = FALSE) {
	x = gsub("=[\n\r]+", "", x);
	### Chars:
	# Note:
	# - NOT robust: quick hack!
	# - 3 or 4 byte chars must be processed first (greedy);
	x = gsub("=C4=83", "a", x); # a\/
	x = gsub("=C3=AE", "i", x);
	x = gsub("=C3=8E", "I", x);
	x = gsub("=C3=A2", "a", x); # a^
	x = gsub("=C8=9B", "t", x);
	x = gsub("=C8=9A", "T", x);
	x = gsub("=C8=99", "s", x);
	x = gsub("=C8=98", "S", x);
	# x = gsub("=C8=\n=9B", "t", x);
	# x = gsub("=C4=\n=83", "a", x); # a
	
	### Other:
	x = gsub("=E2=80=93", "-", x);
	x = gsub("=3D", "=", x);
	return(x);
}
cat.QPhtml = function(x) {
	x = decode.ro(x);
	x = strsplit(x, "(?<=[>])(?=[<])", perl=TRUE);
	tmp = lapply(x, cat, sep = "\n");
	invisible(x);
}

### QP: Quoted Printable
decode.QPcode = function(x) {
	# x = utf8ToInt(x);
	x = strsplit(x, "=")[[1]];
	x = x[nchar(x) > 0];
	if(length(x) == 0) return("");
	x = as.raw(as.hexmode(x));
	return(x);
}
as.QPchar = function(x) {
	rawToChar(decode.QPcode(x));
}
