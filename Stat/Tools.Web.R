

### Web Scraping

# Tables from a dynamic web page;

library(rvest)
# Requires also package chromote;


### Functions

# n = Number of Rows per Table-Page;
# idCheck = Column used to check Load-failure (should contain unique values);
# wait    = Time (in s) to wait to render web page;
# Note:
# - If nrow(TABLE) == exact multiple of n,
#   then function will fail; max.pages must be properly set (to pages - 1);
scrap.byClick = function(url, idTable, nextPage, n = 100, idCheck = 1,
		wait = 1, max.pages = 500) {
	#
	nRowCheck = 1:3;
	#
	if(! inherits(url, "LiveHTML")) {
		x = read_html_live(url);
		Sys.sleep(wait);
	} else x = url;
	# 
	tbl = html_node(x, idTable);
	if(is.na(tbl)) {
		warning("No table! Consider increasing wait-time.");
		return(NULL);
	}
	tbl = html_table(tbl);
	#
	if(nrow(tbl) == 0) {
		# ERROR?
		warning("NO Table!");
		return(NULL);
	} else if(nrow(tbl) < n) {
		# END of table
		return(tbl);
	}
	# Check:
	hasCheck = ! is.null(idCheck);
	prev     = ifelse(hasCheck, tbl[nRowCheck, idCheck], 0);
	### Next Page in Table:
	retryScrap = function(x, hasCheck, prev) {
		Sys.sleep(wait);
		for(iRetr in seq(3)) {
			tmp2 = html_node(x, idTable);
			if(! is.na(tmp2)) break;
		}
		if(is.na(tmp2)) {
			# ERROR:
			warning("Failed to load table!");
			return(list(ERR = 1));
		}
		# Table on Next Page:
		tmp2 = html_table(tmp2);
		if(nrow(tmp2) == 0) {
			# ERROR?
			warning("NO Table!");
			return(list(ERR = 1));
		} else if(nrow(tmp2) < n) {
			# END of table
			return(list(ERR = 2, Tbl = tmp2));
		}
		if(hasCheck) {
			if(all(prev == tmp2[nRowCheck, idCheck])) {
				# Failed to load next page!
				return(list(ERR = 3));
			}
		}
		return(list(ERR = 0, Tbl = tmp2));
	}
	cat("Page: ");
	for(i in seq(max.pages)) {
		cat(i, ",", sep = ""); flush.console();
		x$click(nextPage);
		for(iRetry in seq(3)) {
			tbl2 = retryScrap(x, hasCheck, prev);
			if(tbl2$ERR == 1) {
				# ERROR
				cat("\n");
				return(tbl);
			} else if(tbl2$ERR == 2) {
				# END of Table:
				tbl = rbind(tbl, tbl2$Tbl);
				cat("\n");
				return(tbl);
			} else if(tbl2$ERR == 3) {
				next;
			}
			tbl2 = tbl2$Tbl;
			prev = tbl2[nRowCheck, idCheck];
			tbl  = rbind(tbl, tbl2);
			break;
		}
	}
	cat("\n");
	return(tbl);
}

