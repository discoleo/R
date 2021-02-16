
### package transport
### Improvement

transparent.col = function(col, min=40) {
	len = length(col);
	if(len == 1) { len = col; col = heat.colors(len); }
	apply(
		rbind(col2rgb(col),
			alpha=seq(min, 255, length.out=len))/255, 2,
		function(x) rgb(x[1], x[2], x[3], x[4]))
}

plot.pgrid = function (x, y = NULL, tplan = NULL, col,
	mass = c("colour", "thickness"), length = 0.1, acol, bcol = 4, lwd, rot = FALSE, 
    overlay = FALSE, static.mass = TRUE, ...) 
{
    stopifnot(class(x) == "pgrid")
    if (class(y) == "wpp") {
        if (missing(acol)) {
            acol <- "#996699"
        }
        if (missing(lwd)) {
            lwd <- 1.5
        }
        return(plot_pgrid_wpp(x, y, tplan, pmass = TRUE, cex = 0.8, 
            length = length, acol = acol, col = bcol, lwd = lwd, 
            rot = TRUE, ...))
    }
    a <- x
    if (a$dimension != 2) 
        stop("plot.pgrid is currently only implemented for 2-d grids")
    xi <- a$generator[[1]]
    eta <- a$generator[[2]]
    if (missing(y)) {
        image2(xi, eta, a$mass, rot = rot, col = grey(0:200/200), 
            asp = 1, xlab = "", ylab = "")
        invisible()
    }
    else {
        stopifnot(class(y) == "pgrid")
        b <- y
        if (b$dimension != 2) 
            stop("plot.pgrid is currently only implemented for 2-d grids")
        if (!(a$structure %in% c("square", "rectangular")) || 
            !(b$structure %in% c("square", "rectangular"))) 
            stop("transport.pgrid is currently only implemented for rectangular pixel grids")
        if (missing(tplan)) {
            if (!overlay) {
                par(mfrow = c(1, 2))
                image2(xi, eta, a$mass, rot = rot, col = grey(0:200/200), 
                  asp = 1, xlab = "", ylab = "", 
                  ...)
                image2(b$generator[[1]], b$generator[[2]], b$mass, 
                  rot = rot, col = grey(0:200/200), asp = 1, 
                  xlab = "", ylab = "", ...)
                par(mfrow = c(1, 1))
            }
            else {
                stopifnot(compatible(a, b))
                zeta <- a$mass - b$mass
                image2(xi, eta, zeta, rot = rot, col = grey(0:200/200), 
                  asp = 1, xlab = "", ylab = "", 
                  ...)
            }
            invisible()
        }
        else {
            stopifnot(compatible(a, b))
            mass <- match.arg(mass)
            if (missing(acol)) 
                acol <- 2
            if (missing(lwd)) 
                lwd <- 3
            gg <- expand.grid(xi, eta)
            maxmass <- max(tplan[, 3])
            zeta <- a$mass - b$mass
            if (rot) {
                yco <- rev(gg[, 1])
                xco <- gg[, 2]
            }
            else {
                xco <- gg[, 1]
                yco <- gg[, 2]
            }
            transport:::image2(xi, eta, zeta, rot = rot, col = grey(0:200/200), 
                asp = 1, axes = FALSE, xlab = "", ylab = "", 
                ...)
            if (mass == "colour") {
                stplan <- tplan[order(tplan[, 3]), ]
				if(missing(col)) {
					cc = heat.colors(128)
				} else {
					cc = col;
				}
				col.len = length(cc);
                arrcols <- cc[as.numeric(as.character(cut(stplan[, 3],
					breaks = seq(0, maxmass, length.out = 1 + col.len),
					labels = 1:col.len)))]
                wh <- which(stplan$from != stplan$to)
                arrows(xco[stplan$from[wh]], yco[stplan$from[wh]], 
                  xco[stplan$to[wh]], yco[stplan$to[wh]], angle = 5, 
                  length, col = arrcols[wh], lwd = lwd)
                if (static.mass == TRUE) {
                  nwh <- (1:dim(stplan)[1])[-wh]
                  points(xco[stplan$from[nwh]], yco[stplan$from[nwh]], 
                    pch = 16, cex = 0.5 + lwd * 0.1, col = arrcols[nwh])
                  points(xco[stplan$from[nwh]], yco[stplan$from[nwh]], 
                    cex = 0.5 + lwd * 0.1)
                }
            }
            else {
                wh <- which(tplan$from != tplan$to)
                arrows(xco[tplan$from[wh]], yco[tplan$from[wh]], 
                  xco[tplan$to[wh]], yco[tplan$to[wh]], angle = 5, 
                  length, col = acol, lwd = (8 * tplan[, 3]/maxmass)[wh])
                if (static.mass == TRUE) {
                  nwh <- (1:dim(tplan)[1])[-wh]
                  points(xco[tplan$from[nwh]], yco[tplan$from[nwh]], 
                    pch = 16, cex = 0.5 + (8 * tplan[, 3]/maxmass)[nwh] * 
                      0.1, col = acol)
                  points(xco[tplan$from[nwh]], yco[tplan$from[nwh]], 
                    cex = 0.5 + (8 * tplan[, 3]/maxmass)[nwh] * 
                      0.1)
                }
            }
            invisible()
        }
    }
}
