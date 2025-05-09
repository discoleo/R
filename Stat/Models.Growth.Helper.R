
### Src
# Helper Functions
# for Models.Growth.R


####################

### Helper Functions

library(demodelr)
library(dplyr)
library(ggplot2)

### Init:
xlim = c(0, 15); ylim = c(0, 2.1);

### Palette
col.blue    = c("#0000FFC0", "#72A2FAA0", "#A888FAA0", "#88A8FAA0");
col.magenta = c("#F000FFC0", "#F8A2A2A0", "#F864B4A0", "#F8B4FFA0");
col.green   = c("#00FF00A0", "#90FE64A0", "#32FC90A0", "#64FC90A0");
col.brown   = c("#F8B464A0");
col.grey    = c("#909090B0");
# c("#FF0000A0", "#FC2490A0", "#FC9064A0", "#FE8090A0")


### Eval Formulas:
eval.Fun = function(x, nm, FUN, params) {
	sapply(x, function(x) {
		params[[nm]] = x;
		do.call(FUN, params);
	})
}
eval.Exp = function(x, nm, e, params) {
	sapply(x, function(x) {
		params[[nm]] = x;
		eval(e, params);
	})
}

solve.eq = function(eq, params, Vx = Vmax / 2, Vmax = 2, ..., lim = c(0, 100)) {
	params = c(as.list(params), Vmax=Vmax);
	eqf = function(x)
		eval(eq, c(params, t=x)) - Vx;
	uniroot(eqf, lim, ...);
}

### Plot:
# xlim = defined globally;
curve.ref = function(col = "#F02424E2", lwd = 2, params = NULL, Vmax = 2,
		col.h = "#FF3224B0", xlab = "Time", ylab = "Growth") {
	ylim = c(0, Vmax + 0.1);
	if(is.null(params)) params = list(Vmax = Vmax, b = 1, n = 1);
	curve(eval.Exp(x, "t", MM.eq, params), col=col, lwd=lwd,
		xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab);
	abline(h = Vmax, col = col.h, lwd=lwd);
}
plot.curve = function() {
	# TODO
}


### Michaelis-Menten Models

MM.eq = expression(Vmax * t^n / (b + t^n));

# Vx = solve V - Vx == 0;
curve.MM = function(col, b = c(1, 1/2, 2, 3), n = 1, lwd = 2, legend = TRUE,
		Vmax = 2, Vx = Vmax / 2, xy.txt = c(8, 0.5)) {
	doS = (Vx != 0);
	params = list(Vmax = Vmax, b = b[1], n = n[1]);
	if(doS) sol1 = solve.eq(MM.eq, params, Vx=Vx);
	# curve(eval.Exp(x, "t", MM.eq, params), col = col[1], xlim=xlim, ylim=ylim)
	curve.ref(params = params, col = col[1], lwd=lwd);
	if(length(b) == 1) return(); # Hack: Reference curve;
	params$b = b[2];
	if(doS) sol2 = solve.eq(MM.eq, params, Vx=Vx);
	curve(eval.Exp(x, "t", MM.eq, params), add = T, col = col[2], lwd=lwd)
	params$b = b[3];
	if(doS) sol3 = solve.eq(MM.eq, params, Vx=Vx);
	curve(eval.Exp(x, "t", MM.eq, params), add = T, col = col[3], lwd=lwd)
	params$b = b[4];
	if(doS) sol4 = solve.eq(MM.eq, params, Vx=Vx);
	curve(eval.Exp(x, "t", MM.eq, params), add = T, col = col[4], lwd=lwd)
	if(legend) {
		legend(xy.txt[1], xy.txt[2], legend = paste0("MM", seq_along(b)), fill = col);
	}
	if(doS) {
		sol = rbind(sol1, sol2, sol3, sol4);
		return(sol);
	}
}
curve.MMextn = function(b, n, col, Vx = Vmax / 2, Vmax = 2, lwd = 2, xy.legend) {
	params  = list(Vmax = Vmax, b = b[1], n = n[1]);
	doSolve = (Vx != 0);
	len = length(b);
	if(len == 0) return();
	sol = list();
	for(id in seq(len)) {
		params$b = b[id]; params$n = n[id];
		curve(eval.Exp(x, "t", MM.eq, params), add = TRUE,
			col = col[id], lwd=lwd);
		if(doSolve) sol[[id]] = solve.eq(MM.eq, params, Vmax=Vmax, Vx=Vx);
	}
	xy = xy.legend;
	legend(xy[1], xy[2], legend = paste0("MMv", seq(len)), fill = col);
	if(doSolve) {
		sol = do.call(rbind, sol);
		return(sol);
	}
}


### Saturated Exponential Models

exps.eq = expression(Vmax * (1 - exp(-b*t^n))^k)

curve.ExpS = function(b = 1, n = 1, k = 1,
		Vx = Vmax / 2, Vmax = 2, col = 1, lwd = 2,
		labels = NULL, xy.labels = NULL) {
	len = c(length(b), length(n), length(k));
	if(any(len == 0)) return();
	len = max(len);
	# Params:
	if(length(b) == 1 && len > 1) b = rep(b, len);
	if(length(n) == 1 && len > 1) n = rep(n, len);
	if(length(k) == 1 && len > 1) k = rep(k, len);
	params = list(Vmax = Vmax, b = b[1], n = n[1], k = k[1]);
	sol    = list();
	if(length(col) == 1 && len > 1) col = rep(col, len);
	doSolve = (Vx != 0);
	# Curves:
	for(id in seq(len)) {
		params$b = b[id]; params$n = n[id]; params$k = k[id];
		curve(eval.Exp(x, "t", exps.eq, params), add = TRUE,
			col = col[id], lwd=lwd);
		if(doSolve) sol[[id]] = solve.eq(exps.eq, params, Vmax=Vmax, Vx=Vx);
	}
	# Legend:
	if(! is.null(labels)) {
		xy = xy.labels;
		if(is.null(xy)) xy = c(11.5, Vmax / 4); # Hardcoded!
		legend(xy[1], xy[2], labels, fill=col);
	}
	if(doSolve) {
		sol = do.call(rbind, sol);
		return(sol);
	}
}

# Mix of various Saturated Exponentials
curve.ExpSMix = function(b = c(1, 1, 1/2, 1/5), n = c(1, 1/2, 1, 1), k = 1,
		Vx = Vmax / 2, Vmax = 2, col = 1, lwd = 2,
		xy.labels = c(11.5, Vmax / 4)) {
	doLegend = ! is.null(xy.labels);
	lbls = if(doLegend) {
		paste0("ES", rep(c("n","b"), each=2), 1:2);
	} else NULL;
	curve.ExpS(b=b, n=n, k=k, Vx=Vx, Vmax=Vmax, col=col, lwd=lwd,
		labels = lbls, xy.labels = xy.labels);
}


### ATAN Model

atan.eq = expression(Vmax * (atan(b * t^n) * 2/pi)^k)

curve.Atan = function(b = 1, n = 1, k = 1,
		Vx = Vmax / 2, Vmax = 2,
		col = 1, lwd = 2, labels = NULL, xy.labels = NULL) {
	len = c(length(b), length(n), length(k));
	if(any(len == 0)) return();
	len = max(len);
	# Params:
	if(length(b) == 1 && len > 1) b = rep(b, len);
	if(length(n) == 1 && len > 1) n = rep(n, len);
	if(length(k) == 1 && len > 1) k = rep(k, len);
	params = list(Vmax=Vmax, b = b[1], n = n[1], k = k[1]);
	sol    = list();
	if(length(col) == 1 && len > 1) col = rep(col, len);
	doSolve = (Vx != 0);
	for(id in seq(len)) {
		params$b = b[id]; params$n = n[id]; params$k = k[id];
		curve(eval.Exp(x, "t", atan.eq, params), add = TRUE,
			col = col[id], lwd=lwd);
		if(doSolve) sol[[id]] = solve.eq(atan.eq, params, Vmax=Vmax, Vx=Vx);
	}
	# Legend:
	if(! is.null(xy.labels)) {
		if(is.null(labels)) {
			labels = paste0("AT", seq(len));
		} else if(length(labels == 1) && len > 1) {
			labels = paste0(labels, seq(len));
		}
		xy = xy.labels;
		legend(xy[1], xy[2], labels, fill=col);
	}
	if(doSolve) {
		sol = do.call(rbind, sol);
		return(sol);
	}
}

### Integrals

### I( x^p / (x^n + 1)^k )
frI = function(up, Vmax, p, n, k=1) {
	FUN = \(x) x^p / (x^n + 1)^k;
	tmp = integrate(FUN, lower = 0, upper = up, rel.tol = 1E-8)$value;
	# Vmax * sin(pi*(p+1)/n) * n / pi * tmp; # k == 1;
	scale = gamma((p+1)/n) * gamma(k - (p+1)/n) / gamma(k) / n;
	Vmax * tmp /scale;
}

curve.IntFr = function(p = c(0, 1/3, 1/2, 4/3), n = 3, k = 1,
		col = 1, lwd = 2, Vx = Vmax / 2, Vmax = 2,
		labels = NULL, xy.labels = NULL) {
	len = c(length(p), length(n), length(k));
	if(any(len == 0)) return();
	if(Vx != 0) doSolve = TRUE;
	len = max(len);
	sol = list();
	eq  = expression(frI(up=t, Vmax=Vmax, p=p, n=n, k=k));
	# Params:
	if(length(p) == 1 && len > 1) b = rep(b, len);
	if(length(n) == 1 && len > 1) n = rep(n, len);
	if(length(k) == 1 && len > 1) k = rep(k, len);
	params = list(Vmax=Vmax, p = p[1], n = n[1], k = k[1]);
	for(id in seq(len)) {
		params$p = p[id]; params$n = n[id];
		curve(eval.Fun(x, "up", frI, params), col = col[id], lwd=lwd, add = TRUE);
		if(doSolve) sol[[id]] = solve.eq(eq, params, Vmax=Vmax, Vx=Vx);
	}
	# Legend:
	if(! is.null(labels)) {
		xy = xy.labels;
		if(is.null(xy)) xy = c(11.5, Vmax / 4); # Hardcoded!
		legend(xy[1], xy[2], labels, fill=col);
	}
	if(doSolve) {
		sol = do.call(rbind, sol);
		return(sol);
	}
}

### ODE

### Logistic Growth

log.eq = V ~ k * V^p * (Vmax^n - V^n)^m

curve.odeLog = function(k = 1, p = 1, m = 1, n = 1,
		Vmax = 2, V0 = 0.1, col = 1, lwd = 2,
		labels = NULL, xy.labels = NULL,
		n_steps = 160, dt = 0.1) {
	len = c(length(k), length(p), length(n));
	if(any(len == 0)) return(data.frame(t = numeric(0), V = numeric(0)));
	len = max(len);
	# Parameters:
	if(length(k) == 1 && len > 1) k = rep(k, len);
	if(length(p) == 1 && len > 1) p = rep(p, len);
	if(length(n) == 1 && len > 1) n = rep(n, len);
	params = list(Vmax = Vmax, k=k[1], p=p[1], m=m[1], n=n[1]);
	V0 = c(V = V0);
	for(id in seq(len)) {
		params$k = k[id]; params$p = p[id]; params$n = n[id];
		sol = rk4(log.eq, V0, parameters = params, deltaT=dt, n_steps = n_steps);
		lines(sol, col = col[id], lwd=lwd);
	}
	# Legend:
	if(! is.null(xy.labels)) {
		if(is.null(labels)) {
			labels = paste0("OL", seq(len));
		} else if(length(labels) == 1 && len > 1) {
			labels = paste0(labels, seq(len));
		}
		xy = xy.labels;
		legend(xy[1], xy[2], labels, fill=col);
	}
}
# ODE: Logistic Eq with Initial conditions
# - Effect of Initial condition;
curve.odeLogI0 = function(p = c(1, 0.6, 1.25), V0 = c(0.1, 0.005),
		k = 1, n = 1, lwd = 2,
		xy1 = c(8, 0.5), xy2 = c(11.5, 0.5), adj.xy2 = c(0.25, 0.05)) {
	xy1 = xy1 + V0[1]; xy2 = xy2 + V0[1];
	curve.ref(col.grey[1], lwd=lwd);
	#
	col = col.magenta; lbls = "OL";
	curve.odeLog(k = k, p = p, n = n, V0=V0[1], col=col, lwd=lwd,
		xy=xy1, labels = lbls);
	abline(h = V0, lty = 2, col = "#B06464B2")
	# Initial condition
	col = col.blue; lbls = "OL0";
	curve.odeLog(k = k, p = p, n = n, V0=V0[2], col=col, lwd=lwd,
		xy=xy2, labels = lbls);
	xy2 = xy2 + adj.xy2;
	lblV0 = paste0("V0 = ", V0[2]);
	text(xy2[1], xy2[2], lblV0, col = col[1], adj = c(0, 0));
}

# Variation of Power n
curve.odeLogPn = function(n = c(1, 1/2), p = c(1,1/3,1.5), k = 1,
		V0 = 0.1, Vmax = 2,
		lwd = 2, xy1 = c(8, 0.5), xy2 = c(11.5, 0.5), xy.adj = c(0.5, 0.05)) {
	curve.ref(col.grey[1], lwd=lwd, Vmax=Vmax); # Ref
	xy1 = xy1 + V0; xy2 = xy2 + V0;
	#
	col = col.green;
	curve.odeLog(k = k, p = p, n = n[1], V0=V0, Vmax=Vmax,
		col=col, lwd=lwd, xy=xy1, labels = "OL");
	# Variant: n
	if(length(n) == 1) return();
	col = col.magenta;
	curve.odeLog(k = k, p = p, n = n[2], V0=V0, Vmax=Vmax,
		col=col, lwd=lwd, xy=xy2, labels = "OLn");
	xy2  = xy2 + xy.adj;
	lbln = paste0("n = ", n[2]);
	text(xy2[1], xy2[2], lbln, col = col[1], adj = c(0, 0))
	#
	abline(h = V0, lty = 2, col = "#B06464B2")
}

