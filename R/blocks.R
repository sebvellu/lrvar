blocks <- function(lgth, sgmh, full = TRUE) {
	if (full) {
		lval <- seq(1, floor(lgth/sgmh))
		fnsh <- pmin(lval * sgmh + sgmh, lgth)
	} else {
		lval <- seq(1, floor(lgth/sgmh) - 1)
		fnsh <- lval * sgmh + sgmh
	}
	strt <- lval * sgmh - sgmh + 1
	mval <- length(strt)
	return(list(strt = strt, fnsh = fnsh, mval = mval))
}
