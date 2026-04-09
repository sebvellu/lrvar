#' @export

tvlrvarmult <- function(tsrs, sgmh, krnl = "ba", band = "and", mvng = TRUE) {
	lgth <- NROW(tsrs)
	#
	if (mvng) {
		strt <- c(rep(1, sgmh - 1), sgmh:lgth - sgmh + 1)
		fnsh <- c(seq_len(lgth - sgmh) + sgmh, rep(lgth, sgmh))
		mval <- lgth
		#if (zero) {
		#	strt <- c(1, strt)
		#	fnsh <- c(sgmh, fnsh)
		#}
	} else {
		indx <- blocks(lgth, sgmh, TRUE)
		strt <- indx$strt
		fnsh <- indx$fnsh
		mval <- indx$mval
	}
	#strt <- 1:(lgth - 2  * sgmt + 1)
	#fnsh <- (2 * sgmt):lgth
	#
	rslt <- tvlrvar(tsrs, strt, fnsh, krnl, band)
	rslt <- c(rslt, list(mval = mval))
	return(rslt)
}
