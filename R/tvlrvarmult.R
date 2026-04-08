#' @export

tvlrvarmult <- function(
		tsrs, sgmh, slct = 1, krnl = "ba", band = "and"#, zero = FALSE #segment-half-width
) {
	lgth <- NROW(tsrs)
	#
	if (slct == 1) {
		indx <- blocks(lgth, sgmh, TRUE)
		strt <- indx$strt
		fnsh <- indx$fnsh
		mval <- indx$mval
	} else { #if (slct == 3) {
		strt <- c(rep(1, sgmh - 1), sgmh:lgth - sgmh + 1)
		fnsh <- c(seq_len(lgth - sgmh) + sgmh, rep(lgth, sgmh))
		mval <- lgth
		#if (zero) {
		#	strt <- c(1, strt)
		#	fnsh <- c(sgmh, fnsh)
		#}
	}
	#strt <- 1:(lgth - 2  * sgmt + 1)
	#fnsh <- (2 * sgmt):lgth
	#
	rslt <- tvlrvar(tsrs, strt, fnsh, krnl, band)
	rslt <- c(rslt, list(mval = mval))
	return(rslt)
}
