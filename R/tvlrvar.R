tvlrvar <- function(
		tsrs, strt, fnsh, krnl = "ba", band = "and"
) {
	lgth <- NROW(tsrs)
	nums <- NCOL(tsrs)
	#
	lrbw <- bwclc
	#bwdh <- match.fun(bwclc)(tsrs, krnl, band, bfrc, NULL)
	lrwg <- lrwghts
	#
	stfn <- fnsh - strt + 1
	#
	if (nums == 1) {
		long <- numeric(length(stfn))
		shrt <- long
		half <- long
		delt <- long
		fixb <- long
		#
		comp <- computelrvec
		#
		for (h in seq_along(stfn)) {
			locl <- tsrs[strt[h]:fnsh[h]]
			#
			bwdh <- lrbw(locl, krnl, band, NULL)
			wght <- lrwg(stfn[h], bwdh, krnl)
			lvar <- comp(locl/sqrt(stfn[h]), wght)
			#
			long[h] <- lvar$long
			shrt[h] <- lvar$shrt
			half[h] <- lvar$half
			delt[h] <- lvar$delt
			fixb[h] <- bwdh/stfn[h]
		}
	} else {
		long <- array(0, c(nums, nums, length(stfn)))
		shrt <- long
		half <- long
		delt <- long
		fixb <- numeric(length(stfn))
		#
		comp <- computelrmat
		#
		for (h in seq_along(stfn)) {
			locl <- tsrs[strt[h]:fnsh[h], , drop = FALSE]
			#
			bwdh <- lrbw(locl, krnl, band, NULL)
			wght <- lrwg(stfn[h], bwdh, krnl)
			lvar <- comp(locl/sqrt(stfn[h]), wght)
			#
			long[, , h] <- lvar$long
			shrt[, , h] <- lvar$shrt
			half[, , h] <- lvar$half
			delt[, , h] <- lvar$delt
			fixb[h] <- bwdh/stfn[h]
		}
	}
	#
	return(list(
		fixb = fixb,
		shrtvar = shrt,
		halfvar = half,
		longvar = long,
		deltvar = delt
	))
}
