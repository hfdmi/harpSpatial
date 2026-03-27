harpSpatial_neighborhood_contingency_tables <- function(obfield=NULL, fcfield=NULL, thresholds=NULL, scales=0, probs=seq(0, 1, 0.1)) {

    nprob <- length(probs)
    for (threshold in thresholds) {
        cumsum_ob <- cumsum(obfield, threshold)
        for (scale in scales) { #obmask=1 if threshold is exceeded somewhere in neighbourhood
            obmask <- window_mask_from_cumsum(cumsum_ob, scale)
            for (iprob in 1:nprob) {
                prob <- probs[iprob]
                fcmask <- 1  * (fcfield >= prob)
                hit <- sum(fcmask * obmask)
                fa <- sum(fcmask * (1-obmask))
                miss <- sum((1-fcmask) * obmask)
                cr <- sum((1-fcmask) * (1-obmask))
                if (iprob == 1) { #First entry in contingency_table
                    contingency_table <- data.frame(
                        "threshold" = threshold,
                        "scale" = scale,
                        "prob" = prob,
                        "hit" = hit,
                        "fa" = fa,
                        "miss" = miss,
                        "cr" = cr
                    )
                } else { #Subsequent entries in contingency_table
                    contingency_table[nrow(contingency_table)+1,] <- c(threshold, scale, prob, hit, fa, miss, cr)
                }
            }
        }
    }
    return(contingency_table)
}

cumsum <- function(indat, threshold) {
  outdat <- 1 * (indat >= threshold) #outdat is a geofield with dimensions of indat
  Ni <- attributes(indat)$dim[[1]]
  Nj <- attributes(indat)$dim[[2]]
  for (j in 1:Nj) {
    for (i in seq(2,Ni)) {
      outdat[i,j] <- 1 * (indat[i,j] >= threshold) + outdat[i-1,j]
    }
  }
  for (i in 1:Ni) {
    for (j in seq(2,Nj)) {
      outdat[i,j] <- outdat[i,j] + outdat[i,j-1]
    }
  }
  return(outdat)
}

window_mask_from_cumsum <- function(indat, r) {
  outdat <- 0 * (indat == 0)      #outdat is a geofield with dimensions of indat
  Ni <- attributes(indat)$dim[[1]]
  Nj <- attributes(indat)$dim[[2]]
  for (i in 1:Ni) {
    imax <- min(i+r,Ni)
    for (j in 1:Nj) {
      jmax <- min(j+r,Nj)
      outdat[i,j] <- indat[imax,jmax]
      if (i-1 > r) {
    outdat[i,j] <- outdat[i,j] - indat[i-r-1,jmax]
    if (j-1 > r) {
      outdat[i,j] <- outdat[i,j] + indat[i-r-1,j-r-1] - indat[imax,j-r-1]
    }
      }
      else if (j-1 > r) {
    outdat[i,j] <- outdat[i,j] - indat[imax,j-r-1]
      }
    }
  }
  outdat <- 1 * (outdat > 0)
}
