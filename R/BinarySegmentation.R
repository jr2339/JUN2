

#' Title
#'
#' @param data observations
#' @param lower lower open interval bound
#' @param upper upper open interval bound
#'
#' @return vector containing all CUSUM values for all points in the interval
#' @export

bs_stat <- function(data, lower, upper){

  test.stat <- NULL
  cx <- cumsum(data)
  for (i in (lower+1):(upper-1)){

    #for each k in the interval we create a test score,
    #|bs_stat|

    coeff <- (i - lower)/(upper - lower)
    max_difference <- (cx[upper]-cx[lower])
    const <- coeff*max_difference

    top_fraction <- const - (cx[i]-cx[lower])
    bottom_fraction <- sqrt((i - lower)*(1-(coeff)))
    test.stat <- c(test.stat, abs((top_fraction)/(bottom_fraction)))

  }
  return(test.stat)
}



#' Title
#'
#' @param data observations
#' @param interval intervals to partition observations
#'
#' @return a list of est.signal and est.means
#' @export

bs.model <- function(data, interval){


  est.signal <- NULL
  means <- NULL

  for(i in 1:(length(interval)-1)){
    a <- interval[i]
    b <- interval[i+1]
    est.signal <- c(est.signal, rep(mean(data[a:b]), (b-a)))

    means <- c(means, mean(data[a:b]))
  }
  est.signal <- c(est.signal[1],est.signal)
  return(list(est.signal=est.signal, est.means=means))
}





#' Title
#'
#' @param data observations
#' @param former.segments intervals to be tested
#' @param segments critical value used in testing
#'
#' @return finsal.segments (intervals to test in next iteration)
#' @export

bs.helper <- function(data, former.segments = NULL,segments){
  changepoints <- NULL
  est.changepoints <- NULL
  finsal.segments <- NULL

  if(is.null(former.segments)){
    current.segment <- c(1, segments)
    finsal.segments <- c(1,segments)
  } else {
    current.segment <- sort(former.segments)
    finsal.segments <- former.segments
  }



  #we need to order the estimated change-points so that the largest bs_stat is first
  for (i in 1:(length(current.segment)-1)){ # n changepoints gives n-1 segments
    #we take the maximum as a estimated change point
    # the zero removes the error message, but has no affect on the output of the function
    if(max(c(bs_stat(data, current.segment[i],current.segment[i+1]), 0), na.rm = TRUE)>= 0.05){
      ch <- which.max(bs_stat(data, current.segment[i],current.segment[i+1]))+current.segment[i]
      #we test if the changepoinit is already accounted for
      if(!( ch %in% former.segments)){

        est.changepoints <- c(est.changepoints, ch)
        changepoints <- c(changepoints, max(bs_stat(data, current.segment[i],current.segment[i+1]), na.rm = TRUE))
      }
    }
  }

  if(!(is.null(est.changepoints))){

    bound <- cbind(est.changepoints, changepoints)

    est.changepoints <- bound[order(-bound[,2]), 1]

  }
  finsal.segments <- c(finsal.segments, est.changepoints)
  return(unname(finsal.segments))
}



#' Title
#'
#' @param data datae data
#' @param segments the number of segments
#'
#' @return {estimated signal, estimated means, estimated changepoints}
#' @export

BINSEG <- function(data,segments){

  held_segments <- 0
  former.segments <- NULL
  est.signal <- NULL
  former.segments <- bs.helper(data,NULL,segments)

  while(!setequal(held_segments, former.segments)){
    # binary segmentation runs until the operation doesn't produce a new change-point
    held_segments <- former.segments
    former.segments <- bs.helper(data, held_segments,segments)
  }
  #print("in order of occurence the estimated changepoints are ")
  ordered_finsal.segments <- sort(former.segments)
  print(ordered_finsal.segments)
  output <- bs.model(data, ordered_finsal.segments)
  interval <- list(est.changepoints = ordered_finsal.segments)
  output <- c(output, interval)

  return(output)
}

