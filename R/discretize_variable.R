#' @export
discretize_variable <- function(data, variable, num_cats, breakpoints=NULL){
  vals <- unlist(data[,variable, with = FALSE])
  if(!is.factor(vals)){
    if(is.null(breakpoints)){
      quants <- seq(from=0,to=1,length=num_cats+1)  
      breakpoints <- quantile(vals, quants)
      
      # todo: warn on reducing breakpoint number
      breakpoints <- unique(breakpoints)
    }
    
    
    indexes <- seq_len(length(breakpoints)-1)
    # manual labels because cut too easily defaults to scientific notation
    left_bracket <- "["
    right_bracket <- ifelse(indexes < (length(breakpoints)-1), ")", "]")
    formatted <- format(breakpoints, trim=TRUE, digits=2)
    lb <- formatted[indexes]
    ub <- formatted[indexes+1]
    labels <- paste(left_bracket, lb, ",", ub, right_bracket, sep="")

    #todo: format labels more cleanly
    as_factor <- cut(vals, breakpoints, labels=labels, right = FALSE, include.lowest = TRUE)
    set(data, , variable, as_factor)
  }
}