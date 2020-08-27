#' @importFrom R6 R6Class
#' @import data.table

# Summary measure object that collapses the past history into a finite dimensional vector for each person
# Argument summary_function should be a function that takes a data.table representing the history of a single person,
# where each row corresponds to a time point of the past, and the columns are covariates.
# The returned result should be a vector.
# Note that some elements of the final row may be NA if some covariates are measured at the same present time but are still
# ... part of the past due to time ordering. This should be handled as needed.
# Note this is avoided if the covariates specified by "column_names" all occur at the same time (e.g. are all from the L2 node)
# Data to be summarized must have a column name "id" representing the id for each individual.
# A number of possible summary measure objects are available in constructors below.
# Note that the constructor "make_summary_measure_apply" can be used to safely make
# a summary_measure object. The input function should take a vector of observations of a single covariate over time
# for a single person and spit out a single value.

#' @export
Summary_measure <- R6Class(
  classname = "Summary_measure",
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(column_names, summary_function, name = "Summary", strict_past = F, args_to_pass = NULL){
        # Summary function must return data.table with nrow = 1 ...
       # for self$summarize to work correctly.
        summary_function_wrap <- function(data, time,  ...){
          if(all(is.na(data))) {
            return(data)
          }
          result <- summary_function(data, time, ...)
          if(!is.data.table(result)){
            result <- data.table(matrix(result, nrow =1))
          }
          return(result)
        }
        params <- sl3::args_to_list()
        params$summary_function <- summary_function_wrap
        private$.params <- params
    },
    set_name = function(name){
      private$.params$name <- name
    },
    set_strict_past = function(strict_past){
      private$.params$strict_past <- strict_past
    },
    summarize = function(data, time, add_id = T){

      #data <- private$.process_data(data, time, NULL)
      #ssertthat::assert_that(all(c("id", "t") %in% colnames(data)), msg = "Error: Column 'id' or 't' not found in data.")
      if(!is.data.table(data)){
        data = as.data.table(data)
      }

      if(self$params$strict_past) {
        data <- data[which(data$t < time),]
      } else {
        data <- data[which(data$t <= time),]
      }

      func <- private$.params$summary_function
      # Needed since pass by promise would break next line apparently



      reduced_data <- data[,func(.SD, time, self$params$args_to_pass), by = id,
                           .SDcols = self$params$column_names]

    # This code isn't needed unless func does not return a data.table, which can't happen.
     #  num_sample <- length(unique(reduced_data$id))
     #  num_summary_vars <- nrow(reduced_data) / num_sample
     #  reduced_data$summary_id <- c(1:num_summary_vars, num_sample)
     #  reduced_data <- reshape(reduced_data, idvar = "id", timevar = "summary_id", direction = "wide")

       assertthat::assert_that(is.null(self$params$name) | ncol(reduced_data)-1 == length(self$params$name),
                              msg = "The summary measure names does not match length of summary measure function output.")
      if(!is.null(self$params$name)){
        colnames(reduced_data) <- c("id", self$params$name)
      }
      if(!add_id){
        reduced_data$id = NULL
      }
      return(reduced_data)
    }
  ),
  active = list(
    column_names = function(){
      self$params$column_names
    },
    name = function(){
      self$params$name
    },
    strict_past = function(){
      self$params$strict_past
    },
    params = function(){
      private$.params
    }
  ),
  private = list(
    .params = NULL,
    .process_data = function(data, time, row_index = NULL){
      assertthat::assert_that(all(c("id", "t") %in% colnames(data)), msg = "Error: Column 'id' or 't' not found in data.")
      if(!is.data.table(data)){
        data = as.data.table(data)
      }

      if(self$params$strict_past) {
        data <- data[which(data$t < time), ]
      } else {
        data <- data[which(data$t <= time), ]
      }

      if(is.null(row_index)){
        return(data)
      }
      return(data[row_index,])
    }
  )
)

make_summary_measure_NULL <- function(column_names = ""){
  name =  NULL
  summary_function <- function(data,...){
    return(data.table(NULL))

  }

  return(Summary_measure$new(column_names, summary_function, name))
}


make_summary_measure_FULL <- function(column_names){
  column_names <- union("t", column_names)
  name =  NULL

  summary_function <- function(data,...){
    t <- data$t
    data$t <- NULL
    data <- do.call(cbind, lapply(1:ncol(data), function(i){
      dat <- data.table::transpose(data[,i, with =F])
      colnames(dat) <- paste(colnames(data)[i], t, sep = "_")

      return(dat)}))

    return(data)

  }

  return(Summary_measure$new(column_names, summary_function, name))
}

make_summary_measure_baseline <- function(column_names){
  name = paste( column_names, "baseline", sep = "_")

  summary_function <- function(data,...){
    return(first(data))
  }

  return(Summary_measure$new(column_names, summary_function, name))
}


make_summary_measure_last_value <- function(column_names, strict_past = F){
  name = paste(column_names, "most_recent", sep = "_")

  summary_function <- function(data, time,...){

    if(!all.equal(colnames(data), column_names)){
      if(!(all(column_names %in% colnames(data)))){
        stop("Summary function error: Not all column names found in data object.")
      }
      data <- data[, column_names, with = F]
    }

    last_vals <- data[nrow(data),]

    if(length(which(is.na(last_vals)))>0){
    last_vals[,which(is.na(last_vals)),with=F]  <- data[nrow(data)-1, which(is.na(last_vals)), with = F]
    }

    return(last_vals)

  }
  most_recent <-  function(v){v[length(v)]}
  return(make_summary_measure_apply(column_names,  most_recent, strict_past))
  return(Summary_measure$new(column_names, summary_function, name))
}


make_summary_measure_apply <- function(column_names, FUN, strict_past = T){
  name = as.character(substitute(FUN))
  if(name[1] == "function")
  {
    name = "FUN"
  }

  wrap_FUN <- function(v){
    FUN(as.vector(na.omit(v)))
  }
  summary_function <- function(data,...){
    if(!all.equal(colnames(data), column_names)){
      if(!(all(column_names %in% colnames(data)))){
        stop("Summary function error: Not all column names found in data object.")
      }
      data <- data[, column_names, with = F]
    }
    data <- as.data.table(t(apply(data, 2, wrap_FUN)))
    colnames(data) <- as.character(1:ncol(data))

    return(data)

  }
  return(Summary_measure$new(column_names, summary_function, paste(column_names, name, sep = "_"),strict_past))

}


make_summary_measure_running_average <- function(column_names){
  name = paste(column_names, "avg", sep = "_")
  return(make_summary_measure_apply(column_names, mean))
}

make_summary_measure_running_median <- function(column_names){
  name = paste(column_names, "median", sep = "_")
  return(make_summary_measure_apply(column_names,name ))
}
make_summary_measure_relative_difference_from_t0 <- function(column_names){
  summary_function <- function(data,...){
    if(!all.equal(colnames(data), column_names)){
      if(!(all(column_names %in% colnames(data)))){
        stop("Summary function error: Not all column names found in data object.")
      }
      data <- data[, column_names, with = F]
    }
    diff <- data[nrow(data),] - data[1,]
    change <- which(is.na(diff))
    diff <- data[nrow(data)-1,change, with = F] - data[1,change, with = F]
    return(diff)

  }
  rel_diff_t0 <-  function(v){v[length(v)] - v[1]}
  return(make_summary_measure_apply(column_names,  rel_diff_t0))
  return(Summary_measure$new(column_names, summary_function, paste(column_names, "rel_diff_t0")))

}

make_summary_measure_relative_difference_from_last_t <- function(column_names){
  name = paste(column_names, "rel_diff_last_t", sep = "_")

  summary_function <- function(data,...){
    if(!all.equal(colnames(data), column_names)){
      if(!(all(column_names %in% colnames(data)))){
        stop("Summary function error: Not all column names found in data object.")
      }
      data <- data[, column_names, with = F]
    }
    data <- data[nrow(data) - data[nrow(data)-1,],]
    change <- which(is.na(diff))
    data[,change,with=F] <- data[nrow(data)-1,change, with = F] - data[nrow(data)-2,change, with = F]
  }
  rel_diff_last_t <-  function(v){v[length(v)] - v[length(v)-1]}
  return(make_summary_measure_apply(column_names,  rel_diff_last_t))
  return(Summary_measure$new(column_names, summary_function, name))

}

# takes competing risk columns and returns indicator variable if at risk
make_summary_measure_competing_indicator <- function(competing, strict_past = T){
  column_names <- c(competing)
  name <- paste(paste(competing, collapse = "_"), "at_risk", sep = "_")
  summary_function <- function(data,...){
    if(!all.equal(colnames(data), column_names)){
      if(!(all(column_names %in% colnames(data)))){
        stop("Summary function error: Not all column names found in data object.")
      }
      data <- data[, column_names, with = F]
    }
    # If any of competing risks jumped then at_risk is 0
    at_risk <- as.numeric(all(rowSums(data[,competing, with = F])==0))
    return(at_risk)
  }
  return(Summary_measure$new(column_names, summary_function, name, strict_past))

}

make_summary_measure_slope <- function(column_names){
  name = paste(column_names, "slope_in_t", sep = "_")

  summary_function <- function(data,...){

    if(!all.equal(colnames(data), column_names)){
      if(!(all(column_names %in% colnames(data)))){
        stop("Summary function error: Not all column names found in data object.")
      }
      data <- data[, column_names, with = F]
    }
   if("t" %in% colnames(data)){
     t = data[,"t",with = F]
   }
    else{
      t = 1:nrow(data)
      data = cbind(t, data)
    }


    slopes = sapply(colnames(data)[-1], function(name){
      return(as.vector(coef(lm(as.formula(paste(name, "~ t")), data.frame(data)[, c("t", "name")]))[2]))
    })


  }
  return(Summary_measure$new(column_names, summary_function, name))
}
