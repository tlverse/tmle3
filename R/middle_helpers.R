#' Helper Functions for Longitudinal Mediation
#'
#' Handles the ARZL(Y) structure
#'
#' @param data a \code{data.frame}, or \code{data.table} containing data for use in estimation
#' @param node_list a list of character vectors, listing the variables that comprise each node
#' @param variable_types a list of variable types, one for each node. If missing, variable types will be guessed
#' @param tmle_task a \code{\link{tmle3_Task}} as constructed via \code{point_tx_task}
#' @param learner_list a list of sl3 learners, one for A and one for Y to be used for likelihood estimation
#' @param ... extra arguments.
#' @export
#' @rdname point_tx
middle_npsem <- function(node_list, variable_types = NULL) {
  # make tmle_task
  npsem <- c(define_node("L_0", node_list$L_0, variable_type = variable_types$L_0),
             lapply(2:length(node_list), function(k) {
               if (k < length(node_list)) {
                 define_node(names(node_list)[k],
                             node_list[[k]],
                             names(node_list)[1:(k-1)],
                             variable_type = variable_types[[ names(node_list)[k] ]])
               } else {
                 define_node(names(node_list)[k],
                             node_list[[k]],
                             names(node_list)[1:(k-1)],
                             variable_type = variable_types[[ names(node_list)[k] ]],
                             scale = TRUE)
               }
             })
  )

  # npsem <- list(
  #   define_node("W", node_list$W, variable_type = variable_types$W),
  #   define_node("A", node_list$A, c("W"), variable_type = variable_types$A),
  #   define_node("Y", node_list$Y, c("A", "W"), variable_type = variable_types$Y, scale = TRUE)
  # )

  return(npsem)
}

#' @export
#' @rdname point_tx
middle_task <- function(data, node_list, variable_types = NULL, ...) {
  setDT(data)

  npsem <- middle_npsem(node_list, variable_types)

  if (!is.null(node_list$id)) {
    tmle_task <- tmle3_Task$new(data, npsem = npsem, id = node_list$id, ...)
  } else {
    tmle_task <- tmle3_Task$new(data, npsem = npsem, ...)
  }

  return(tmle_task)
}

#' @export
#' @rdname point_tx
middle_likelihood <- function(tmle_task, learner_list) {
  factor_list <- list()

  # covariates
  L_0_factor <- define_lf(LF_emp, "L_0")
  factor_list[[1]] <- L_0_factor

  # treatment (bound likelihood away from 0 (and 1 if binary))
  # A_type <- tmle_task$npsem[["A"]]$variable_type
  A_type <- tmle_task$npsem[[ grep("A", names(tmle_task$npsem))[1] ]]$variable_type  # evaluate the first A node
  if (A_type$type == "continous") {
    A_bound <- c(1 / tmle_task$nrow, Inf)
  } else if (A_type$type %in% c("binomial", "categorical")) {
    A_bound <- 0.025
  } else {
    A_bound <- NULL
  }

  temp_names <- names(tmle_task$npsem)
  loc_A <- grep("A", temp_names)
  factor_list[loc_A] <- lapply(loc_A, function(k) {
    LF_fit$new(temp_names[k], learner = learner_list[[ temp_names[k] ]], bound = A_bound, type = "density")
  })
  # A_factor <- define_lf(LF_fit, "A", learner = learner_list[["A"]], bound = A_bound)

  # others
  loc_others <- (1:length(temp_names))[-c(grep("A", temp_names), 1)]
  factor_list[loc_others] <- lapply(loc_others, function(k) {
    LF_fit$new(temp_names[k], learner = learner_list[[ temp_names[k] ]], type = "density")
  })

  # # outcome
  # Y_factor <- define_lf(LF_fit, "Y", learner = learner_list[["Y"]], type = "mean")
  #
  # # construct and train likelihood
  # factor_list <- list(W_factor, A_factor, Y_factor)

  likelihood_def <- Likelihood$new(factor_list)
  likelihood <- likelihood_def$train(tmle_task)
  temp <- likelihood$list_all_predicted_lkd
  return(likelihood)
}




#' @export
expand_values <- function(variables, to_drop = NULL, values = NULL, rule_variables = NULL, rule_values = NULL, ...) {
  # variables is the vector of variable names
  # drop is either a vector of which variables to drop, or their indices
  # values are the possible values to expand; if null, generate binary values
  # ... are other value rules

  # input_list <- list(A = 1)

  input_list <- list(...)
  rules_list <- lapply(names(input_list), function(eachName) {
    if (length(grep("_", eachName)) > 0) eachName else variables[grep(eachName, map_chr(variables, ~strsplit(.x, "_")[[1]][1]))]
  })
  if(!is.null(rule_variables)) {
    temp_list <- as.list(rule_values)
    names(temp_list) <- rule_variables
    input_list <- c(input_list, temp_list)
    rules_list <- c(rules_list, rule_variables)
  }

  # drop variables that don't want to generate
  # to_drop <- c("A", "L1_0")
  if (is.null(to_drop)) variables_to_generate <- variables else
    if(is.numeric(to_drop)) variables_to_generate <- variables[-to_drop] else
      if (is.character(to_drop)) {
        ind_to_drop <- map(to_drop, ~which(variables == .x)) %>% compact %>% unlist
        if (length(ind_to_drop) == 0) variables_to_generate <- variables else
          variables_to_generate <- variables[-ind_to_drop]
      } else {
        # consider how to identify error later
        print("invalid variables to drop")
        variables_to_generate <- variables
      }

  all_possible_values <- map(variables_to_generate, function(eachVar) {
    test_rules <- which(map_dbl(rules_list, ~length(grep(eachVar, .x))) != 0)
    if (length(test_rules) == 0) return(0:1) else {
      return(input_list[test_rules] %>% as.numeric)
    }
  }) %>% expand.grid() %>% data.frame
  colnames(all_possible_values) <- variables_to_generate

  # to add values
  # to add variable name rules

  return(all_possible_values)
}

# values <- expand_values(variables, A = 1, rule_variables = last(variables), rule_values = 1)
# merge(values[1:5, ], values)

#' @export
which_variable <- function(variables, target_variable, timepoint) {
  grep(paste0(target_variable, "_", timepoint), variables)
}

#' @export
# get the orders of variables after droping by name or time
which_variable_drop <- function(variables, to_drop_variable = NULL, to_drop_time = NULL, to_drop = NULL) {
  if (!is.null(to_drop_variable)) {
    temp_drop_variable <- lapply(to_drop_variable, function(s) grep(s, map_chr(variables, ~strsplit(.x, "_")[[1]][1]))) %>% unlist
  } else
    temp_drop_variable <- NULL
  if (!is.null(to_drop_time)) {
    temp_drop_time <- lapply(to_drop_time, function(s) grep(s, map_chr(variables, ~strsplit(.x, "_")[[1]][2]))) %>% unlist
  }  else temp_drop_time <- NULL
  if (!is.null(to_drop)) {
    temp_drop_both <- lapply(to_drop, function(s) grep(s, variables)) %>% unlist
  } else temp_drop_both <- NULL
  temp_drop <- c(temp_drop_variable, temp_drop_time, temp_drop_both)
  if (!is.null(temp_drop)) (1:length(variables))[-temp_drop] %>% return else 1:length(variables)
}

# which_variable_drop(variables, "L")

#' @export
# take variables by name or time
which_variable_take <- function(variables, to_take_variable = NULL, to_take_time = NULL, logic = "or") {
  if (!is.null(to_take_variable)) {
    temp_take_variable <- lapply(to_take_variable, function(s) grep(s, map_chr(variables, ~strsplit(.x, "_")[[1]][1]))) %>% unlist
  } else
    temp_take_variable <- NULL
  if (!is.null(to_take_time)) {
    temp_take_time <- lapply(to_take_time, function(s) grep(s, map_chr(variables, ~strsplit(.x, "_")[[1]][2]))) %>% unlist
  }  else temp_take_time <- NULL
  # default is to take all selected, by name or by time
  # and is to take the intersection
  if (logic == "or") temp_take <- c(temp_take_variable, temp_take_time) else
    if (logic == "and") temp_take <- intersect(temp_take_variable, temp_take_time)
    if (!is.null(temp_take)) sort(unique(temp_take)) else NULL
}




#' @export
get_current_newH <- function(loc_node,
                             tmle_task, obs_data,
                             intervention_variables, intervention_levels_treat, intervention_levels_control,
                             list_all_predicted_lkd
) {
  obs_variable_names <- colnames(obs_data)
  temp_node_names <- names(tmle_task$npsem)
  loc_A <- grep("A", temp_node_names)
  loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
  loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
  intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))

  # ind_var is the order among variables
  ind_var <- loc_current_var <- which(obs_variable_names == tmle_task$npsem[[loc_node]]$variables)
  intervention_variables_loc_needed <- intervention_variables_loc[intervention_variables_loc < ind_var]

  # only update t!=0 and non-A nodes
  if (loc_node %in% c(loc_Z, loc_RLY)) {
    temp_current <- list_all_predicted_lkd[[loc_node]]
    Xt_input <- temp_current[ind_var]
    loc_to_update <- Xt_input == 1
    df_to_update <- temp_current[loc_to_update, ]  # the possible input where L_t = 1
    df_to_update_0 <- df_to_update
    df_to_update_0[ind_var] <- 0  # replacing the above input with Lt = 0


    # get Q1 current
    {
      data_temp <- df_to_update[1:ind_var]
      # generate all possible 0, 1 valued rlz; set A to 0 first; set y_tau to always 1;
      all_possible_rlz_1 <- expand_values(obs_variable_names, to_drop = c(1:(ind_var) ),
                                          A = 1,
                                          rule_variables = c(last(obs_variable_names)), rule_values = c(1))
      all_possible_rlz_0 <- expand_values(obs_variable_names, to_drop = c(1:(ind_var) ),
                                          A = 0,
                                          rule_variables = c(last(obs_variable_names)), rule_values = c(1))
      # for each observed L_0 vector, generate all needed combinations, one version for A = 1, one version for A = 0
      unique_input <- data_temp[1:(ind_var)] %>% unique
      library_output <- data.frame(unique_input, output =
                                     map_dbl(1:nrow(unique_input), function(which_row) {
                                       # probs in the integrals, A=1 or A=0 is inserted
                                       temp_all_comb_0 <- temp_all_comb_1 <- unique_input[which_row, ]
                                       if (length(all_possible_rlz_0) > 0) {
                                         temp_all_comb_0 <- suppressWarnings(cbind(temp_all_comb_0, all_possible_rlz_0))
                                         temp_all_comb_1 <- suppressWarnings(cbind(temp_all_comb_1, all_possible_rlz_1))
                                       }
                                       if (length(intervention_variables_loc_needed) > 0) {
                                         for (i in 1:length(intervention_variables_loc_needed)) {
                                           temp_all_comb_0[, intervention_variables_loc[i]] <- intervention_levels_control[i]
                                           temp_all_comb_1[, intervention_variables_loc[i]] <- intervention_levels_treat[i]
                                         }
                                       }
                                       # for all non-A, non-0 variables, calculate the variable by rule
                                       # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
                                       loc_Z_need <- loc_Z[loc_Z > loc_node]  # only integrate previous variables; debug: not including current
                                       temp_list_0 <- lapply(loc_Z_need,
                                                             function(each_t) {
                                                               left_join(temp_all_comb_0, list_all_predicted_lkd[[each_t]])$output
                                                             })
                                       loc_other_needed <- loc_RLY[loc_RLY > loc_node]
                                       temp_list_1 <- lapply(loc_other_needed,
                                                             function(each_t) {
                                                               left_join(temp_all_comb_1, list_all_predicted_lkd[[each_t]])$output
                                                             })
                                       temp_list <- c(temp_list_0, temp_list_1)
                                       pmap_dbl(temp_list, prod) %>% sum %>% return
                                     })
      )
      Q1_current <- left_join(data_temp, library_output)$output
      if (loc_node == length(temp_node_names)) Q1_current <- rep(1, length(Q1_current))
    }

    # get Q0 current
    {
      data_temp <- df_to_update_0[1:ind_var]
      # generate all possible 0, 1 valued rlz; set A to 0 first; set y_tau to always 1;
      all_possible_rlz_1 <- expand_values(obs_variable_names, to_drop = c(1:(ind_var) ),
                                          A = 1,
                                          rule_variables = c(last(obs_variable_names)), rule_values = c(1))
      all_possible_rlz_0 <- expand_values(obs_variable_names, to_drop = c(1:(ind_var) ),
                                          A = 0,
                                          rule_variables = c(last(obs_variable_names)), rule_values = c(1))
      # for each observed L_0 vector, generate all needed combinations, one version for A = 1, one version for A = 0
      unique_input <- data_temp[1:(ind_var)] %>% unique
      library_output <- data.frame(unique_input, output =
                                     map_dbl(1:nrow(unique_input), function(which_row) {
                                       # probs in the integrals, A=1 or A=0 is inserted
                                       temp_all_comb_0 <- temp_all_comb_1 <- unique_input[which_row, ]
                                       if (length(all_possible_rlz_0) > 0) {
                                         temp_all_comb_0 <- suppressWarnings(cbind(temp_all_comb_0, all_possible_rlz_0))
                                         temp_all_comb_1 <- suppressWarnings(cbind(temp_all_comb_1, all_possible_rlz_1))
                                       }
                                       if (length(intervention_variables_loc_needed) > 0) {
                                         for (i in 1:length(intervention_variables_loc_needed)) {
                                           temp_all_comb_0[, intervention_variables_loc[i]] <- intervention_levels_control[i]
                                           temp_all_comb_1[, intervention_variables_loc[i]] <- intervention_levels_treat[i]
                                         }
                                       }
                                       # for all non-A, non-0 variables, calculate the variable by rule
                                       # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
                                       loc_Z_need <- loc_Z[loc_Z > loc_node]  # only integrate previous variables; debug: not including current
                                       temp_list_0 <- lapply(loc_Z_need,
                                                             function(each_t) {
                                                               left_join(temp_all_comb_0, list_all_predicted_lkd[[each_t]])$output
                                                             })
                                       loc_other_needed <- loc_RLY[loc_RLY > loc_node]
                                       temp_list_1 <- lapply(loc_other_needed,
                                                             function(each_t) {
                                                               left_join(temp_all_comb_1, list_all_predicted_lkd[[each_t]])$output
                                                             })
                                       temp_list <- c(temp_list_0, temp_list_1)
                                       pmap_dbl(temp_list, prod) %>% sum %>% return
                                     })
      )
      Q0_current <- left_join(data_temp, library_output)$output
      if (loc_node == length(temp_node_names)) Q0_current <- rep(0, length(Q0_current))
    }

    # get H current
    {
      data_temp <- df_to_update[1:ind_var]  # take =1 probs
      temp_ind <- ind_var

      all_observed_1 <- all_observed_0 <- data_temp
      if (length(intervention_variables_loc_needed) > 0) {
        for (i in 1:length(intervention_variables_loc_needed)) {
          all_observed_0[, intervention_variables_loc[i]] <- intervention_levels_control[i]
          all_observed_1[, intervention_variables_loc[i]] <- intervention_levels_treat[i]
        }
      }

      # R and L (including Y before t = tau) (or not Z not A not t=0, not Y_tau) can be calculated in the same way
      if (loc_node %in% loc_RLY) {
        loc_A_needed <- loc_A[loc_A < loc_node]  # all needed A nodes
        loc_Z_needed <- loc_Z[loc_Z < loc_node]  # all needed Z nodes

        # all predicted probs at these locs are needed
        # needed inputs are in either all_observed_1 or all_observed_0

        # these A probs will be taken as product
        part_A <- lapply(loc_A_needed, function(k) {
          left_join(all_observed_1, list_all_predicted_lkd[[k]])$output
        }) %>% pmap_dbl(prod)

        # these ratios Z probs will be taken as product
        part_Z <- ifelse_vec(length(loc_Z_needed) == 0, rep(1, length(part_A)),
                             lapply(loc_Z_needed, function(k) {
                               left_join(all_observed_0, list_all_predicted_lkd[[k]])$output /
                                 left_join(all_observed_1, list_all_predicted_lkd[[k]])$output
                             }) %>% pmap_dbl(prod))

        H_current <- ifelse(
          # data_temp[last(intervention_variables_loc_needed)] == 1
          apply(data_temp[intervention_variables_loc_needed] == 1, 1, prod) == 1
          , 1/part_A*part_Z, 0) %>% as.vector
      }
      # Z nodes
      if (loc_node %in% loc_Z) {
        loc_A_needed <- loc_A[loc_A < loc_node]  # all needed A nodes
        loc_RLY_needed <- loc_RLY[loc_RLY < loc_node]

        # these A probs will be taken as product
        part_A <- lapply(loc_A_needed, function(k) {
          left_join(all_observed_0, list_all_predicted_lkd[[k]])$output
        }) %>% pmap_dbl(prod)

        # these ratios Z probs will be taken as product
        part_LR <- lapply(loc_RLY_needed, function(k) {
          left_join(all_observed_1, list_all_predicted_lkd[[k]])$output /
            left_join(all_observed_0, list_all_predicted_lkd[[k]])$output
        }) %>% pmap_dbl(prod)

        H_current <- ifelse(
          # ZW todo: stochastic intervention
          # data_temp[last(intervention_variables_loc_needed)] == 0
          apply(data_temp[intervention_variables_loc_needed] == 0, 1, prod) == 1
          , 1/part_A*part_LR, 0)
      }
    }

    # get newH current
    newH_current <- H_current * (Q1_current - Q0_current)
    return(newH_current)
  } else {
    return(NULL)
  }
}




#' @export
get_obs_Q <- function(tmle_task, obs_data, list_H,
                      intervention_variables, intervention_levels_treat, intervention_levels_control,
                      list_all_predicted_lkd,
                      lt) {
  obs_variable_names <- colnames(obs_data)
  temp_node_names <- names(tmle_task$npsem)
  loc_A <- grep("A", temp_node_names)
  loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
  loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
  intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
  list_Q <- list()
  for (loc_node in 1:length(list_H)) {
    if(!is.null(list_H[[loc_node]])) {
      # the Q integral at the previous variable; current inserted as lt; remove the lt setting if current Q is wanted
      # for the order of ARZL, at each RZL node, integrate out all children and set current = lt
      # note that no A probs will be involved in the product
      # debug: set all A probs as 1, because of Q definition

      # generate all possible 0, 1 valued rlz; set A to 0 first; set y_tau to always 1;
      loc_current_var <- which(obs_variable_names == tmle_task$npsem[[loc_node]]$variables)
      all_possible_RZLY_1 <- expand_values(obs_variable_names, to_drop = c(1:(loc_current_var-1) ),
                                           rule_variables = c(last(obs_variable_names), obs_variable_names[loc_current_var], intervention_variables),
                                           rule_values = c(1, lt, intervention_levels_treat))
      all_possible_RZLY_0 <- expand_values(obs_variable_names, to_drop = c(1:(loc_current_var-1) ),
                                           rule_variables = c(last(obs_variable_names), obs_variable_names[loc_current_var], intervention_variables),
                                           rule_values = c(1, lt, intervention_levels_control))

      # for each observed L_0 vector, generate all needed combinations, one version for A = 1, one version for A = 0
      unique_input <- obs_data[, 1:(loc_current_var-1)] %>% unique
      library_output <- data.frame(unique_input, output =
                                     map_dbl(1:nrow(unique_input), function(which_row) {
                                       # probs in the integrals, A=1 or A=0 is inserted
                                       temp_all_comb_0 <- data.frame(unique_input[which_row, ], all_possible_RZLY_0)
                                       temp_all_comb_1 <- data.frame(unique_input[which_row, ], all_possible_RZLY_1)
                                       for (i in 1:length(intervention_variables_loc)) {
                                         temp_all_comb_0[, intervention_variables_loc[i]] <- intervention_levels_control[i]
                                         temp_all_comb_1[, intervention_variables_loc[i]] <- intervention_levels_treat[i]
                                       }
                                       # for all non-A, non-0 variables, calculate the variable by rule
                                       # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
                                       loc_Z_needed <- loc_Z[loc_Z > loc_node]  # only product children variables
                                       temp_list_0 <- lapply(loc_Z_needed,
                                                             function(each_t) {
                                                               left_join(temp_all_comb_0, list_all_predicted_lkd[[each_t]])$output
                                                             })
                                       loc_RLY_needed <- loc_RLY[loc_RLY > loc_node]
                                       temp_list_1 <- lapply(loc_RLY_needed,
                                                             function(each_t) {
                                                               left_join(temp_all_comb_1, list_all_predicted_lkd[[each_t]])$output
                                                             })
                                       temp_list <- c(temp_list_0, temp_list_1)
                                       pmap_dbl(temp_list, prod) %>% sum %>% return
                                     })
      )
      list_Q[[loc_node]] <- left_join(obs_data[, 1:(loc_current_var-1)], library_output)$output
      if (loc_node == length(list_H))  list_Q[[loc_node]] <- rep(lt, nrow(obs_data))
    }
  }

  return(list_Q)
}





#' @export
get_obs_H <- function(tmle_task, obs_data, current_likelihood,
                      cf_task_treatment, cf_task_control,
                      intervention_variables, intervention_levels_treat, intervention_levels_control
) {
  obs_variable_names <- colnames(obs_data)
  temp_node_names <- names(tmle_task$npsem)
  loc_A <- grep("A", temp_node_names)
  loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
  loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
  intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))

  # get a list of corresponding H covariates; ordered by nodes, not variables
  list_H <- list()
  # calculate RLY nodes
  for (temp_ind in loc_RLY) {
    loc_A_needed <- loc_A[loc_A < temp_ind]  # all needed A nodes
    loc_Z_needed <- loc_Z[loc_Z < temp_ind]  # all needed Z nodes
    # this is the At indicators for H_RLY; now
    A_ind <-
      # obs_data[[tmle_task$npsem[[last(loc_A_needed)]]$variables]] == intervention_levels_treat[tmle_task$npsem[[last(loc_A_needed)]]$variables]
      apply(sapply(loc_A_needed, function(k) {
        obs_data[[tmle_task$npsem[[k]]$variables]] == intervention_levels_treat[tmle_task$npsem[[k]]$variables]
      }), 1, prod) == 1
    # A_ind <- obs_data[[temp_node_names[last(loc_A_needed)]]]  # using variable names (rather than node names) to inquire obs_data
    # these A probs will be taken as product
    part_A <- lapply(loc_A_needed, function(k) current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k])) %>% pmap_dbl(prod)  # this is the likelihood of being 1
    part_Z <- lapply(loc_Z_needed, function(k) {
      current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k]) /
        current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k])
    }) %>% pmap_dbl(prod)
    if(length(part_Z) == 0) part_Z <- 1

    list_H[[temp_ind]] <- ifelse(A_ind, 1/part_A*part_Z, 0) %>% as.vector
  }
  # calculate Z nodes
  for (temp_ind in loc_Z) {
    loc_A_needed <- loc_A[loc_A < temp_ind]  # all needed A nodes
    loc_RLY_needed <- loc_RLY[loc_RLY < temp_ind]
    A_ind <-
      # obs_data[[tmle_task$npsem[[last(loc_A_needed)]]$variables]] == intervention_levels_control[tmle_task$npsem[[last(loc_A_needed)]]$variables]
      apply(sapply(loc_A_needed, function(k) {
        obs_data[[tmle_task$npsem[[k]]$variables]] == intervention_levels_control[tmle_task$npsem[[k]]$variables]
      }), 1, prod) == 1
    part_A <- lapply(loc_A_needed, function(k) current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k])) %>% pmap_dbl(prod)
    part_RLY <- lapply(loc_RLY_needed, function(k) {
      current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k]) /
        current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k])
    }) %>% pmap_dbl(prod)
    list_H[[temp_ind]] <- ifelse(A_ind, 1/part_A*part_RLY, 0) %>% as.vector
  }
  return(list_H)
}





#' @export
get_obs_H_full <- function(tmle_task, obs_data, current_likelihood,
                           cf_task_treatment, cf_task_control,
                           intervention_variables, intervention_levels_treat, intervention_levels_control
) {
  obs_variable_names <- colnames(obs_data)
  temp_node_names <- names(tmle_task$npsem)
  loc_A <- grep("A", temp_node_names)
  loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
  loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
  intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))

  # get a list of corresponding H covariates; ordered by nodes, not variables
  list_H <- list()
  # calculate RLY nodes
  for (temp_ind in loc_RLY) {
    loc_A_needed <- loc_A
    loc_Z_needed <- loc_Z
    # this is the At indicators for H_RLY; now
    A_ind <-
      # obs_data[[tmle_task$npsem[[last(loc_A_needed)]]$variables]] == intervention_levels_treat[tmle_task$npsem[[last(loc_A_needed)]]$variables]
      apply(sapply(loc_A_needed, function(k) {
        obs_data[[tmle_task$npsem[[k]]$variables]] == intervention_levels_treat[tmle_task$npsem[[k]]$variables]
      }), 1, prod) == 1
    # A_ind <- obs_data[[temp_node_names[last(loc_A_needed)]]]  # using variable names (rather than node names) to inquire obs_data
    # these A probs will be taken as product
    part_A <- lapply(loc_A_needed, function(k) current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k])) %>% pmap_dbl(prod)  # this is the likelihood of being 1
    part_Z <- lapply(loc_Z_needed, function(k) {
      current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k]) /
        current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k])
    }) %>% pmap_dbl(prod)
    if(length(part_Z) == 0) part_Z <- 1

    list_H[[temp_ind]] <- ifelse(A_ind, 1/part_A*part_Z, 0) %>% as.vector
  }
  # calculate Z nodes
  for (temp_ind in loc_Z) {
    loc_A_needed <- loc_A
    loc_RLY_needed <- loc_RLY
    A_ind <-
      # obs_data[[tmle_task$npsem[[last(loc_A_needed)]]$variables]] == intervention_levels_control[tmle_task$npsem[[last(loc_A_needed)]]$variables]
      apply(sapply(loc_A_needed, function(k) {
        obs_data[[tmle_task$npsem[[k]]$variables]] == intervention_levels_control[tmle_task$npsem[[k]]$variables]
      }), 1, prod) == 1
    part_A <- lapply(loc_A_needed, function(k) current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k])) %>% pmap_dbl(prod)
    part_RLY <- lapply(loc_RLY_needed, function(k) {
      current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k]) /
        current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k])
    }) %>% pmap_dbl(prod)
    list_H[[temp_ind]] <- ifelse(A_ind, 1/part_A*part_RLY, 0) %>% as.vector
  }
  return(list_H)
}




#' @export
ifelse_vec <- function(condition, out1, out2) {
  if (condition) return(out1) else return(out2)
}

#' @export
logit <- function(x) {
  log(x/(1-x))
}

#' @export
expit <- function(x) {
  exp(x) / (1 + exp(x))
}

#' @export
scale_01 <- function(x) scale(x, center = min(x), scale = max(x) - min(x))





#' @export
ipw_middle <- function(task, lik, ipw_args, fold_number){

  cf_likelihood_control = ipw_args$cf_likelihood_control
  cf_likelihood_treatment = ipw_args$cf_likelihood_treatment
  intervention_list_treatment <- ipw_args$intervention_list_treatment
  intervention_list_control <- ipw_args$intervention_list_control
  # todo: extend for stochastic
  cf_task_treatment <- cf_likelihood_treatment$enumerate_cf_tasks(task)[[1]]
  cf_task_control <- cf_likelihood_control$enumerate_cf_tasks(task)[[1]]

  intervention_nodes <- union(names(intervention_list_treatment), names(intervention_list_control))

  temp_node_names <- names(task$npsem)
  loc_A <- grep("A", temp_node_names)
  loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
  loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
  if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)

  Y <- task$get_tmle_node(last(temp_node_names), format = T)[[1]]

  # get list of all possible predicted lkds
  obs_data <- task$data %>% as.data.frame %>% select(-c(id, t))
  obs_variable_names <- colnames(obs_data)
  # ZW todo: to handle long format and wide format
  # ZW todo: see if observed_likelihood needs to change to targeted likelihood

  intervention_variables <- map_chr(task$npsem[intervention_nodes], ~.x$variables)
  intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
  intervention_levels_treat <- map_dbl(intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
  intervention_levels_control <- map_dbl(intervention_list_control, ~.x$value %>% as.character %>% as.numeric)

  list_H <- get_obs_H_full(task, obs_data, current_likelihood = lik,
                           cf_task_treatment, cf_task_control,
                           intervention_variables, intervention_levels_treat, intervention_levels_control)


  list_newH <- list()
  for (ind_var in 1:length(list_H)) {
    if(!is.null(list_H[[ind_var]])) {
      if (ind_var %in% loc_Z) list_newH[[ind_var]] <- (list_H[[ind_var]] * Y) %>% as.matrix
      if (ind_var %in% loc_RLY) list_newH[[ind_var]] <- (list_H[[ind_var]] * Y) %>% as.matrix
    }
  }
  names(list_newH) <- temp_node_names

  return(list_newH)

  # ZW todo: for categorical variables
}


#' @export
gradient_generator_middle <- function(tmle_task, lik,  node, include_outcome = T, ipw_args = NULL, fold_number){

  task <- tmle_task$get_regression_task(node)
  IC <- ipw_middle(tmle_task, lik,  ipw_args, fold_number)[[node]] %>% as.vector
  new_data <- data.table(IC = IC)
  set(new_data, , node, task$Y)
  cols <- task$add_columns(new_data)
  task <- task$clone()
  nodes <- task$nodes
  nodes$outcome <- "IC"
  nodes$covariates <- c(nodes$covariates, node)

  task$initialize(
    task$internal_data,
    nodes = nodes,
    folds = task$folds,
    column_names = cols,
    row_index = task$row_index,
    outcome_type = "continuous"
  )
  return(task)
}
