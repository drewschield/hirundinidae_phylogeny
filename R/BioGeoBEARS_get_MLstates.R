#### See BioGeoBEARS_results_figure.R

# To see the plot_BioGeoBEARS_results() function:
plot_BioGeoBEARS_results

# Set the colors
palette <- c(254,224,139, # A-Afrotropical
             50,136,189, # B-Palearctic
             171,221,164, # C-Oriental
             253,174,97, # D-Australian
             213,62,79, # E-Oceanian
             102,194,165, # F-Nearctic
             230,245,152, # G-Panamanian
             244,109,67) # H-Neotropical
             

temp.matrix <- matrix(palette, nrow = 3)
rownames(temp.matrix) <- c("red", "green", "blue")
colors_matrix <- temp.matrix # colors_matrix must changed to colors_matrix = colors_matrix on line 160

# Define some stuff, otherwise the function won't work.
results_object <- resDIVALIKEj
addl_params=list("j")
juststats = FALSE
analysis_titletxt = NULL
dej_params_row = NULL
if_ties = "takefirst"
colors_list_for_states = NULL

tmp_fg = par("fg")
par(fg = "black")
BioGeoBEARS_run_object = results_object$inputs
if (is.null(tr)) {
  tr = check_trfn(trfn = BioGeoBEARS_run_object$trfn)
}
tr_pruningwise = reorder(tr, "pruningwise")
tips = 1:length(tr_pruningwise$tip.label)
nodes = (length(tr_pruningwise$tip.label) + 1):(length(tr_pruningwise$tip.label) + 
                                                  tr_pruningwise$Nnode)
if (is.null(tipranges)) {
  if (BioGeoBEARS_run_object$use_detection_model == FALSE) {
    tipranges = getranges_from_LagrangePHYLIP(lgdata_fn = np(BioGeoBEARS_run_object$geogfn))
  }
  if (BioGeoBEARS_run_object$use_detection_model == TRUE) {
    if (BioGeoBEARS_run_object$use_detection_model == 
        TRUE) {
      tipranges = tipranges_from_detects_fn(detects_fn = BioGeoBEARS_run_object$detects_fn)
    }
  }
}
areas = getareas_from_tipranges_object(tipranges)
areas
numareas = length(areas)
numareas
if (!is.na(results_object$inputs$max_range_size)) {
  max_range_size = results_object$inputs$max_range_size
} else {
  max_range_size = length(areas)
}
max_range_size
if (is.null(results_object$inputs$states_list)) {
  numstates = numstates_from_numareas(numareas = length(areas), 
                                      maxareas = max_range_size, include_null_range = results_object$inputs$include_null_range)
  numstates
  states_list_areaLetters = areas_list_to_states_list_new(areas, 
                                                          maxareas = max_range_size, include_null_range = results_object$inputs$include_null_range)
  states_list_0based_index = rcpp_areas_list_to_states_list(areas, 
                                                            maxareas = max_range_size, include_null_range = results_object$inputs$include_null_range)
} else {
  states_list_0based_index = results_object$inputs$states_list
}
param_ests = extract_params_from_BioGeoBEARS_results_object(results_object, 
                                                            returnwhat = "table", addl_params = addl_params, paramsstr_digits = 4)
if (juststats == TRUE) {
  return(param_ests)
} else {
  paramstr = extract_params_from_BioGeoBEARS_results_object(results_object, 
                                                            returnwhat = "string", addl_params = addl_params, 
                                                            paramsstr_digits = 4)
}
param_names = extract_params_from_BioGeoBEARS_results_object(results_object, 
                                                             returnwhat = "param_names", addl_params = addl_params, 
                                                             paramsstr_digits = 4)
if (is.null(analysis_titletxt)) {
  tmptxt = results_object$inputs$description
  if (any(is.null(tmptxt), tmptxt == "", tmptxt == "defaults", 
          tmptxt == "default")) {
    analysis_titletxt = ""
  } else {
    analysis_titletxt = results_object$inputs$description
  }
}
if (is.null(dej_params_row)) {
  analysis_titletxt = paste(analysis_titletxt, "\n", "ancstates: global optim, ", 
                            max_range_size, " areas max. ", paramstr, sep = "")
  analysis_titletxt
} else {
  dej_params_row
  brate_col_TF = names(dej_params_row) == "brate"
  brate_col = (1:length(dej_params_row))[brate_col_TF]
  biogeog_params = dej_params_row[1:(brate_col - 1)]
  biogeog_param_names = names(dej_params_row)[1:(brate_col - 
                                                   1)]
  equals_col = "="
  tmpcols = cbind(biogeog_param_names, equals_col, unlist(biogeog_params))
  tmpcols
  txtrows = apply(X = tmpcols, MARGIN = 1, FUN = paste, 
                  sep = "", collapse = "")
  txtrows
  biogeog_params_txt = paste(txtrows, sep = "", collapse = "; ")
  biogeog_params_txt
  titletxt2 = bquote(paste(.(max_range_size), " areas max., ", 
                           .(biogeog_params_txt), "; ", lambda, "=", .(dej_params_row$brate), 
                           "; ", mu, "=", .(dej_params_row$drate), "; ", alpha, 
                           "=", .(dej_params_row$rangesize_b_exponent), "; ", 
                           omega, "=", .(dej_params_row$rangesize_d_exponent), 
                           "", sep = ""))
}
leftright_nodes_matrix = get_leftright_nodes_matrix_from_results(tr_pruningwise)
marprobs = results_object$ML_marginal_prob_each_state_at_branch_bottom_below_node
left_ML_marginals_by_node = marprobs[leftright_nodes_matrix[, 
                                                            2], ]
right_ML_marginals_by_node = marprobs[leftright_nodes_matrix[, 
                                                             1], ]
right_ML_marginals_by_node
if (is.null(dim(left_ML_marginals_by_node))) {
  left_ML_marginals_by_node = matrix(data = left_ML_marginals_by_node, 
                                     nrow = 1)
}
if (is.null(dim(right_ML_marginals_by_node))) {
  right_ML_marginals_by_node = matrix(data = right_ML_marginals_by_node, 
                                      nrow = 1)
}
relprobs_matrix = results_object$ML_marginal_prob_each_state_at_branch_top_AT_node
if (length(nodes) > 1) {
  relprobs_matrix_for_internal_states = relprobs_matrix[nodes, 
                                                        ]
} else {
  relprobs_matrix_for_internal_states = relprobs_matrix[nodes, 
                                                        ]
  relprobs_matrix_for_internal_states = matrix(data = relprobs_matrix_for_internal_states, 
                                               nrow = 1, ncol = ncol(relprobs_matrix))
}
relprobs_matrix
if (is.null(states_list_0based_index)) {
  statenames = areas_list_to_states_list_new(areas, maxareas = max_range_size, 
                                             include_null_range = results_object$inputs$include_null_range, 
                                             split_ABC = FALSE)
  ranges_list = as.list(statenames)
  statenames
} else {
  ranges_list = states_list_0based_to_ranges_txt_list(state_indices_0based = states_list_0based_index, 
                                                      areanames = areas)
  ranges_list
  statenames = unlist(ranges_list)
  statenames
}
MLprobs = get_ML_probs(relprobs_matrix)
MLprobs
MLstates = get_ML_states_from_relprobs(relprobs_matrix, statenames, 
                                       returnwhat = "states", if_ties = if_ties)
if (is.null(colors_list_for_states)) {
  colors_matrix = colors_matrix
  colors_list_for_states = mix_colors_for_states(colors_matrix, 
                                                 states_list_0based_index, plot_null_range = results_object$inputs$include_null_range)
  colors_list_for_states
}
if (is.null(ranges_list)) {
  possible_ranges_list_txt = areas_list_to_states_list_new(areas, 
                                                           maxareas = max_range_size, split_ABC = FALSE, include_null_range = results_object$inputs$include_null_range)
} else {
  possible_ranges_list_txt = ranges_list
}
cols_byNode = rangestxt_to_colors(possible_ranges_list_txt, 
                                  colors_list_for_states, MLstates)
