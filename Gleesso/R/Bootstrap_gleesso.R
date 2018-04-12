#' @title bootstrap Gleesso_pipeline
#' @description
#' Apply the bootstrap pipeline to a fraction of the cohort.
#' A factor vector can be supplied to stratify the different samples
#' @param N_bootstrap : number of different bootstrap samples that should be drawn
#' @param fraction : fraction of the initial dataset that should be drawn to form each bootstrap samples
#' @param stratifying_vector : a factor vector that represent a class that should be evenly scattered between bootstrap samples
#' @param Graphs_folder : folder to output all graphs and all bootstrap samples (to keep track of which individual was used in each iteration)
#' @param tag_model :
#' @param ... : parameters to pass to the
#' @export

Gleesso_bootstrap <- function(N_bootstrap,
                              fraction,
                              tag_model,
                              variability_treshold,
                              community_table_folder,
                              model_folder, graph_folder,
                              MGS_file, taxo_file,
                            stratifying_vector = NULL, ...
                          )
                          {
  N_samples = dim(readRDS(MGS_file))[2]
  # test if dimension of inputs are compatible
  if(!is.null(stratifying_vector))
  {
     if(length(stratifying_vector) != N_samples)
    {
      stop("the vector of class to stratify by doesn't correspond to the dimension of the abundance matrix")
     }

     # Convert stratifying vector in character to avoid that the table
     # function create empty cases for non observed factor levels
     nsamp_tab = round(table(as.character(stratifying_vector))*fraction)
     which.factor <- function(x){ return(which(stratifying_vector==x))}
     label_classes <- sort(unique(stratifying_vector))
     list_ids  = lapply(label_classes, which.factor)
     names(list_ids) = label_classes
  }

  # iterate bootstraps
  for(i in 1:N_bootstrap)
  {
    if(is.null(variability_treshold))
    {
      tag_graph = paste(c("_", 0.05), collapse = "")

    }
     else{
       tag_graph = paste(c("_", variability_treshold), collapse = "")
     }
    tag_model_m = paste(c(tag_model,"_", i ), collapse = "")

    # compute sample to work on :
    if(!is.null(stratifying_vector))
    {
      boot_samples = c()
      for(x in names(list_ids))
      {
        boot_samples = c(boot_samples, sample(list_ids[[x]], nsamp_tab[x]))
      }
    }
    else
    {
      boot_samples = sample(N_samples, round(N_samples*fraction))
    }

    boot_vec = rep(FALSE, N_samples)
    boot_vec[boot_samples] = TRUE
    print(tag_graph)
    Gleesso_pipeline(community_table_folder,  MGS_file, taxo_file,
       model_folder, graph_folder, boot_vec, tag_model= tag_model_m,
        tag_graph=tag_graph,
        variability_treshold=variability_treshold, analysis_step=0, ... )
      }
      # After generating the bootstraps we compute the graph on all samples
      # The graph on all samples (with the maximum statistical )
      tag_model_m = paste(c(tag_model,"_all_samples"), collapse = "")
      tag_graph = paste(c("_", variability_treshold), collapse = "")
      boot_vec = rep(TRUE, N_samples)
      Gleesso_pipeline(community_table_folder,  MGS_file, taxo_file,
         model_folder, graph_folder, boot_vec, tag_model= tag_model_m,
          tag_graph=tag_graph,
          variability_treshold=variability_treshold, analysis_step=0, ... )
}
