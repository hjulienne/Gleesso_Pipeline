
#' @title Community attribution stability table
#' @description Compute the stability of
#' species community attribution from bootstraped graphs and community
#' table abondance for robust attribution.
#' @param graphs_folder : folder where all graphs are placed (with the / at the end please ^^)
#' @param alluvial_diagnostic : file name for the alluvial graph
#' @param N_alluvial : number of graph to represent on a graph
#' @param join_type : Should we work on the union of species or the intersection
#' @return a list object containing the following elements:
#' staby : table of species stability with community assignation for all graphs
#' stab_n_taxo : table of species stability and taxonomy
#' Robust_community_stability_[...] : abundance of communities for species with a stability above the specified treshold
#' Robust_community_stability_[...]_silhouette_[...] : idem but species also have a silhouette above a specified treshold
#' @export

Robust_table_community <- function(graphs_folder,
    alluvial_diagnostic_file, taxo,
     N_alluvial = 10,
     join_type = "outer",
     stability_treshold = 0.6,
     silhouette_treshold = 0.1
)
 {
  files_all = system(paste("ls ", graphs_folder), intern = TRUE)
  graph_names = files_all[grepl("_\\d+_community_tagged$", files_all)]
  graph_batch = list()
  MGS_abund = readRDS(paste0(graphs_folder, "abund_by_species.rds"))
  # Retrieve all graphs in the list :
  for(gn in graph_names)
  {
    graph_batch[[gn]] = readRDS(paste0(graphs_folder, gn))
  }
  Ngraphs = length(graph_names)
  name_ref =  paste0(graphs_folder, files_all[grepl("_all_samples_community_tagged$", files_all)])
  graph_ref = readRDS(name_ref)
  graphbatch_converted = batch_converter(graph_batch, graph_ref)
  # plot an alluvial plot

  id_repr_allu = sample(Ngraphs, N_alluvial)
  graph_labels = c("graph reference", paste("bootstrap", id_repr_allu))
  svg(alluvial_diagnostic_file)
  print("#########################")
  print("# Plotting alluvial plot")
  parrallel_coord_community(graphbatch_converted[c("graph_reference", names(graphbatch_converted)[id_repr_allu])], graph_labels, join_type = join_type)
  dev.off()

  staby = stability_index_converter(graphbatch_converted, join_type = join_type)
  stab_n_taxo = cbind(staby[!grepl("graph",names(staby))],
                taxo[row.names(staby), c("genus","family","order","class", "superkingdom", "superkingdom")])

  # Compute silhouette based on the walktrap distance between nodes
  graph_all_fn = files_all[grep("all_samples.gl.RDS$", files_all)]
  gphref_igraph = readRDS(paste0(graphs_folder, graph_all_fn))
  pos.grph = delete.edges(gphref_igraph, which(E(gphref_igraph)$weight < 0))
  print("walktrap distance computing")
  Dwkt = walktrap_distance(pos.grph, 9)
  Swkt = Silhouette_to_community(Dwkt, graphbatch_converted[["graph_reference"]])
  stab_n_taxo["Silhouette_walktrap"] = Swkt[row.names(stab_n_taxo),3]
  Nodes_all = graphbatch_converted[["graph_reference"]]

  # Select nodes that are attributed to the same community more than stab treshold
  Stab_Nodes = Nodes_all[row.names(stab_n_taxo)[which(stab_n_taxo$Stability_score >= stability_treshold)],]
  community_table_stab = Compute_community_abondance(Stab_Nodes, MGS_abund, taxo)
  name_stab_com = paste0("Robust_community_stability_", stability_treshold)


  Stab_n_fat_Nodes = Nodes_all[row.names(stab_n_taxo)[which((stab_n_taxo$Stability_score >= stability_treshold) & (stab_n_taxo$Silhouette_walktrap > silhouette_treshold))],]
  community_table_stab_n_fat = Compute_community_abondance(Stab_n_fat_Nodes, MGS_abund, taxo)
  name_stab_n_fat_com = paste0("Robust_community_stability_", stability_treshold, "_silhouette_", silhouette_treshold)

  result_bundle = list(stability = staby)
  result_bundle[["stab_n_taxo"]] = stab_n_taxo
  result_bundle[[name_stab_com]] = community_table_stab
  result_bundle[[name_stab_n_fat_com]] = community_table_stab_n_fat

  return(result_bundle)
}


#' @title create_graph_robust_community_tags
#' @description
#' create a gephi file for nodes with a robust community attribution
#' @param fout : where to save the csv tables with nodes with the robust community column
#' @param abund_by_species : abundance mean group by species or
#' @param taxo_by_species : taxo grouped at the species level
#' @param Robust_table_community : Community attribution
#' @param Nodes tables computed on all samples
#'

create_graph_robust_community_tags <- function(model_folder,
  fout,
  abund_by_species,
  taxo_by_species,
  Robust_table_community,
  Nodes_table_on_all_samples,
  variability_treshold = NULL
 )
{
    fin_glasso <- model_folder + "/glasso_model_all_samples.gl"
    Nodes_table_on_all_samples["Robust_community"] = ""

    list_of_species = lapply(commu_table[1,], strsplit, split = "-")
    # Attribute robust community to the right species
    # in the node table

    for(com in names(list_of_species))
    {
      Nodes_table_on_all_samples[list_of_species[[com]][[1]], "Robust_community"] = com
    }

    fout_gexf = paste0(fout, ".gexf")
    fout_csv = paste0(fout, ".csv")

    Nodes <- create_graph(
                fin_glasso,
                fout_gexf,
                abund_by_species,
                taxo_by_species,
                community=FALSE,
                variability_treshold = variability_treshold,
                additional_info = Nodes_table_on_all_samples$Robust_community)

       write.csv(Nodes, file = fout_csv)
}

#' @title concordance_table
#' @description
#' Enable one to assess if community found in diverse cohort are the same
#' Used to generate the alluvial plot
#' @param join_type : how to join graph row (outer joins or inner join). Outer join means that all species present in at least one graph will be taken into account. Inner join means that only species present in all graphs will treated.
#' @param nlist : list of graphs nodes table with the walktrap_community information Available
#' @param Graph_tags : list of graph labels
#' @export

concordance_table <- function(nlist, Graph_tags, join_type= "outer")
{
  # join species to create row index
    union_speciest = rownames(nlist[[1]])
    if(join_type=="inner")
    {
        for(n in nlist[-1])
        {
            union_speciest = intersect(union_speciest, rownames(n))
        }
    }
    else
    {
        if(join_type=="outer")
        {
            for(n in nlist[-1]){
            union_speciest = union(union_speciest, rownames(n))}
        }
        else{stop("the ways to merge you entered is not valid, you moron. Available option are: 'inner' or 'outer'")}
    }

    community_cohort = data.frame(row.names=union_speciest)
    # Retrieve community information on all graphs
    for(i in 1:length(nlist))
    {
      community_cohort[Graph_tags[i]] = as.character(nlist[[i]][row.names(community_cohort), "walktrap_community"])
    }

    # Add mean abundance on all bootstraps
    species_abondance = data.frame(row.names = union_speciest)
    for(i in 1:length(nlist))
    {
        species_abondance[i] = nlist[[i]][union_speciest, "abondance"]
    }

    v_abondance  = rowMeans(species_abondance, na.rm=TRUE)
    community_cohort['mean_abundance'] = v_abondance
    species_community_cohort = data.table(community_cohort)
    species_community_cohort_aggregated = species_community_cohort[,.(sum_ab= sum(mean_abundance), count = .N) , by= Graph_tags]

    res_comparison = list(species_community_cohort = species_community_cohort,  species_community_cohort_aggregated = species_community_cohort_aggregated)

    return(res_comparison)
}

#' @title parrallel_coord_community
#' @description Alluvial plot of concordance of community
#' Draw a parrallel coord graph of community belonging for
#' different graph object. Enable one to compare and understand the stability or discrepancy
#' between graph community
#' @param graph_node_list: a sequence of nodes tables with community annotated
#' @param Graph_tags: a sequence of str which are the name of nodes tables of graph_node_list
#' @param measure: the weight attributed to each CAG either "sum" of abundance or "count" of objects
#' @param color_graph : index of the graph used to color the parallel coordiante plot
#' @param join_type : take the intersection or the union of species?
#' @export

parrallel_coord_community <- function(graph_node_list,
   Graph_tags,
     measure= "sum_ab",
      color_graph = 1,
       join_type= "inner" )
{
    N_graph <- length(graph_node_list)
    pcd_agg <- concordance_table(graph_node_list, Graph_tags, join_type= join_type)[['species_community_cohort_aggregated']]

    pcd_agg = as.data.frame(pcd_agg, stringsAsFactors = FALSE)
    selec_not_all_NA <- !(apply(is.na(pcd_agg[, 1:N_graph]), 1, all))

    pcd_agg = pcd_agg[selec_not_all_NA,]

    pcd_agg[is.na(pcd_agg)] = ''

    large_palette <- c(brewer.pal(9,'Set1'), brewer.pal(7,'Set2'), brewer.pal(12,'Set3'))

    colors= sample(large_palette)[as.numeric(as.factor(pcd_agg[, Graph_tags[color_graph]]))]

    alluvial(pcd_agg[,Graph_tags], freq=pcd_agg[, measure], border='lightgrey', col= colors, cex=0.5, gap.width=0.1, cw=0.18, alpha=0.5)
                                        #hide = pcd_agg$count > 100,

}




#' @title walktrap_distance
#' @description reproduce the distance used in the walktrap community detection algorithm
#' @param pos.graph: an igraph object that contain only positive edges
#' @param n_steps: number of steps of the random walk on the graph


walktrap_distance <- function(pos.grph, n_steps)
{
# #
    mat_pos <- as_adjacency_matrix(pos.grph, sparse = FALSE, attr= "weight")

    inv_weighted_degre = 1.0/apply(mat_pos, 2, sum)
    inv_weighted_degre[is.infinite(inv_weighted_degre) | is.na(inv_weighted_degre)] = 0
    inv_degree_mat = diag( inv_weighted_degre)

    trans_mat = inv_degree_mat %*% mat_pos
    trans_7_step = trans_mat%^%n_steps
    Nbac = dim(mat_pos)[1]
    Dist_wt <- matrix(0, ncol=Nbac ,nrow=Nbac)

    dimnames(Dist_wt) <- dimnames(mat_pos)

    for( i in 1:dim(mat_pos)[1])
    {
        for(j in 1:dim(mat_pos)[2])
        {
            Dist_wt[i, j] = norm((sqrt(inv_degree_mat)%*%trans_7_step[i,])-(sqrt(inv_degree_mat)%*%trans_7_step[j,]))
        }
    }
    return(Dist_wt)
}



#' @title Silhouette_to_community
#' @description
#' Compute the silhouette cluster metric for all species
#' @param my_dist : distance matrix of species to all species
#' @param Nodes_with_com : Nodes table with a walktrap_community attribution


Silhouette_to_community <- function( my_dist, Nodes_with_com)
{
    cl_community = as.numeric(as.factor(as.character(Nodes_with_com$walktrap_community)))

    SDist = silhouette(cl_community, dmatrix = my_dist)
    S_distdt = as.data.frame(SDist[,1:3], row.names=row.names(Nodes_with_com))

    return(S_distdt)
}

#' @title community_converter
#' @description
#' Community converter function take to two nodes table as argument
#' and give a translation of each community to the other graph based on jacquart distance
#' If I may, it is automated community translation
#' @param nodes_graph1 is the table of nodes with the walktrap column properly filled
#' @param nodes_graph2 is the table of nodes of the second graph

community_converter <- function(nodes_graph1,
                                nodes_graph2,
                                join_type = "inner")
{

    if(join_type == "inner")
    {
        spc = intersect(row.names(nodes_graph1), row.names(nodes_graph2))
    }
    else{
        if(join_type == "outer")
        {
            spc = union(row.names(nodes_graph1), row.names(nodes_graph2))
        }
        else{
            print("your join_type argument doesn't make sense")
        }
    }

    nodes_graph1 = nodes_graph1[spc,]
    nodes_graph2 = nodes_graph2[spc,]

    nodes_graph1_com = na.omit(unique(as.character(nodes_graph1$walktrap_community, na.omit)))
    nodes_graph2_com = na.omit(unique(as.character(nodes_graph2$walktrap_community, na.omit)))

    vec1 = as.character(nodes_graph1$walktrap_community)
    vec1[is.na(vec1)] =""
    vec2 = as.character(nodes_graph2$walktrap_community)
    vec2[is.na(vec2)] =""
    spc = intersect(row.names(nodes_graph1), row.names(nodes_graph2))

    M = matrix(0 ,nrow=length(nodes_graph1_com), ncol=length(nodes_graph2_com))# ,dimnames = )
    dimnames(M) = list(nodes_graph1_com, nodes_graph2_com)
    for(com in  nodes_graph1_com)
    {
        indic = as.numeric(vec1==com)
        for(com2 in nodes_graph2_com)
        {
            indic2 = as.numeric(vec2== com2)
            M[com, com2] = dist.binary(rbind(indic, indic2), method=1)
        }
    }

    converter2to1= as.data.frame(t(rbind(nodes_graph2_com, nodes_graph1_com[apply(M, 2, which.min)])), stringsAsFactors= FALSE)
    converter2to1['numeric_codage'] = as.numeric(factor(converter2to1$V2, levels = levels(as.factor(nodes_graph1$walktrap_community))))

    row.names(converter2to1) = converter2to1$nodes_graph2_com
    converter2to1$nodes_graph2_com = NULL
    convertedtonum1 = as.numeric(as.factor(nodes_graph1$walktrap_community))
    convertedtonum2 = converter2to1[as.character(nodes_graph2$walktrap_community), 'numeric_codage']
    return(list(species_intersection= spc,
                dist_mat = M,
                converter2to1 =converter2to1,
                convertedtonum1 = as.numeric(as.factor(nodes_graph1$walktrap_community)) ,
                convertedtonum2 = convertedtonum2))
}


#' @title batch community converter
#' @description
#' convert community (give the same indice) of a batch
#' of graphs to the community of a reference graph
#' @param graph_batch : a list of graph to convert
#' @param graph_ref : the reference used for conversion

batch_converter <- function(graph_batch, graph_ref)
{
    N_graph <- length(graph_batch)

    for(cgraph in  1:N_graph)
    {
        converting_info = community_converter(graph_ref, graph_batch[[cgraph]])
        spc_to_convert =  which(rownames(graph_batch[[cgraph]]) %in% converting_info$species_intersection)
        graph_batch[[cgraph]]$walktrap_community_bkup = graph_batch[[cgraph]]$walktrap_community
        graph_batch[[cgraph]]$walktrap_community = 0
        graph_batch[[cgraph]]$walktrap_community[spc_to_convert] = converting_info$convertedtonum2
        graph_batch[[cgraph]]$walktrap_community[is.na(graph_batch[[cgraph]]$walktrap_community)] = 0
    }

    graph_ref$walktrap_community_bkup = graph_ref$walktrap_community
    graph_ref$walktrap_community = as.numeric(as.factor(graph_ref$walktrap_community))
    graph_ref$walktrap_community[is.na(graph_ref$walktrap_community)] = 0

    # Putting graph reference at the list of graph
    graph_batch[['graph_reference']] = graph_ref
    return(graph_batch)
}



#' @title stability_index_converter
#' @description
#' stability_index function look at the walktrap community of each species in a list of graph
#' then compute the number of graph where the species as been attributed to the same community
#' @param graph_batch : a list of graph with walktrap community converted to a ref


stability_index_converter <- function(graph_list, join_type = "outer")
{
    union_speciest = rownames(graph_list[[1]])
    N_graph = length(graph_list)

    for(grp in 2:N_graph)
    {
        union_speciest = union(union_speciest, rownames(graph_list[[grp]]))
    }

    community_accross_graphs = data.frame(row.names=union_speciest)
    stability_index= data.frame(row.names=union_speciest)
    for(grp in 1:N_graph)
    {
        name_col = "graph_"+ as.character(grp)
        community_accross_graphs[name_col] = -1
        community_accross_graphs[row.names(graph_list[[grp]]), name_col] = graph_list[[grp]]$walktrap_community
    }
    print(dim(stability_index))
    print(dim(community_accross_graphs))

    stability_index['Stability_score'] = -1
    stability_index['Stability_score'] = unlist(lapply(apply(community_accross_graphs, 1, table), max)) / N_graph

    stability_index['Stability_index'] = -1
    stability_index['Stability_index'] = unlist(lapply(apply(community_accross_graphs, 1, table), max))


    name.which.max <- function(x)
    {
        id= which.max(x)
        return(names(x)[id])
    }
    name.second.max <- function(x)
    {
        x = x[order(x)]
        return(names(x)[2])
    }

    stability_index['most_stab_attribution'] = -1
    stability_index['most_stab_attribution'] = unlist(lapply(apply(community_accross_graphs, 1, table), name.which.max))
    stability_index['second_most_stab_attribution'] = -1
    stability_index['second_most_stab_attribution'] = unlist(lapply(apply(community_accross_graphs, 1, table), name.second.max))

    return(cbind(stability_index, community_accross_graphs))
}
