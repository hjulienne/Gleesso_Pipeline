

#' @title Pipeline launcher
#' @description
#'Launching the complete glasso analysis (data prep,  graph inference, community detection, community abundances computation)
#' Warning : folder shouldn't be indicated with an / at the end 
#' @param data_folder : where the community abundance table will be written
#' @param MGS_file: The Metagenomic species abundance file in the RDS format
#' @param model_folder: where to save the graphical Lasso model object (output of the spiec.easi function)
#' @param graph_folder: where the graphical representation of the model will be saved (in the gephi format)
#' @param contrast_vector : a boolean vector to select a subset of the cohort (the model will be infered on samples with a TRUE value)
#' @param tag_model : a tag that will be inserted in output file to recognize the model parameter
#' @param tag_graph : a tag that will be inserted in output file to recognize the graph
#' @param community: Should the community structure be calculated?
#' @param nlambda: Number of regularisation parameter that will be tested (see huge::huge() documentation)
#' @param lambda.min.ratio: the smallest value of lambda as a fraction of its maximum (see huge::huge() documentation)
#' @param occurence_treshold: minimum fraction of samples where a species must be present to be taken into account in the analysis
#' @param abundance_threshold: minimum mean abundances for a species to be included in the analysis
#' @param variability_treshold: The maximum mean variability for graph edge presence. If null, the optimal covariance matrix will correspond to a variability of 0.05
#' @param analysis_step: At which step the analysis should be started (0: from scratch, 1: model inferences, 2: save gephi network, 3: Community detection). If NULL (default), the step will be infered from the files present in the output folders. Use analysis_step=0 to force computation from scratch. 
#' @param species_mode : should the graph inference be done on MGS (FALSE) or with MGS of the same specied merged togethere (TRUE)
#' @export
Gleesso_pipeline <- function(data_folder,
                             MGS_file, taxo_file,
                             model_folder,
                             graph_folder,
                             contrast_vector,
                             tag_model, tag_graph,
                             community = TRUE ,
                             nlambda = 20,
                             lambda.min.ratio = 0.1,
                             occurence_treshold = 0.05,
                             abundance_treshold = 10^-7,
                             variability_treshold=NULL,
                             analysis_step = NULL,
                             species_mode = TRUE
                             )
{
    # determine step of the analysis:
    if(is.null(analysis_step)){
        if(!file.exists(data_folder+'/'+"taxo_by_species.rds"))
        {
            analysis_step = 0
        }
        else
        {
            if(!file.exists(model_folder+'/glasso_model_'+ tag_model+'.gl'))
            {
                analysis_step = 1
            }
            else
            {
                if(!file.exists(graph_folder+'/graph_'+tag_model+ tag_graph+'.gexf'))
                {
                    analysis_step = 2
                }
                else
                {
                    if(!file.exists(data_folder+'/community_abundance_'+tag_model+ tag_graph+'.gexf'))
                    {
                        analysis_step = 3
                    }
                }
            }
        }
    }
    ##Compute steps if undone#
    ## 
    print("Starting Gleesso analysis for:")
    print("#################   "+tag_model+"_"+ tag_graph+"   #################")
    print("At analysis step :"+ as.character(analysis_step))
    ## ##################################

    if(analysis_step == 0){  #data_prep
        abund_mgs <- readRDS(MGS_file)
        taxo_mgs <- readRDS(taxo_file)

        if(species_mode == TRUE)
        {
        taxo_by_species <- get_taxo_by_species(taxo_mgs)
        abund_by_species <-  get_abundance_by_taxo_species(abund_mgs, taxo_mgs)
       }
    else
    {
        taxo_by_species <- generate_annotation_CAG_level(taxo_mgs)
        abund_by_species <-  abund_mgs
        
    }
       fout = data_folder+ "/taxo_by_species.rds"
        saveRDS(object = taxo_by_species, file = fout)
        fout = data_folder+ "/abund_by_species.rds"
        saveRDS(object = abund_by_species, file = fout)
  
    }
    else
    {
        fout = data_folder+ "/taxo_by_species.rds"
        taxo_by_species <- readRDS(fout)
        fout = data_folder+ "/abund_by_species.rds"
        abund_by_species <- readRDS(fout)
    }
    
    
    # Compute Glasso model
    if(analysis_step <= 1)
    {
        fout = model_folder + "/glasso_model_" + tag_model+ ".gl"
        Compute_graph(abund_by_species, contrast_vector, fout, nlambda=nlambda, lambda.min.ratio=lambda.min.ratio, abundance_treshold = abundance_treshold, occurence_treshold = occurence_treshold)       
    }

    # save and write network
    if(analysis_step <= 3)
    {
        fin_glasso <- model_folder + "/glasso_model_" + tag_model+ ".gl"
        fout_gexf <- graph_folder + "/graph_" + tag_model+ tag_graph
        Nodes <- create_graph(fin_glasso, fout_gexf, abund_by_species, taxo_by_species, community=community, variability_treshold = variability_treshold)
    }
    # save and write network with the community structure
    if(analysis_step <= 3 & community==TRUE)
    {
        fout = data_folder + "/community_abondance_" +tag_model+ tag_graph+ ".rds"
        community_abundance <- Compute_community_abondance(Nodes, abund_by_species, taxo_by_species)
        saveRDS( community_abundance, file= fout)
#### Add community tagging ####
        
        Nodes_tg = Tag_community_names(Com_tab = community_abundance, Nodes)
        fout = graph_folder + "/graph_" + tag_model + tag_graph+ "_community_tagged"
        saveRDS(Nodes_tg, file = fout )
        generate_graph_from_tables(fout, Nodes_tg, fin_glasso,  variability_treshold = variability_treshold)
       
    }
    print("########Ending Gleesso analysis##########")
}





#' @title Compute Glasso
#' @description
#' Compute the Glasso model from the MGS abundances on individuals with a true value in the
#' contrast vector.
#' @param MGS_abundance: The Metagenomic species abundance table
#' @param fout: where to save the model object
#' @param contrast_vector : a boolean vector to select a subset of the cohort (the model will be infered on samples with a TRUE value)
#' @param community: Should the community structure be calculated?
#' @param nlambda: Number of regularisation parameter that will be tested (see huge::huge() documentation)
#' @param lambda.min.ratio: the smallest value of lambda as a fraction of its maximum (see huge::huge() documentation)
#' @param occurence_treshold: minimum fraction of samples where a species must be present to be taken into account in the analysis
#' @param abundance_threshold: minimum mean abundances for a species to be included in the analysis
#' @param rep.num: Number of subsampling to compute the edge stability with "StarS" (see huge::huge.select documentation)
#' @param lambda: A sequence of regularisation parameter. If not null, it will override the automatic computation of the lambda sequence (with nlambda and lambda.min.ratio)
#' @export

Compute_graph <- function(MGS_abundance,
                          contrast_vector,
                          fout,
                          abundance_treshold=10^-7,
                          occurence_treshold = 0.05,
                          nlambda=20,
                          lambda.min.ratio=0.1,
                          lambda = NULL,
                          rep.num =20)
{

    mgs_tp = as.data.frame(MGS_abundance[, which(contrast_vector)])
  
    print("Number of patient in the contrast vector:")
    nb_pat <- dim(mgs_tp)[2]
    print(nb_pat)
    print("mgs_tp")
   
    #Species with sufficient abondance and prevalence to be taken into account
    
    detected_species <- rownames(mgs_tp)[(rowSums(mgs_tp) > (abundance_treshold*nb_pat)) & (apply((mgs_tp) > (abundance_treshold),1,mean ) > occurence_treshold)]
    mgs_tp <- mgs_tp[detected_species,]

    print("Number of species:")
    print(dim(mgs_tp)[1])
    
    if(is.null(lambda))
    {
        se.gl.tp = spiec.easi(round(t(mgs_tp)*10^6), method='glasso', icov.select=TRUE, sel.criterion='stars', icov.select.params=list(rep.num=rep.num), nlambda=nlambda, lambda.min.ratio=0.1)
    }
    else
    {
        se.gl.tp = spiec.easi(round(t(mgs_tp)*10^6), method='glasso', icov.select=TRUE, sel.criterion='stars', icov.select.params=list(rep.num=rep.num), lambda=seq_lambda)#nlambda=20, lambda.min.ratio=0.001)    
    }
    print(se.gl.tp$opt.lambda)
    print(se.gl.tp$opt.sparsity)
    plot(se.gl.tp)
    huge::plot.huge(se.gl.tp)
    save(se.gl.tp, file=fout)
    return(0)
}






#' @title Create a gephi format graph from the graphical Lasso  model
#' @description
#' The Function also compute the community structure of the graph with various algorithms (betweeness community, walktrap community...)
#' community specified by the user
#' @param file_input: emplacement of the GLASSO model object
#' @param file_output: where to save the network representation file
#' @param MGS_by_taxo_species: The Metagenomic species abundance table
#' @param species_taxo: The Metagenomic species taxonomy table
#' @param nspins: number of spin for the spinglass community detection algorithm
#' @param community: Should the community structure be calculated?
#' @param additional_info: a vector or data.frame containing information to add to the nodes table of the network
#' @param spinglass_opt : Should the number of spin of the spinglass community be optimized?
#' @param variability_treshold: The maximum mean variability for graph edge presence. If null, the optimal covariance matrix will correspond to a variability of 0.05
#' @export
create_graph <- function(
                         file_input,
                         file_output,
                         MGS_by_taxo_species,
                         species_taxo,
                         nspins = 20,
                         mod_rep = 10,
                         community=FALSE,
                         additional_info=NULL,
                         spinglass_opt= FALSE,
                         variability_treshold = NULL)
{
    se.gl = get(load(file_input))
    #Compute percentage where the MGS is present
    occurence <- round(apply(MGS_by_taxo_species > 10^-7, 1, mean)*100)
    # Compose abondance of the MGS only in individuals where it is present 
    prevalence <- round(apply(MGS_by_taxo_species , 1, mean_when_present)*100*10^6)
    # remove diagonal terms
    if(is.null(variability_treshold))
    {
        icov_mat = se.gl$opt.icov
    }
    else
    {
        id_lbd <- max(which(se.gl$variability < variability_treshold))
        var_graph <- se.gl.tp$variability[id_lbd]
        lbd <- se.gl.tp$lambda[id_lbd]
        icov_mat <- se.gl.tp$icov[[id_lbd]] 
        print("|||||||||||||| Graph parameters |||||||||||||||" )
        print("Lambda: " + as.character(lbd))
        print("Variability: " + as.character(var_graph))
        print('|||||||||||||||||||||||||||||||||||||||||||||||')
    }
    
    icov_mat = as.matrix(icov_mat)
    diag(icov_mat) = 0
    nodes_cag = colnames(se.gl$data)# useful only in CAG mode (not in species aggregation mode)
    
    cnames <- sapply(species_taxo[colnames(se.gl$data), "annot"], remove_xml_char)

        dimnames(icov_mat) = list(cnames, cnames)
 
    tnames <- sapply(rownames(species_taxo), remove_xml_char)
    rownames(species_taxo) <- tnames
    
    ig.gl <- igraph::graph.adjacency(as.matrix(-icov_mat), weighted=TRUE, mode='min')
    igraph::E(ig.gl)[weight < 0]$color <- adjustcolor('coral', alpha.f=0.6) 
    igraph::E(ig.gl)[weight > 0]$color <- adjustcolor('darkslateblue', alpha.f=0.6) 

    #"Out" rule on a bacteria
    Nodes <- data.frame(ID = c(1:igraph::vcount(ig.gl)), NAME = igraph::V(ig.gl)$name)
    edges <- as.data.frame(igraph::get.edges(ig.gl, 1:igraph::ecount(ig.gl)))
    
    if(community==TRUE)
    {
        con.comp <- igraph::clusters(ig.gl)
        con.grph <- igraph::induced_subgraph(ig.gl, con.comp$membership==con.comp$membership[which.max(con.comp$csize)])

        pos.grph <- igraph::delete.edges(ig.gl, which(igraph::E(ig.gl)$weight < 0 ))

        print("###### Computing other community algorithm ######")
        
        btw_com <- igraph::edge.betweenness.community(pos.grph,weights = igraph::E(pos.grph)$weight)
        walktrap.comm <- opt_walktrap(pos.grph) #walktrap.community(pos.grph, weight = E(pos.grph)$weight)

        
        print("##### Modularity Walktrap:")
        print(igraph::modularity(walktrap.comm), weights =  igraph::E(pos.grph)$weight)

        print("Modularity edge betweness:")
        print(igraph::modularity(btw_com), weights =  igraph::E(pos.grph)$weight)

        igraph::V(ig.gl)[walktrap.comm$names]$walk_com <- walktrap.comm$membership
        walktrap.comm <- igraph::walktrap.community(pos.grph, weight = -igraph::E(pos.grph)$weight)
        
        nodes_viz_att <- data.frame(
            occurence = occurence[colnames(se.gl$data)],
            abondance = prevalence[colnames(se.gl$data)],
            species =  sapply(species_taxo[nodes_cag,]$species, remove_xml_char), 
            genus = species_taxo[nodes_cag,]$genus,
            family = species_taxo[nodes_cag,]$family,
            phylum = species_taxo[nodes_cag,]$phylum,
            class = species_taxo[nodes_cag,]$class,
            order = species_taxo[nodes_cag,]$order,
                                        # connected_component = con.comp$membership,
                                        # spinglass_community = V(ig.gl)$spin_com,
                                        # betweness_community = V(ig.gl)$bet_com,
            walktrap_community = igraph::V(ig.gl)$walk_com  
        )
         # get optimized spinglass community
        if(spinglass_opt == TRUE)
        {
            spin_com <- opt_spinglass_com(con.grph, mod_rep, nspins)
            V(ig.gl)[spin_com$names]$spin_com <- spin_com$membership
            
            nodes_viz_att['spinglass_community'] = igraph::V(ig.gl)$spin_com
        }
    }
    else
    {
        nodes_viz_att <- data.frame(
            occurence = occurence[colnames(se.gl$data)],
            abondance = prevalence[colnames(se.gl$data)],
            species =  sapply(species_taxo[nodes_cag,]$species, remove_xml_char), 
            genus =  species_taxo[nodes_cag,]$genus,
            phylum = species_taxo[nodes_cag,]$phylum,
            class = species_taxo[nodes_cag,]$class,
            order = species_taxo[nodes_cag,]$order
        )

    }
    # Add example of profile
    for( c in sample(colnames(MGS_by_taxo_species),10))
    {
        nodes_viz_att[c]= round(MGS_by_taxo_species[nodes_cag,c]*10^8)
    }

    if(!is.null(additional_info))
    {
        nodes_viz_att <- cbind(nodes_viz_att, additional_info)
    }
    #  E(ig.gl)[weight < 0]$color <- adjustcolor('coral', alpha.f=0.6) 
    #  E(ig.gl)[weight > 0]$color <- adjustcolor('darkslateblue', alpha.f=0.6)
    
    edges_viz_att <- data.frame(WGH = abs(igraph::E(ig.gl)$weight), COLOR = igraph::E(ig.gl)$color)
    row.names(nodes_viz_att) <- remove_xml_char(row.names(nodes_viz_att))
    Nodes$NAME <- remove_xml_char(Nodes$NAME)
    print('allmost here')
    write.csv(nodes_viz_att, file=paste(file_output,"_nodes.csv", sep=''))
    write.csv(edges_viz_att, file=paste(file_output,"_edges.csv", sep=''))
    saveRDS(ig.gl, file = paste(file_output, '.gl.RDS', sep=""))
    
    rgexf::write.gexf(output= paste(file_output, '.gexf', sep=""), nodes = Nodes, edges = edges, edgesWeight =round(abs(E(ig.gl)$weight)*100+1), nodesAtt = nodes_viz_att, edgesAtt = edges_viz_att, defaultedgetype = "undirected")
    
    return(nodes_viz_att)
}


#' @title compute community abundances
#' @description Compute the sum of abundance of species in each community given the abundance matrix by species and the graph community
#' Also compute the p-value of difference of abundance given a contrast vector with to class (2 classes)
#' @param Nodes: The network node table with the a community attribution column
#' @param abundance: The Metagenomic species abundance table
#' @param species_taxo: The Metagenomic species taxonomy table
#' @param contrast: a boolean vector to form two group of samples. for each community the rank test difference of abundance p-value is calculated between the two groups.
#' @param community kind: the algorithm of used to compute the community : "spinglass_community", "walktrap_community"
#' @return a table of community abundance and composition
#' @export
Compute_community_abondance <- function(Nodes,
                                        abundance,
                                        species_taxo,
                                        contrast=NULL,
                                        community_kind="walktrap_community")
{
# Beware that the row names of the node table must be species name
    abundance_by_comm <- data.frame(row.names = c( 'Community_composition',"Community_pval", colnames(abundance), "community_index"))
                                      
    abundance = as.data.frame(abundance)
    names(abundance) = sapply(names(abundance), remove_xml_char)

    for(spc in unique(Nodes[,community_kind]))
    {
        if (!is.na(spc)){
            MGS_in_c <- rownames(Nodes[which(Nodes[ ,community_kind]==spc),])
        
            
            if(length(MGS_in_c) > 1)
            {
                
                print('#############################')
                print('#######  Community  #########')
                                
                MGS_max <- MGS_in_c[which.max(apply(abundance[MGS_in_c, ], 1, mean, na.rm=TRUE))]
                c <- sort(apply(abundance[MGS_in_c, ], 1, mean), decreasing=TRUE)
                                
                MGS_rpr <- get_rpr(c)
                print(MGS_rpr)
                abundance_by_comm[MGS_rpr] = 0
                abundance_by_comm["community_index", MGS_rpr] = spc
                #_ get track of the integer index of the community for further calculation
                
                abundance_by_comm[setdiff(row.names(abundance_by_comm), c('Community_composition', "Community_pval", "community_index") ), MGS_rpr] <- apply(abundance[MGS_in_c,], 2, sum, na.rm=TRUE)
                
                abundance_by_comm[1, MGS_rpr] <- paste(MGS_in_c, collapse='-')
                
                print("Community number")
                print(spc)
                print("represented by:")
                 print(MGS_rpr)
                 print("complete list")
                 print(MGS_in_c)
                if(!is.null(contrast)){  
                    ab_h <- as.numeric(as.character(abundance_by_comm[(which(contrast)+2),  MGS_rpr]))
                    ab_l <- as.numeric(as.character(abundance_by_comm[(which(!contrast)+2), MGS_rpr]))
                
                    t <- wilcox.test(ab_h,ab_l) # log10(ab_h+10^(-8)), log10(ab_l+10^(-8)))
                                        # p_val_by_comm[MGS_rpr] <- signif(t$p.value, 3)
                 #   print(t)
                                        # print(wilcox.test(abundance[which(contrast), MGS_max],abundance[which(contrast),MGS_max]))
                                        # if(t$p.value < 0.05){ sign_com <- c(sign_com, MGS_rpr)}
                    abundance_by_comm["Community_pval", MGS_rpr] <- t$p.value
                }       
            }
        }
    }
    print("Community abundance is computed")
    return(abundance_by_comm)
}







#' generate_graph_from_tables
#' create a gephi graph from a Nodes table and Glasso model object
generate_graph_from_tables <- function(fout, nodes_viz_att, fgraph_model, variability_treshold=NULL) 
{
    se.gl = get(load(fgraph_model))
    if(is.null(variability_treshold))
    {
        icov_mat = se.gl$opt.icov
    }
    else
    {
       id_lbd <- max(which(se.gl$variability < variability_treshold))
       var_graph <- se.gl$variability[id_lbd]
       lbd <- se.gl$lambda[id_lbd]
       print(lbd)
       icov_mat <- se.gl$icov[[id_lbd]] 
       print("|||||||||||||| Graph parameters |||||||||||||||" )
       print(paste("Lambda:",as.character(lbd)))
       print(paste("Variability: " , as.character(var_graph)))
       print('|||||||||||||||||||||||||||||||||||||||||||||')
    }

        icov_mat = as.matrix(icov_mat)
    diag(icov_mat) = 0

    cnames <- sapply(colnames(se.gl$data), remove_xml_char)
    dimnames(icov_mat) = list(cnames, cnames)
    
    ig.gl <- graph.adjacency(Matrix(-icov_mat), weighted=TRUE, mode='min')
    E(ig.gl)[weight < 0]$color <- adjustcolor('coral', alpha.f=0.6) 
    E(ig.gl)[weight > 0]$color <- adjustcolor('darkslateblue', alpha.f=0.6) 
    
# "Out" rule on a bacteria
    Nodes <- data.frame(ID = c(1:vcount(ig.gl)), NAME = V(ig.gl)$name)
    edges <- as.data.frame(get.edges(ig.gl, 1:ecount(ig.gl)))

    edges_viz_att <- data.frame(WGH = abs(E(ig.gl)$weight), COLOR = E(ig.gl)$color)
    
    write.gexf(output= paste(fout, '.gexf', sep=""), nodes = Nodes, edges = edges, edgesWeight =round(abs(E(ig.gl)$weight)*100+1), nodesAtt = nodes_viz_att, edgesAtt = edges_viz_att, defaultedgetype = "undirected")
}


                                       #grepl("^[^_]+_2",s
get_rpr <- function(c)
{ 
    s= ""
    i = 1
    for(n in names(c))
    {
        if(!grepl("\\d",n)){
            s=paste(s,n)
            i = i+1
        }
        if(i>2) break
    }    
    if(nchar(s)==0)
    {s=names(c)[1]}
    return(s)
}



#' Computing spinglass communities and their modularity for a range of number of spin
#'  We then retrieve the optimal number of spin according to modularity
opt_spinglass_com <- function(con.grph, mod_rep, nspins)
{
    print('#####################################')
    print('#Optimization of spinglass community#')
    print("#####################################")
    seq_spin <- seq(nspins-10, nspins+10)
    mod_spin <- rep(0, length(seq_spin))
    err_mod <- rep(0, length(seq_spin))
    i <- 1
    for(ns in seq_spin)
    {
        mod_tmp <- rep(0, mod_rep)
        for(k in 1:mod_rep)
        {
            spin_com <- spinglass.community(con.grph, update.rule='config', implementation='neg', weights=E(con.grph)$weight, gamma=1, gamma.minus=0,spins=ns)
            # spin_com <- spinglass.community(pos.grph, weights = E(pos.grph)$weight, spins=ns)
            #print(table((E(pos.grph)$weight > 0)))
            # mod_tmp[k] <- modularity(spin_com, weights = E(pos.grph)$weight)
            mod_tmp[k] <- modularity(spin_com, weights=E(con.grph)$weight)
        }
        mod_spin[i]<-mean(mod_tmp)
        #print("number of spin:")
        #print(ns)
        #print("modularity")
        #print(mod_spin[i])
        #print("standart deviation on modularity:")
        #print(sd(mod_tmp))
        err_mod[i] <- sd(mod_tmp)/sqrt(mod_rep)
        i = i+1
    }
    
    p = ggplot2::qplot(seq_spin,mod_spin) + geom_errorbar(aes(x=seq_spin, ymin=mod_spin-err_mod, ymax=mod_spin+err_mod), width=0.25)
    #print(p)
    best_ns <-  seq_spin[which.max(mod_spin)]
    # spin_com <- spinglass.community(pos.grph, weights = E(pos.grph)$weight, spins=best_ns)
    # compute spinglass com for the best setting
    spin_com <- spinglass.community(con.grph, update.rule='config',implementation = 'neg', weights = E(con.grph)$weight, spins=best_ns, gamma=1, gamma.minus=0)
    print("best number of spin:")
    print(best_ns)
    print("best modularity")
    print(modularity(spin_com), weights =  E(con.grph)$weight)
    return(spin_com)
}


opt_walktrap <- function(pos.grph, lpath_range=seq(2,25,1), mod_rep=10)
{
    print('#####################################')
    print('# Optimization of walktrap community#')
    print("#####################################")
    mod_spin <- rep(0, length(lpath_range))
    err_mod <- rep(0, length(lpath_range))
    i <- 1
    for(ns in lpath_range)
    {
        walk_com <- walktrap.community(pos.grph, weight = E(pos.grph)$weight, steps=ns)
        mod_spin[i] <- modularity(walk_com, weights=E(pos.grph)$weight)
                
        # print("number of steps:")
        # print(ns)       
        # print("modularity")
        # print(mod_spin[i])
        
        i = i+1
    }
    
    p = qplot(lpath_range,mod_spin) #+ geom_errorbar(aes(x=lpath_range, ymin=mod_spin-err_mod, ymax=mod_spin+err_mod), width=0.25)
                                        #  print(p)
    best_ns <-  lpath_range[which.max(mod_spin)]
                                        # compute  walktrap community for the best setting in term of modularity
    walk_com <-  walktrap.community(pos.grph, weight = E(pos.grph)$weight, steps=ns)
  
    print("best number of steps:")
    print(best_ns)
    print("best modularity")
    print(modularity(walk_com), weights =  E(pos.grph)$weight)
    return(walk_com)

}


remove_xml_char <- function(x){return(gsub('[&"<>/-]', ' ', x))}

get_abundance_by_taxo_species <- function(abundance, taxo)
{    
    taxo[which(as.character(taxo$Taxo_level) !='species' & !is.na(taxo$Taxo_level)),'annot'] = paste(taxo[which(as.character(taxo$Taxo_level) !='species' & !is.na(taxo$Taxo_level)),'annot'], rownames(taxo)[which(as.character(taxo$Taxo_level) !='species' & !is.na(taxo$Taxo_level))])
    taxo[is.na(taxo$Taxo_level),'annot'] = rownames(taxo)[is.na(taxo$Taxo_level)]

    MGS_by_taxo_species <- as.data.frame(abundance)
    MGS_by_taxo_species[,'species'] <- as.character(taxo[rownames(MGS_by_taxo_species),"annot"])

    MGS_by_taxo_species[MGS_by_taxo_species$species!="unclassified", ]
    MGS_by_taxo_species <- data.table::data.table(MGS_by_taxo_species[MGS_by_taxo_species$species!="unclassified", ])
    MGS_by_taxo_species <- MGS_by_taxo_species[, lapply(.SD,mean), by =species]

    setkey(MGS_by_taxo_species, species)
    MGS_by_taxo_species <- data.frame(MGS_by_taxo_species)
    rownames(MGS_by_taxo_species) <- MGS_by_taxo_species$species
    MGS_by_taxo_species <- subset(MGS_by_taxo_species, select = - species)
    rnames <- sapply(rownames(MGS_by_taxo_species), remove_xml_char)
    rownames(MGS_by_taxo_species) <- rnames
    return(MGS_by_taxo_species)
}

generate_annotation <- function(taxo)
{
    
    taxo[which(as.character(taxo$Taxo_level) !='species' & !is.na(taxo$Taxo_level)),'annot'] = paste(
        taxo[which(as.character(taxo$Taxo_level) !='species' & !is.na(taxo$Taxo_level)),'annot'],
        rownames(taxo)[which(as.character(taxo$Taxo_level) !='species' & !is.na(taxo$Taxo_level))])
    taxo[is.na(taxo$Taxo_level), 'annot'] = rownames(taxo)[is.na(taxo$Taxo_level)]
    taxo$annot = sapply(taxo$annot, remove_xml_char)
    
    return(taxo)
}

generate_annotation_CAG_level <- function(taxo)
{
    
    taxo[,'annot'] = paste(taxo[,'annot'], rownames(taxo))
    taxo$annot = sapply(taxo$annot, remove_xml_char)
    #row.names(taxo) = taxo$annot
    return(taxo)
}



get_taxo_by_species <- function(taxo)
{
    taxo <- generate_annotation(taxo)
    species_taxo <- data.table::data.table(taxo)
    setkey(species_taxo, "annot")

    species_taxo <- unique(species_taxo[,.(annot,species,genus,family,order, class, phylum, superkingdom)])

    species_taxo <- as.data.frame(species_taxo)
    species_taxo <- species_taxo[which(!is.na(species_taxo$annot)),]

    rownames(species_taxo) <- species_taxo$annot  
    rnames <- sapply(rownames(species_taxo), remove_xml_char)

    rownames(species_taxo) <- rnames
    return(as.data.frame(species_taxo))
}


mean_when_present <- function(x)
{

    return(mean(x[x> 10^-7]))
}

remove_first <- function(x)
{
    if(substr(x, 1, 1)=='X')
    {
        x = substr(x, 3, nchar(x))
    }
    return(substr(x, 1, 23))
}

Tag_community_names <-function(Com_tab, Nodes)
{
    code = data.frame(row.names = as.numeric(unlist(Com_tab['community_index',])), Community = sapply(names(Com_tab), remove_first))
    Nodes$walktrap_community = code[as.character(Nodes$walktrap_community),] # Cast in character to avoid indexing by position
    return(Nodes)
}


#' overload '+' operator to allow character strings concatenation
`+` <- function(e1, e2){ 
	if (is.character(e1) && is.character(e2)) { 
		paste(e1,e2,sep="") 
	}
	else { 
		base::`+`(e1,e2) 
	} 
}
