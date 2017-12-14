
suppressMessages(library(Gleesso))

taxo_9.9_file = "./data/taxo_hs_9.9_3463_cleaned_CAG_sup50_freeze1.rds"

mgs_igc_ch_da_file = './data/MGS_IGC_healthy_CH_Da.rds'

meta_IGC = read.csv("./data/meta_IGC_healthy.csv")

selection_vector = rep(TRUE, dim(readRDS(mgs_igc_ch_da_file))[2])

Gleesso_pipeline('./outputs/',
                                 mgs_igc_ch_da_file,
                                 taxo_9.9_file,
                                 './outputs/',
                                 './outputs/',
                                    selection_vector,
                                 "Gleesso_test",
                                 "default_param"
                                )

community_table = readRDS("./outputs/community_abondance_Gleesso_testdefault_param.rds")

head(community_table)

Gleesso_bootstrap(2, 0.8, "healthy_Da_Chi", 
                  NULL, "./outputs/",
                  "./outputs/", "./outputs/", 
                  "./data/MGS_IGC_healthy_CH_Da.rds",
                  "./data/taxo_hs_9.9_3463_cleaned_CAG_sup50_freeze1.rds",
                  stratifying_vector = as.character(meta_IGC[,'country'])
                     )

taxo_by_species = readRDS("./outputs/taxo_by_species.rds")

robust_tables = Robust_table_community("./outputs/",
                    "./IGC_community_100_graph",
                      taxo_by_species
                      )

saveRDS(robust_tables, "./Robust_tables_100_graph.rds")

p_value = community_contrast_dashboard(robust_tables$Robust_community_stability_0.6,
                                "country",
                            meta_IGS$country, 
                                      nrow=5
                                      )

options(repr.plot.height= 6, repr.plot.width= 10)
draw_community_total_abundance(abund = extract_community_abundance_table(robust_tables[["Robust_community_stability_0.6"]]))

options(repr.plot.height= 6, repr.plot.width= 10)
draw_community_total_species_count(robust_tables[["Robust_community_stability_0.6"]])

Nodes_all_samples = 
read.csv("./outputs/graph_healthy_Da_Chi_all_samples__nodes.csv", 
         row.names=1)

abund_by_species= readRDS("./outputs/abund_by_species.rds")

node_robust  = create_graph_robust_community_tags(
model_folder = "./outputs",
                                      "./Robust_tagged_100_graph",                                    
                                  abund_by_species,
                                      taxo_by_species,
                        "healthy_Da_Chi",
                                     robust_tables[["Robust_community_stability_0.6"]],
                            Nodes_all_samples
                                  )

com_tab = robust_tables[["Robust_community_stability_0.6"]]

abund_com = extract_community_abundance_table(com_tab)

trim_tab = community_taxa_abundance(Nodes =  node_robust,
                        abund_by_species,
                         "./images/Taxa_relative_abundance_",                         
                        community_kind = "additional_info"
                        )
