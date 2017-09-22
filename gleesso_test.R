library(Gleesso)

data_folder = "."
model_folder = "."
graph_folder = "."

MGS_file = "./MGS_abundances.RDS"
taxo_file = "./MGS_taxonomy.RDS"

boot_vec = rep(TRUE, 200)

    
Gleesso_pipeline(data_folder,
                 MGS_file,
                 taxo_file,
                 model_folder,
                 graph_folder,
                 boot_vec,
                 tag_model= "Gleesso_test",
                 tag_graph="_0.05_",
                 lambda.min.ratio= 0.01,
                 nlambda=50,
                 variability_treshold=0.2,
                 analysis_step =0,
                 occurence_treshold=0.5)
