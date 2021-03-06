% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Bootstrap_community.R
\name{create_graph_robust_community_tags}
\alias{create_graph_robust_community_tags}
\title{create_graph_robust_community_tags}
\usage{
create_graph_robust_community_tags(model_folder, fout, abund_by_species,
  taxo_by_species, model_tag, Robust_table_community,
  Nodes_table_on_all_samples, variability_treshold = NULL)
}
\arguments{
\item{fout}{: where to save the csv tables with nodes with the robust community column}

\item{abund_by_species}{: abundance mean group by species or}

\item{taxo_by_species}{: taxo grouped at the species level}

\item{Robust_table_community}{: Community attribution}

\item{Nodes}{tables computed on all samples}
}
\description{
create a gephi file for nodes with a robust community attribution
}
