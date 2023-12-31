\name{getGene}
\alias{getGene}
\title{Retrieves gene annotation information given a vector of identifiers}
\description{This function retrieves gene annotations from Ensembl given a vector of identifiers.  Annotation includes chromsome name, band, start position, end position, gene description and gene symbol.  A wide variety of identifiers is available in Ensembl, these can be found with the listFilters function.}
\usage{getGene( id, type, mart)}

\arguments{
  \item{id}{vector of gene identifiers one wants to annotate}
  \item{type}{type of identifier, possible values can be obtained by the listFilters function.  Examples are  entrezgene_id, hgnc_symbol (for hugo gene symbol), ensembl_gene_id, unigene, agilentprobe, affy_hg_u133_plus_2,  refseq_dna, etc.}
  \item{mart}{object of class Mart, containing connections to the BioMart databases.  You can create such an object using the function useMart.}
}

\author{Steffen Durinck} 

\examples{

if(interactive()){

mart = useMart("ensembl", dataset="hsapiens_gene_ensembl")

#example using affy id

g = getGene( id = "1939_at", type = "affy_hg_u95av2", mart = mart)
show(g)

#example using Entrez Gene id

g = getGene( id = "100", type = "entrezgene_id", mart = mart)
show(g)
}
}
\keyword{methods}

