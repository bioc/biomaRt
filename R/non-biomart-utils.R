## The Ensembl BioMart isn't always fit for purpose, and some operations can also be achieved using the 
## Ensembl REST API.  This file contains some functions that may be useful to biomaRt users,
## but which don't actually interact with BioMart.

createNameToAliasMap <- function() {
  url <- "https://rest.ensembl.org/info/species?"

  req <- httr2::request(url) |>
    httr2::req_headers("Accept" = "application/json")
  json <- httr2::req_perform(req) |>
    httr2::resp_body_json()
  
  names <- vapply(json$species, FUN = \(x) { x$name }, FUN.VALUE = character(1))
  
  aliases <- lapply(json$species, 
                    FUN = \(x) { 
                      all <- c(x$display_name, unlist(x$aliases), x$common_name) 
                      return(unique(tolower(all)))
                    })
  
  names(aliases) <- names
  return(aliases)
}

pickReferenceStrain <- function(genomes_to_choose_from) {
  url <- "https://rest.ensembl.org/info/species?"
  
  req <- httr2::request(url) |>
    httr2::req_headers("Accept" = "application/json")
  json <- req_perform(req) |>
    resp_body_json()
  
  names <- vapply(json$species, FUN = \(x) { x$name }, FUN.VALUE = character(1))
  
  idx <- which(names %in% genomes_to_choose_from)
  
  strains <- vapply(json$species[idx], 
                    FUN = \(x) { 
                      unique(tolower(x$strain))
                    },
                    FUN.VALUE = character(1))
  
  idx2 <- which(strains == "reference")
  if(length(idx2) == 0)
    idx2 <- which(grepl(pattern = "reference", x = strains))
  
  if(length(idx2) != 1) {
    stop("unable to determine reference strain")
  }
  
  return(names[idx][idx2])
}

pickReferenceStrain_msg <- function(input_term) {
  message("Your search term was ambigous and multiple strains matching '", input_term, "' were found.\n",
          "Selecting the reference genome for this organism.\n",
          "Use a more specific search term if this is inappropriate.")
}

## Given an input string, try to identify the genome name for the organism
findGenomeName <- function(input) {
  
  input <- tolower(input)
  aliases <- createNameToAliasMap()
  
  if(input %in% names(aliases)) {
    res <- input
  } else {
    search <- vapply(aliases, FUN = \(x, input) { 
      any(grepl(pattern = input, x = x))
    }, input = paste0("^", input, "$"), 
    FUN.VALUE = logical(1))
    if(any(search) && (sum(search) == 1)) {
      res <- names(aliases)[which(search)]
    } else if (any(search) && (sum(search) > 1)) {
      res <- pickReferenceStrain(genomes_to_choose_from = names(aliases)[which(search)])
      pickReferenceStrain_msg(input)
    } else {
      stop('Unable to match the search term to a genome')
    }
  }
  return(res)
}

## use the Ensembl Rest API to get the current Ensembl version
getCurrentEnsemblRelease <- function() {
  
  req <- httr2::request("https://rest.ensembl.org/info/data/") |>
    httr2::req_headers("Accept" = "application/json")
  
  version <- req |> 
    httr2::req_perform() |>
    httr2::resp_body_json() |>
    unlist() |>
    as.integer()
  
  return(version)
}

## list files in Ensembl FTP TSV directory
## this is used to identify if there is an ensembl<->entrez mapping file
listFilesInEnsemblFTP <- function(species, release, dir) {
  
  loadNamespace("curl")
  
  ftp_dir <- sprintf("ftp://ftp.ensembl.org/pub/release-%s/%s/%s/", release, dir, species)
  
  list_files <- curl::new_handle()
  curl::handle_setopt(list_files, ftp_use_epsv = TRUE, dirlistonly = TRUE)
  con <- curl::curl(url = ftp_dir, open = "r", handle = list_files)
  files <- readLines(con)
  close(con)
  
  return(files)
  
}

## helper to make sure that we *have* entrez gene IDs to map to at ensembl!
.ensemblMapsToEntrezId <- function(species, release=NULL) {
    files <- listFilesInEnsemblFTP(species = species, release = release, dir = 'tsv')
    res <- any(grepl('entrez.tsv.gz$', files))
    return(res)
}

.shrinkDatasetName <- function(input) {

  parts <- stringr::str_split_fixed(input, pattern = "_", n = Inf)[1,]
  short <- paste(stringr::str_sub(head(parts, -1), 1, 1), collapse = "")
  final <- paste0(short, tail(parts, 1))
  return(final)
}

getHomologs <- function(ensembl_gene_ids, species_from, species_to) {
  
  from <- findGenomeName(species_from)
  to <- findGenomeName(species_to)
  
  dataset_from <- paste0(
    .shrinkDatasetName(from),
    "_gene_ensembl"
  )
  
  homolog_attribute <- paste0(
    .shrinkDatasetName(to),
    "_homolog_ensembl_gene"
  )
  
  mart <- useEnsembl(biomart = "genes",
                     dataset = dataset_from)
  
  q1 <- getBM(attributes=c("ensembl_gene_id", homolog_attribute),
              mart = mart,
              values = ensembl_gene_ids, 
              filters= "ensembl_gene_id")
  q1
}
