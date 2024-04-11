## The Ensembl BioMart isn't always fit for purpose, and some operations can also be achieved using the 
## Ensembl REST API.  This file contains some functions that may be useful to biomaRt users,
## but which don't actually interact with BioMart.

createNameToAliasMap <- function() {
  url <- "https://rest.ensembl.org/info/species?"

  req <- httr2::request(url) |>
    httr2::req_headers("Accept" = "application/json")
  json <- req_perform(req) |>
    resp_body_json()
  
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
  
  message("Your search term was ambigous and multiple strains matching this term were found\n",
          "Selecting the reference geneome for this organism.")
  
  return(names[idx][idx2])
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

