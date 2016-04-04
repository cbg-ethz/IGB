#rm(list=ls())

RTdrugs <- c("ABC", "X3TC", "AZT", "D4T", "DDI", "EFV", "NVP", "TDF", "DDC", "DLV", "FTC")
PRdrugs <- c("IDV", "LPV", "NFV", "RTV", "SQV", "TPV", "ATV")
AllHIVdrugnames <- HIVdrugnames <- c(RTdrugs, PRdrugs)


read_mutations_group <- function(ds_name)
{
  mutation_groups = list()
  
  con<- file(paste("ds_definitions/", ds_name, ".txt", sep=''), 'r') 
  drug_name = NA
  input <- readLines(con, n=1)
  
  while (length(input) >  0 ){ 
    
    line = unlist(strsplit(input, " "))
    
    if(length(line) == 1 && line[1] %in% AllHIVdrugnames)
    {
      if(is.na(drug_name) == FALSE)
        mutation_groups[[drug_name]] = mut_list
      
      drug_name = line[1]
      mut_list = list()
    } else if(length(line) > 1)
    {
      group_name = line[1]
      mutations = line[2:length(line)]
      mut_list[[group_name]] = c(mutations)
    }
    input <- readLines(con, n=1)
  }
  
  mutation_groups[[drug_name]] = mut_list
  close(con)
  mutation_groups
}

generate_X_vector <- function(mut_group, obs_mutations) {
  
  group_names <-   names(mut_group)
  
  # number of samples
  Xnew <- rep(0, length(group_names))
  names(Xnew) <- group_names
  
  g_name = group_names[11]
  for(g_name in group_names)
  {
    
    if(length(mut_group[[g_name]]) > 1)
    {
      mutations <- mut_group[[g_name]]
      if(any(mutations%in%obs_mutations) )
        Xnew[g_name] = 1
    }
    else {
      if(length(grep(mut_group[[g_name]] , obs_mutations )) != 0) {
        Xnew[g_name] = 1
      }
    }
  }
  Xnew
}

computeIGB <-function(obs_mutations, drug_name)
{   
  mut_groups <- read_mutations_group("StanfordSubset")
  mut_group <- mut_groups[[drug_name]]
  
  load(paste("icbn_models/GenBarDrug-", drug_name, ".Rdata", sep=""))
  X <- generate_X_vector(mut_group, obs_mutations) 

  index <- apply(GenBarDrug, 1, function(x) { all(x[1:length(X)] == X) })
  GenBarDrug[index, "GenBar"]
}

# Example
obs_mutations <- c("RT184V", "RT67N", "RT65R")
drug_name = "AZT"

computeIGB(obs_mutations, drug_name)
