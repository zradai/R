
#####################################################
### process kraken2 reports into community matrix ###
#####################################################

# Takes kraken2 report files and extracts species-level (or other, user-specified) taxa, into sample-by-taxa community matrix, 
# where the values are the raw (or standardized: only TMM for now) read counts.

# SEE REQUIRED PACKAGES BELOW ('list.of.packages')

# dir: path to report files
# min.percentage: percentage value of fragments covered by the clade rooted at this taxon, below which records are ignored
# min.cladereads: number of reads covered by the clade rooted at this taxon, below which records are ignored
# file.pattern: optional string used in 'list.files' as a pattern for the selection of report files
# taxrank: taxonomic rank for which to acquire clade read data; possible levels are:
## (U)nclassified, (R)oot, (D)omain, (K)ingdom, (P)hylum, (C)lass, (O)rder, (F)amily, (G)enus, (S)pecies
# lineage: name of taxon for which to constrain the retained findings;
## e.g. if 'lineage="bacteria"' then only bacterial taxa will be considered;
## can be any taxonomic label that can be found in the reports' 'taxLineage' column, but should be compatible with 'taxrank'
# nc: number of CPU threads to use in parallelization; default is '0' when the half of the available cores will be used

# TMM: trimmed mean of M values -> a read count normalization method (see https://doi.org/10.1186/gb-2010-11-3-r25)
# lib.sizes: a named vector of library sizes, with length of number of samples (report files), and names being the report files' names
## (can be omitted, but then the filtered (i.e. "after-criteria-matching") sample-wise sums of read counts are used, which is not recommended...)
# NOTE: edgeR is likely to throw errors for sparse and/or low-variability matrices -> a filtering option will be added later...

# when used on desktop PC / laptop be mindful of resource usage: memory usage at 16 cores was 3.8 GB on my laptop

#--- examples ---

# basic usage
kraken.to.matrix(dir = "/path/to/kraken2_reports")

# get only bacteria
kraken.to.matrix(dir = "kraken2_reports", lineage = "bacteria")

# get bacteria, with TMM-standardized read counts (from 10 reports)
libs = runif(10, 5e3, 1e6)
names(libs) = list.files("kraken_reports/")
kraken.to.matrix(dir = "kraken2_reports", lineage = "bacteria", TMM = TRUE, lib.sizes = libs)

# get only species-level bacterial records which contain at least 200 reads, using 12 CPU cores
kraken.to.matrix(dir = "kraken2_reports", min.cladereads = 200, taxrank = "S", lineage = "bacteria", nc = 12)
  

#--- function ---

kraken.to.matrix = function(dir=getwd(), min.percentage=0, min.cladereads=0, file.pattern=NULL, taxrank="S", lineage=NULL, nc=0, TMM=FALSE, lib.sizes=NULL){
  
  if(substr(dir, nchar(dir), nchar(dir))!="/"){dir = paste0(dir, "/")}
  
  list.of.packages <- c(
    "parallel",
    "foreach",
    "doParallel",
    "ranger",
    "palmerpenguins",
    "tidyverse",
    "kableExtra",
    "pavian",
    "edgeR"
  )
  for(package.i in list.of.packages){
    suppressPackageStartupMessages(
      library(
        package.i, 
        character.only = TRUE
      )
    )
  }
  
  if(nc==0){nc = floor(detectCores()/2)}else if(nc > detectCores()){
    nc = floor(detectCores()/2)
    warning("Specified 'nc' value is larger than the number of available cores - 'nc' set to ",nc,"!")
  }
  
  if(TMM & is.null(lib.sizes)){
    warning("No library size information was provided! Using sum of reads per sample...")
  }
  
  fls = list.files(path = dir, pattern = file.pattern)
  
  cl = makeCluster(nc)
  registerDoParallel(cl)
  fls = fls
  ee = foreach(
    f = fls
  ) %dopar% {
    
    if(file.size(paste0(dir, f))>0){
      
      kr.f = pavian::read_report2(myfile = paste0(dir,f))
      if(!is.null(lineage)){kr.f = kr.f[grepl(pattern = lineage, x = kr.f$taxLineage, ignore.case = TRUE),]}
      kr.f = kr.f[kr.f$taxRank==taxrank,]
      kr.f = kr.f[kr.f$percentage>=min.percentage & kr.f$cladeReads>=min.cladereads,]
      
      if(nrow(kr.f)>0){
        ss = kr.f$cladeReads
        names(ss) = sub(pattern = paste0(tolower(taxrank), "_"), replacement = "", x = kr.f$name)
        
        ll = list()
        ll[[f]] = ss
        
        ll
      }else{
        paste0("Report file '", f, "' contains no records with the current specifications, and therefore was omitted!")
      }
      
    }else{
      paste0("Report file '", f, "' was empty!")
    }
    
  }
  stopCluster(cl)
  
  ll = list()
  for(i in 1:length(ee)){
    if(!is.null(ee[[i]]) & is.list(ee[[i]])){
      ll[[names(ee[[i]])]] = ee[[i]][[1]]
    }else{
      warning(ee[[i]])
    }
  }
  
  taxa = sort(unique(unlist(lapply(ll, names))))
  
  mm = matrix(0, nrow = length(ll), ncol = length(taxa), dimnames = list(names(ll),taxa))
  for(i in names(ll)){
    mm[i,names(ll[[i]])] = ll[[i]]
  }
  
  
  if(TMM){
    
    mm.t = t(mm)
    if(is.null(lib.sizes)){
      lib.sizes = colSums(mm.t)
    }
    dge = edgeR::DGEList(counts = mm.t, lib.size = lib.sizes[colnames(mm.t)])
    tmm = edgeR::calcNormFactors(dge, method = "TMM")
    mm.t.tmm = edgeR::cpm(tmm)
    mm = t(mm.t.tmm)
    
  }
  
  return(mm)
  
}
