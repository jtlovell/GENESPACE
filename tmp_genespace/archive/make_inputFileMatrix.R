#' @title Make input metadata for pipe_Diamond2MCScanX
#'
#' @description
#' \code{make_inputFileMatrix} Utility function to build metadata
#'
#' @param peptide.dir The location of the Diamond program executable
#' @param mappings.dir ID for genome 1
#' @param gff.dir ID for genome 2
#' @param mcscan.dir Peptide file for genome 1
#' @param ploidy.dict Named list of ploidies, where each list element is named
#' for the genomeID and the one-element numeric vector is the ploidy
#' @param abbrev.dict Same as ploidy.dict, but the content is an appreciation for
#' each genomeID
#' @param ref.id Identifier of the reference genome. This is the only genome that
#' blastst to itself.
#' @param ... Not currently in use
#' @details See pipe_Diamond2MCScanX for more information.
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' none yet
#' }
#' @export
make_inputFileMatrix<-function(peptide.dir,
                               remap.dir,
                               mappings.dir,
                               gff.dir,
                               mcscan.dir,
                               abbrev.dict,
                               ploidy.dict,
                               ref.id){
  ids = gsub(".pep.fa","",dir(peptide.dir), fixed = T)
  id.mat = expand.grid(ids,ids)
  id.mat = t(apply(id.mat,1,function(x) x[order(x)]))
  id.mat = cbind(ifelse(id.mat[,1]==ref.id, id.mat[,1], id.mat[,2]),
                 ifelse(id.mat[,1]==ref.id, id.mat[,2], id.mat[,1]))
  id.mat = data.frame(id.mat[!duplicated(id.mat),],
                      stringsAsFactors = F)
  colnames(id.mat)<-c("id1","id2")
  id.mat = id.mat[id.mat$id1 == ref.id | id.mat$id1 != id.mat$id2,]

  make_abbrev = function(x,abbrev.dict){
    for(i in 1:length(abbrev.dict)){
      x[x == names(abbrev.dict)[i]]<-abbrev.dict[[i]]
    }
    return(x)
  }

  make_ploidy = function(x,ploidy.dict){
    for(i in 1:length(ploidy.dict)){
      x[x == names(ploidy.dict)[i]]<-ploidy.dict[[i]]
    }
    return(x)
  }

  id.mat$ploidy1 = as.numeric(make_ploidy(id.mat$id1, ploidy.dict = ploidy.dict))
  id.mat$ploidy2 = as.numeric(make_ploidy(id.mat$id2, ploidy.dict = ploidy.dict))
  id.mat$abbrev1 = make_abbrev(id.mat$id1, abbrev.dict = abbrev.dict)
  id.mat$abbrev2 = make_abbrev(id.mat$id2, abbrev.dict = abbrev.dict)
  id.mat$pep1 = file.path(peptide.dir,paste0(id.mat$id1,".pep.fa"))
  id.mat$pep2 = file.path(peptide.dir,paste0(id.mat$id2,".pep.fa"))
  id.mat$gff1 = file.path(gff.dir, paste0(id.mat$id1,".gff3"))
  id.mat$gff2 = file.path(gff.dir, paste0(id.mat$id2,".gff3"))
  id.mat$remap.pep1 = file.path(remap.dir,paste0(id.mat$id1,".pep.fa"))
  id.mat$remap.pep2 = file.path(remap.dir,paste0(id.mat$id2,".pep.fa"))

  id.mat$id2 <- ifelse(id.mat$id1 == id.mat$id2, paste0(id.mat$id2,"_1"),id.mat$id2)
  id.mat$blast1 = file.path(mappings.dir, paste0(id.mat$id1,"_",id.mat$id2,".blast8"))
  id.mat$blast2 = file.path(mappings.dir, paste0(id.mat$id2,"_",id.mat$id1,".blast8"))

  # id.mat = id.mat[id.mat$id1 == ref.id,]
  return(id.mat)
}
