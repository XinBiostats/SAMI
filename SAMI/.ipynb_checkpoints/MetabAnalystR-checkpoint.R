.set.mSet <- function(mSetObj=NA){
  return(mSetObj);
}

.get.mSet <- function(mSetObj=NA){
  return(mSetObj);
}

InitDataObjects <- function(data.type, anal.type, paired=FALSE){

  dataSet <- list();
  dataSet$type <- data.type;
  dataSet$design.type <- "regular"; # one factor to two factor
  dataSet$cls.type <- "disc"; # default until specified otherwise
  dataSet$format <- "rowu";
  dataSet$paired <- paired;
  dataSet$pair.checked <- FALSE;
  analSet <- list();
  analSet$type <- anal.type;
  Sys.setenv("OMP_NUM_THREADS" = 2); # to control parallel computing for some packages
  Sys.setenv("OPENBLAS_NUM_THREADS" = 2);
  mSetObj <- list();
  mSetObj$dataSet <- dataSet;
  mSetObj$analSet <- analSet;
  mSetObj$imgSet <- list();
  mSetObj$msgSet <- list(); # store various message during data processing
  mSetObj$msgSet$msg.vec <- vector(mode="character");     # store error messages
  mSetObj$cmdSet <- vector(mode="character"); # store R command
  metaboanalyst_env <<- new.env(); # init a marker env for raw data processing
  
  if (anal.type == "mummichog") {
    # Define this parameter set to avoid global variable
    # Author: Zhiqiang
    mSetObj$paramSet$mumRT <- NA;
    mSetObj$paramSet$mumRT.type <- NA;
    mSetObj$paramSet$version <- NA;
    mSetObj$paramSet$mumDataContainsPval <- 1;
    mSetObj$paramSet$mode <- NA;
    mSetObj$paramSet$adducts <- NA;
    mSetObj$paramSet$peakFormat <- "mpt";
  } else if (anal.type == "metapaths") {
    # Define this parameter set to avoid global variable
    # Author: Zhiqiang
    paramSet <- list();
    paramSet$mumRT <- NA;
    paramSet$mumRT.type <- NA;
    paramSet$version <- NA;
    paramSet$mumDataContainsPval <- 1;
    paramSet$mode <- NA;
    paramSet$adducts <- NA;
    paramSet$peakFormat <- "mpt";
    paramSet$metaNum <- 0;
    mSetObj$paramSet <- paramSet;
    # This is an empty paramSet, and will be copied for multiple datasets
    dataNMs <- names(mSetObj)[grepl("MetaData",names(mSetObj))];
    if(length(dataNMs)>0){
      for(n in dataNMs){
        mSetObj[[n]] <- NULL;
      }
    }
  }
  
  anal.type <<- anal.type;
  print("MetaboAnalyst R objects initialized ...");
  return(.set.mSet(mSetObj));
}

Setup.MapData <- function(mSetObj=NA, qvec){
  mSetObj <- .get.mSet(mSetObj);
  mSetObj$dataSet$cmpd <- qvec;
  return(.set.mSet(mSetObj));
}

CrossReferencing <- function(mSetObj=NA, q.type, hmdb=T, pubchem=T, 
                             chebi=F, kegg=T, metlin=F, lipid=F){
  
  mSetObj <- .get.mSet(mSetObj);
  
  # record the filter for 8 major databases
  mSetObj$return.cols <- c(hmdb, pubchem, chebi, kegg, metlin);
  mSetObj$lipid.feats <- lipid
  
  # record all the data
  if(!exists("name.map", where = mSetObj)){
    mSetObj$name.map <- list();
  }
  
  # distribute job
  mSetObj$dataSet$q.type <- q.type;
  
  mSetObj <- MetaboliteMappingExact(mSetObj, q.type, lipid);

  # do some sanity check
  todo.inx <- which(is.na(mSetObj$name.map$hit.inx));
  if(length(mSetObj$name.map$hit.inx) == 0){
    mSetObj$msgSet$nmcheck.msg <- c(0, "No hits found for the given compound ID. Please make 
                                    sure that correct compound IDs or common compound names are used.");
  }else if(length(todo.inx)/length(mSetObj$name.map$hit.inx) > 0.5){
    mSetObj$msgSet$nmcheck.msg <- c(0, "Over half of the compound IDs could not be matched to our database. Please make 
                                    sure that correct compound IDs or common compound names are used.");
  }else if (length(todo.inx) > 15){
    mSetObj$msgSet$nmcheck.msg <- c(2, "There are >15 compounds without matches. You can either proceed or if necessary, update these compound IDs and upload again.");        
  }else{
    mSetObj$msgSet$nmcheck.msg <- c(1, "Name matching OK, please inspect (and manual correct) the results then proceed.");   
  }
  

  print(mSetObj$msgSet$nmcheck.msg)
  
  if(length(todo.inx) == length(mSetObj$name.map$hit.inx)){
    AddErrMsg("Name matching failed! Please make sure that correct standardized feature names are used!")
    return(0)
  }
  
  return(.set.mSet(mSetObj));
}



MetaboliteMappingExact <- function(mSetObj=NA, q.type, lipid = F){
  
  mSetObj <- .get.mSet(mSetObj);
  
  if(lipid & anal.type == "msetqea"){
    qvec <- names(mSet$dataSet$url.var.nms);
  }else{
    qvec <- mSetObj$dataSet$cmpd;
  }
  
  # variables to record results
  hit.inx <- vector(mode='numeric', length=length(qvec)); # record hit index, initial 0
  names(hit.inx) <- qvec;
  match.values <- vector(mode='character', length=length(qvec)); # the best matched values (hit names), initial ""
  match.state <- vector(mode='numeric', length=length(qvec));  # match status - 0, no match; 1, exact match; initial 0 
  
  if(anal.type %in% c("msetora", "msetssp", "msetqea") & lipid){
    cmpd.db <- qs::qread("../lib/lipid_compound_db.qs");
  }else if(anal.type == "utils"){
    cmpd.db <- qs::qread("../lib/master_compound_db.qs");
  }else{
    cmpd.db <- qs::qread("../lib/compound_db.qs");
  }
  
  if(q.type == "hmdb"){
    n <- 5 # Number of digits for V3 of HMDB
    hmdb.digits <- as.vector(sapply(cmpd.db$hmdb, function(x) strsplit(x, "HMDB", fixed=TRUE)[[1]][2]))
    hmdb.v3.ids <- paste0("HMDB", substr(hmdb.digits, nchar(hmdb.digits)-n+1, nchar(hmdb.digits)))
    hit.inx.v3 <- match(tolower(qvec), tolower(hmdb.v3.ids));
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$hmdb));
    hit.inx[is.na(hit.inx)] <- hit.inx.v3[is.na(hit.inx)]
    match.values <- cmpd.db$name[hit.inx];
    match.state[!is.na(hit.inx)] <- 1;
  }else if(q.type == "pubchem"){
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$pubchem));
    match.values <- cmpd.db$name[hit.inx];
    match.state[!is.na(hit.inx)] <- 1;
  }else if(q.type == "chebi"){
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$chebi));
    match.values <- cmpd.db$name[hit.inx];
    match.state[!is.na(hit.inx)] <- 1;
  }else if(q.type == "metlin"){
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$metlin));
    match.values <- cmpd.db$name[hit.inx];
    match.state[!is.na(hit.inx)] <- 1;
  }else if(q.type == "kegg"){
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$kegg));
    #hit.inx2 <- match(tolower(qvec), rev(tolower(cmpd.db$kegg)));
    
    # unique hits
    #nonuniq.hits <- hit.inx + hit.inx2 != nrow(cmpd.db) + 1;
    #hit.inx[nonuniq.hits] <- NA;
    match.values <- cmpd.db$name[hit.inx];
    match.state[!is.na(hit.inx)] <- 1;
    
  }else if(q.type == "name"){
    # first find exact match to the common compound names
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$name));
    match.values <- cmpd.db$name[hit.inx];
    match.state[!is.na(hit.inx)] <- 1;
    
    # then try to find exact match to synonyms for the remaining unmatched query names one by one
    if(anal.type %in% c("msetora", "msetssp", "msetqea") & lipid){
      syn.db <- qs::qread("../lib/lipid_syn_nms.qs")
    }else if(anal.type == "utils"){
      syn.db <- qs::qread("../lib/master_syn_nms.qs")
    }else{
      syn.db <- qs::qread("../lib/syn_nms.qs")
    }
    
    syns.list <-  syn.db$syns.list;
    todo.inx <- which(is.na(hit.inx));
    
    if(length(todo.inx) > 0) {
      for(i in 1:length(syns.list)){
        syns <-  syns.list[[i]];
        hitInx <- match(tolower(qvec[todo.inx]), tolower(syns));
        
        hitPos <- which(!is.na(hitInx));
        if(length(hitPos)>0){
          # record matched ones
          orig.inx<-todo.inx[hitPos];
          hit.inx[orig.inx] <- i;                  
          # match.values[orig.inx] <- syns[hitInx[hitPos]];  # show matched synnames
          match.values[orig.inx] <- cmpd.db$name[i];    # show common name
          match.state[orig.inx] <- 1;
          
          # update unmatched list
          todo.inx<-todo.inx[is.na(hitInx)];
        }
        if(length(todo.inx) == 0) break;
      }
    }
  } else {
    print(paste("Unknown compound ID type:", q.type));save(qvec, file = "qvec__checking.rda");save(cmpd.db, file = "cmpd.db__checking.rda")
    # guess a mix of kegg and hmdb ids
    
    n <- 5 # Number of digits for V3 of HMDB
    hmdb.digits <- as.vector(sapply(cmpd.db$hmdb, function(x) strsplit(x, "HMDB", fixed=TRUE)[[1]][2]))
    hmdb.v3.ids <- paste0("HMDB", substr(hmdb.digits, nchar(hmdb.digits)-n+1, nchar(hmdb.digits)))
    hit.inx.v3 <- match(tolower(qvec), tolower(hmdb.v3.ids));
    hit.inx <- match(tolower(qvec), tolower(cmpd.db$hmdb));
    hit.inx[is.na(hit.inx)] <- hit.inx.v3[is.na(hit.inx)]
    
    #    hit.inx <- match(tolower(qvec), tolower(cmpd.db$hmdb));
    hit.inx2 <- match(tolower(qvec), tolower(cmpd.db$kegg));
    nohmdbInx <- is.na(hit.inx);
    hit.inx[nohmdbInx]<-hit.inx2[nohmdbInx]
    match.values <- cmpd.db$name[hit.inx];
    match.state[!is.na(hit.inx)] <- 1;
    
  }
  # empty memory
  gc();
  
  mSetObj$name.map$query.vec <- qvec; 
  mSetObj$name.map$hit.inx <- hit.inx;
  mSetObj$name.map$hit.values <- match.values;
  mSetObj$name.map$match.state <- match.state;
  
  return(.set.mSet(mSetObj));
}

SetMetabolomeFilter<-function(mSetObj=NA, TorF){
  
  mSetObj <- .get.mSet(mSetObj);
  
  mSetObj$api$filter <- TorF
  
  mSetObj$dataSet$use.metabo.filter <- TorF;
  return(.set.mSet(mSetObj));
}


SetCurrentMsetLib <- function(mSetObj=NA, libname, excludeNum=0){
  
  mSetObj <- .get.mSet(mSetObj);
  
  if(libname=="self"){
    ms.list <- mSetObj$dataSet$user.mset;
    ms.list <- lapply(ms.list, function(x) unique(unlist(strsplit(x, "; ", fixed=TRUE))));
    current.msetlib <- vector("list", 3)
    names(current.msetlib) <- c("name", "member", "reference")
    mSetObj$analSet$msetlibname <- libname;
  } else {
    if(grepl("kegg", libname)){ # api only for KEGG msets
      mSetObj$api$libname <- libname
      mSetObj$api$excludeNum = excludeNum
      mSetObj$analSet$msetlibname <- libname
      return(.set.mSet(mSetObj));
    }
    
    # feature enhancement https://omicsforum.ca/t/error-in-setcurrentmsetlib-function-in-r/2058
    if(!exists("current.msetlib") || is.null(mSetObj$analSet$msetlibname) || mSetObj$analSet$msetlibname != libname) {
      destfile <- paste('../lib/',libname, ".qs", sep = "");
      my.qs <- paste("https://www.metaboanalyst.ca/resources/libs/msets/", destfile, sep="");
      if(!file.exists(destfile)){
        download.file(my.qs, destfile, method = "curl");
      }
      current.msetlib <- qs::qread(destfile);
      mSetObj$analSet$msetlibname <- libname;
    }
    # create a named list, use the ids for list names
    # https://github.com/xia-lab/MetaboAnalystR/issues/172
    ms.list <- iconv(current.msetlib[, 3], from = 'utf8', to = 'utf8');
    ms.list <- lapply(ms.list, function(x) unique(unlist(strsplit(x, "; ", fixed=TRUE))));
    names(ms.list) <- current.msetlib[,2];
  }
  
  if(excludeNum > 0){
    cmpd.count <- lapply(ms.list, length);
    sel.inx <- cmpd.count >= excludeNum;
    ms.list <- ms.list[sel.inx];
    
    if(libname!="self"){
      current.msetlib <- current.msetlib[sel.inx,];
    }
  }
  
  # total uniq cmpds in the mset lib
  mSetObj$dataSet$uniq.count <- length(unique(unlist(ms.list, use.names = FALSE)));
  
  # update current.mset and push to global env
  current.msetlib$member <- ms.list;
  
  if(libname=="self"){
    current.msetlib$name <- names(ms.list)
    current.msetlib$reference <- rep("User-uploaded", length(ms.list))
  }
  
  current.msetlib <<- current.msetlib;
  return(.set.mSet(mSetObj));
}

CalculateHyperScore <- function(mSetObj=NA){
  
  mSetObj <- .get.mSet(mSetObj);
  
  # make a clean dataSet$cmpd data based on name mapping
  # only valid hmdb name will be used
  nm.map <- GetFinalNameMap(mSetObj);
  valid.inx <- !(is.na(nm.map$hmdb)| duplicated(nm.map$hmdb));
  ora.vec <- nm.map$hmdb[valid.inx];
  
  q.size<-length(ora.vec);
  
  if(all(is.na(ora.vec)) || q.size==0) {
    AddErrMsg("No valid HMDB compound names found!");
    return(0);
  }
  
  # move to api only if R package + KEGG msets
  if(grepl("kegg", mSetObj$analSet$msetlibname)){
    
    # make this lazy load
    if(!exists("my.hyperscore.kegg")){ # public web on same user dir
      .load.scripts.on.demand("util_api.Rc");    
    }
    
    mSetObj$api$oraVec <- ora.vec; 
    
    if(mSetObj$api$filter){
      mSetObj$api$filterData <- mSetObj$dataSet$metabo.filter.kegg
      toSend <- list(mSet = mSetObj, libNm = mSetObj$api$libname, filter = mSetObj$api$filter,
                     oraVec = mSetObj$api$oraVec, filterData = mSetObj$api$filterData,
                     excludeNum = mSetObj$api$excludeNum)
    }else{
      toSend <- list(mSet = mSetObj,libNm = mSetObj$api$libname, 
                     filter = mSetObj$api$filter, oraVec = mSetObj$api$oraVec, excludeNum = mSetObj$api$excludeNum)
    }
    saveRDS(toSend, "tosend.rds")
    return(my.hyperscore.kegg());
  }
  
  current.mset <- current.msetlib$member;
  
  # make a clean metabilite set based on reference metabolome filtering
  # also need to update ora.vec to the updated mset
  if(mSetObj$dataSet$use.metabo.filter && !is.null(mSetObj$dataSet$metabo.filter.hmdb)){
    current.mset <- lapply(current.mset, function(x){x[x %in% mSetObj$dataSet$metabo.filter.hmdb]})
    mSetObj$dataSet$filtered.mset <- current.mset;
  }
  
  set.size<-length(current.mset);
  if(set.size ==1){
    AddErrMsg("Cannot perform enrichment analysis on a single metabolite set!");
    return(0);
  }
  
  # now perform enrichment analysis
  
  # update data & parameters for ORA stats, based on suggestion
  # https://github.com/xia-lab/MetaboAnalystR/issues/168
  # https://github.com/xia-lab/MetaboAnalystR/issues/96
  # https://github.com/xia-lab/MetaboAnalystR/issues/34
  
  # the universe based on reference metabolome (should not matter, because it is already filtered based on the reference previously)
  
  my.univ <- unique(unlist(current.mset, use.names=FALSE));
  if(!is.null(mSetObj$dataSet$metabo.filter.hmdb)){
    my.univ <- unique(mSetObj$dataSet$metabo.filter.hmdb);
  }
  uniq.count <- length(my.univ);
  ora.vec <- ora.vec[ora.vec %in% my.univ];
  q.size <- length(ora.vec); 
  
  hits<-lapply(current.mset, function(x){x[x %in% ora.vec]});  
  hit.num<-unlist(lapply(hits, function(x) length(x)), use.names = FALSE);
  
  if(sum(hit.num>0)==0){
    AddErrMsg("No match was found to the selected metabolite set library!");
    return(0);
  }
  
  set.num<-unlist(lapply(current.mset, length), use.names = FALSE);
  
  # prepare for the result table
  res.mat<-matrix(NA, nrow=set.size, ncol=6);        
  rownames(res.mat)<-names(current.mset);
  colnames(res.mat)<-c("total", "expected", "hits", "Raw p", "Holm p", "FDR");
  for(i in 1:set.size){
    res.mat[i,1]<-set.num[i];
    res.mat[i,2]<-q.size*(set.num[i]/uniq.count);
    res.mat[i,3]<-hit.num[i];
    
    # use lower.tail = F for P(X>x)
    # phyper("# of white balls drawn", "# of white balls in the urn", "# of black balls in the urn", "# of balls drawn")
    res.mat[i,4]<-phyper(hit.num[i]-1, set.num[i], uniq.count-set.num[i], q.size, lower.tail=F);
  }
  
  # adjust for multiple testing problems
  res.mat[,5] <- p.adjust(res.mat[,4], "holm");
  res.mat[,6] <- p.adjust(res.mat[,4], "fdr");
  
  res.mat <- res.mat[hit.num>0,];
  
  ord.inx<-order(res.mat[,4]);
  mSetObj$analSet$ora.mat <- signif(res.mat[ord.inx,],3);
  mSetObj$analSet$ora.hits <- hits;
  
  #fast.write.csv(mSetObj$analSet$ora.mat, file="msea_ora_result.csv");
  return(.set.mSet(mSetObj));
}

GetFinalNameMap <- function(mSetObj=NA, lipid = FALSE){
  
  mSetObj <- .get.mSet(mSetObj);
  
  lipid = mSetObj$lipid.feats
  
  if (is.null(lipid)) {
    lipid = FALSE
  }
  
  hit.inx <- mSetObj$name.map$hit.inx;
  hit.values <- mSetObj$name.map$hit.values;
  match.state <- mSetObj$name.map$match.state;
  
  qvec <- mSetObj$dataSet$cmpd;
  nm.mat <- matrix(nrow=length(qvec), ncol=4);
  colnames(nm.mat) <- c("query", "hmdb",  "kegg", "hmdbid");
  
  if(anal.type %in% c("msetora", "msetssp", "msetqea") & lipid){
    cmpd.db <- qs::qread("../lib/lipid_compound_db.qs");
  }else if(anal.type == "utils"){
    cmpd.db <- qs::qread("../lib/master_compound_db.qs");
  }else{
    cmpd.db <- qs::qread("../lib/compound_db.qs");
  }
  
  for (i in 1:length(qvec)){
    
    hit <-cmpd.db[hit.inx[i], ,drop=FALSE];
    if(match.state[i]==0){
      hmdb.hit <- NA;
      hmdb.hit.id <- NA;
      kegg.hit <- NA;
    }else{
      hmdb.hit <- ifelse(nchar(hit.values[i])==0, NA, hit.values[i]);
      hmdb.hit.id <- ifelse(nchar(hit$hmdb_id)==0, NA, hit$hmdb_id);
      kegg.hit <- ifelse(nchar(hit$kegg_id)==0, NA, hit$kegg_id);
    }
    nm.mat[i, ]<-c(qvec[i], hmdb.hit, kegg.hit, hmdb.hit.id);
  }
  return(as.data.frame(nm.mat));
}

fast.write.csv <- function(dat, file, row.names=TRUE){
  tryCatch(
    {
      if(is.data.frame(dat)){
        # there is a rare bug in data.table (R 3.6) which kill the R process in some cases 
        data.table::fwrite(dat, file, row.names=row.names);
      }else{
        write.csv(dat, file, row.names=row.names);  
      }
    }, error=function(e){
      print(e);
      write.csv(dat, file, row.names=row.names);   
    }, warning=function(w){
      print(w);
      write.csv(dat, file, row.names=row.names); 
    });
}

PathwayEnrichment <- function(markers,clusters,lipid,merge){
  
  ora <- NULL
  
  if (lipid==FALSE){
    for (i in 1:length(clusters)){
      tryCatch(
        {input_list <- unique(markers[markers$cluster==clusters[i],]$Name)
        if (length(input_list)==0){
          next
        }
        rm(mSet,current.msetlib)
        mSet <- InitDataObjects("list","msetora",FALSE)
        mSet <- Setup.MapData(mSet,input_list)
        mSet <- CrossReferencing(mSet,"name")
        mSet <- SetMetabolomeFilter(mSet,F)
        mSet <- SetCurrentMsetLib(mSet, "smpdb_pathway",2)
        mSet <- CalculateHyperScore(mSet)
        metab_ora <- as.data.frame(mSet$analSet$ora.mat)
        metab_ora['cluster'] <- clusters[i]
        metab_ora['sample'] <- 'metab&glycan'
        metab_ora['pathway'] <- row.names(metab_ora)
        row.names(metab_ora) <- NULL
        
        ora <- rbind(ora,metab_ora)
        },error=function(err){
          cat('########## In cluster',clusters[i],"An Error occured:",conditionMessage(err),'\n')
        }
      )
    }
    if (merge==TRUE){
      write.csv(ora,'../results/pathway/ora_metgly_merge.csv')
    }else{
      write.csv(ora,'../results/pathway/ora_metgly.csv')
    }
  } else {
    for (i in 1:length(clusters)){
      tryCatch(
        {input_list <- markers[markers$cluster==clusters[i],]$Name
        if (length(input_list)==0){
          next
        }
        rm(mSet,current.msetlib)
        mSet <- InitDataObjects("list","msetora",FALSE)
        mSet <- Setup.MapData(mSet,input_list)
        mSet <- CrossReferencing(mSet,"name", lipid=TRUE)
        mSet <- SetMetabolomeFilter(mSet,F)
        mSet <- SetCurrentMsetLib(mSet, "sub_class",2)
        mSet <- CalculateHyperScore(mSet)
        lipid_ora <- as.data.frame(mSet$analSet$ora.mat)
        lipid_ora['cluster'] <- clusters[i]
        lipid_ora['sample'] <- 'lipidomics'
        lipid_ora['pathway'] <- row.names(lipid_ora)
        row.names(lipid_ora) <- NULL
        
        ora = rbind(ora,lipid_ora)
        },error=function(err){
          cat('########## In cluster',clusters[i],"An Error occured:",conditionMessage(err),'\n')
        }
      )
    }
    if (merge==TRUE){
      write.csv(ora,'../results/pathway/ora_lipid_merge.csv')
    }else{
      write.csv(ora,'../results/pathway/ora_lipid.csv')
    }
  }
}
