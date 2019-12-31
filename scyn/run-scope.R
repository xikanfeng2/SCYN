library(SCOPE)
library(WGSmapp)
args <- commandArgs(TRUE)
# bamFile <- list.files(args[1], pattern = paste("*", args[6],"bam$", sep=""))
bamFile <- list.files(args[1], pattern = args[8])
bamdir <- file.path(args[1], bamFile)
sampname_raw <- sapply(strsplit(bamFile, ".", fixed = TRUE), "[", 1)
bambedObj <- get_bam_bed(bamdir = bamdir, sampname = sampname_raw, resolution=strtoi(args[11]))
bamdir <- bambedObj$bamdir
sampname_raw <- bambedObj$sampname
ref_raw <- bambedObj$ref


if (args[9] == 'hg19') {
    data("mapp_hg19")
    mapp <- get_mapp(ref_raw)
    gc <- get_gc(ref_raw)
    values(ref_raw) <- cbind(values(ref_raw), DataFrame(gc, mapp))
}else{
    library(BSgenome.Hsapiens.UCSC.hg38)
    data("mapp_hg38")
    mapp <- get_mapp(ref_raw, hgref = "hg38")
    gc <- get_gc(ref_raw, genome = BSgenome.Hsapiens.UCSC.hg38)
    values(ref_raw) <- cbind(values(ref_raw), DataFrame(gc, mapp))
}


coverageObj <- get_coverage_scDNA(bambedObj, mapqthres = strtoi(args[10]), seq = args[7])
Y_raw <- coverageObj$Y

QCmetric_raw <- get_samp_QC(bambedObj)
qcObj <- perform_qc(Y_raw = Y_raw, 
    sampname_raw = sampname_raw, ref_raw = ref_raw, 
    QCmetric_raw = QCmetric_raw)
Y <- qcObj$Y
sampname <- qcObj$sampname
ref <- qcObj$ref
QCmetric <- qcObj$QCmetric

write.csv(Y, args[2])
write.csv(ref, args[3])

Gini <- get_gini(Y)
write.csv(Gini, args[4])
normObj.sim <- normalize_codex2_ns_noK(Y_qc = Y,
                                        gc_qc = ref$gc,
                                        norm_index = which(Gini<=0.12))

Yhat.noK.sim <- normObj.sim$Yhat
beta.hat.noK.sim <- normObj.sim$beta.hat
fGC.hat.noK.sim <- normObj.sim$fGC.hat
N.sim <- normObj.sim$N

# Ploidy initialization
ploidy.sim <- initialize_ploidy(Y = Y, Yhat = Yhat.noK.sim, ref = ref)
write.csv(ploidy.sim, args[5])
# If using high performance clusters, parallel computing is 
# easy and improves computational efficiency. Simply use 
# normalize_scope_foreach() instead of normalize_scope(). 
# All parameters are identical. 
normObj.scope.sim <- normalize_scope_foreach(Y_qc = Y, gc_qc = ref$gc,
    K = 1, ploidyInt = ploidy.sim,
    norm_index = which(Gini<=0.12), T = 1:7,
    beta0 = beta.hat.noK.sim, nCores = 2)
# normObj.scope.sim <- normalize_scope(Y_qc = Y_sim, gc_qc = ref_sim$gc,
#     K = 1, ploidyInt = ploidy.sim,
#     norm_index = which(Gini<=0.12), T = 1:7,
#     beta0 = beta.hat.noK.sim)
Yhat.sim <- normObj.scope.sim$Yhat[[which.max(normObj.scope.sim$BIC)]]
fGC.hat.sim <- normObj.scope.sim$fGC.hat[[which.max(normObj.scope.sim$BIC)]]

write.csv(Yhat.sim, args[6])

# chrs <- unique(as.character(seqnames(ref)))
# segment_cs <- vector('list',length = length(chrs))
# names(segment_cs) <- chrs
# for (chri in chrs) {
#     message('\n', chri, '\n')
#     segment_cs[[chri]] <- segment_CBScs(Y = Y,
#                                     Yhat = Yhat.sim,
#                                     sampname = colnames(Y),
#                                     ref = ref,
#                                     chr = chri,
#                                     mode = "integer", max.ns = 1)
# }
# iCN_sim <- do.call(rbind, lapply(segment_cs, function(z){z[["iCN"]]}))
# row.names(iCN_sim) <- row.names(Y)
# write.csv(iCN_sim, args[2])
