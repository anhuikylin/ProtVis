# Reproducible B73/Y12 validation from the verified quantitative matrix.
# Usage: Rscript dev/validate_analysis_repairs.R matrix.csv terms.csv output_directory
# terms.csv must contain protein_id and term (one mapping per row).
# Requires jsonlite, data.table, tibble, limma and statmod.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3L) stop("Provide matrix.csv, terms.csv and output_directory.")
script <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1L])
root <- normalizePath(file.path(dirname(script), '..'))
matrix_path <- args[[1L]]; terms_path <- args[[2L]]; output_path <- args[[3L]]
files <- c('protvis_dataset.R','protvis_architecture.R','protvis_runs_v4.R',
'protvis_analysis_scale.R','protvis_checkpoint.R','protvis_pipeline.R',
'maxquant_preprocessing_defaults.R','aaa_sample_grouping.R','data_source_parsers.R',
'protvis_import.R','correct_noise.R','data_transformed.R','DEP_analysis.R',
'sage_search.R','sage_search_module.R','zzzz_normalization_methods.R',
'zzzzz_pipeline_normalization.R')
for (file in files) source(file.path(root,'R',file))

source(file.path(root, 'R', 'protvis_workflow_status.R'))
raw <- read.csv(matrix_path, check.names=FALSE)
flags <- c('Only identified by site','Reverse','Potential contaminant')
keep <- rowSums(sapply(raw[flags], function(x) !is.na(x) & trimws(x)=='+'))==0
raw <- raw[keep, ,drop=FALSE]
x <- as.matrix(raw[,grepl('^(B73|Y12)_',names(raw)),drop=FALSE])
storage.mode(x)<-'numeric'; x[x<=0]<-NA_real_;rownames(x)<-raw$Protein_ID
info <- data.frame(sample_id=colnames(x), group=sub('_[0-9]+$','',colnames(x)))
terms <- read.csv(terms_path)
# Protein groups map to their members once, then duplicate term assignments collapse.
members <- strsplit(rownames(x), ';', fixed=TRUE)
mapping <- do.call(rbind,lapply(seq_along(members), function(i) {
 t <- unique(terms$term[terms$protein_id %in% members[[i]]])
 if(!length(t)) return(NULL)
 data.frame(protein_id=rownames(x)[i],term=t)
}))
cat('Input:',nrow(x),'proteins;',ncol(x),'samples;',nrow(mapping),'mapped terms\n')
d <- create_protvis_dataset(x,sample_info=info,annotation=list(GO_annotation=mapping),metadata=list(source='verified quantitative matrix',expression_scale='raw'))
checkpoint <- output_path
dir.create(checkpoint, showWarnings=FALSE)
p <- run_protvis_pipeline(d, stages=c('noise_correction','transformation','imputation','normalization','dimensionality_reduction'),params=list(noise_correction=list(max_missing=.99),transformation=list(method='log2',pseudocount=1e-8),imputation=list(method='median'),normalization=list(method='median')), checkpoint_dir=checkpoint,stop_on_error=TRUE)
stopifnot(validate_protvis_dataset(p),identical(p$metadata$expression_scale,'log2'))
# The observed matrix used by Shiny's Recommended DEP has no imputed values.
observed <- log2(x+1e-8)
observed <- observed[rownames(p$expression_data),,drop=FALSE]
centered <- .protvis_dep_prepare_recommended_matrix(observed)
summary <- list()
for (stage in c('Root_VE','Root_V1.V2','Root_V4','Leaf_VE.V1.V2','Leaf_V4.V6.V8')) {
 g1 <- paste0('B73_',stage);g2 <- paste0('Y12_',stage)
 s1<-info$sample_id[info$group==g1];s2<-info$sample_id[info$group==g2]
 ids <- .protvis_dep_shared_ids(observed,s1,s2,mode='both_genotypes',min_detected=2L)
 presence <- .protvis_dep_presence_absence(observed,s1,s2,min_detected=2L)
 limma <- .protvis_dep_run_limma_recommended(centered[ids,,drop=FALSE],s1,s2,g1,g2)
 limma <- .protvis_dep_classify(limma,logfc=1,cutoff=.05,p_metric='adj.P.Val')
 # Standard ORA uses the proteins actually tested in this comparison.
 native <- create_protvis_dataset(centered[ids,,drop=FALSE], sample_info=info,
   annotation=p$annotation, metadata=list(expression_scale='log2'))
 native <- .protvis_enrichment(native, list(selected_ids=limma$ID[limma$regulation!='Not significant']))
 stopifnot(native$analysis_results$enrichment$status=='success',
   all(is.finite(native$analysis_results$enrichment$table$p_value)))
 saveRDS(native$analysis_results$enrichment, file.path(checkpoint,
   paste0('recommended_enrichment_',stage,'.rds')))
 # Compare effect signs on identical complete observed data, independently of p-value engine.
 complete_ids <- ids[rowSums(!is.finite(centered[ids,c(s1,s2),drop=FALSE]))==0]
 equal <- create_protvis_dataset(centered[complete_ids,c(s1,s2),drop=FALSE],sample_info=info[match(c(s1,s2),info$sample_id),],metadata=list(expression_scale='log2'))
 equal <- .protvis_differential_analysis(equal,list(group1=g1,group2=g2))
 tt <- equal$analysis_results$differential_analysis$table
 stopifnot(max(abs(tt$log2FC-limma$logFC[match(tt$protein_id,limma$ID)]))<1e-8)
 result <- run_protvis_pipeline(p,stages=c('differential_analysis','enrichment','network'),params=list(differential_analysis=list(group1=g1,group2=g2),network=list(max_nodes=50)),checkpoint_dir=file.path(checkpoint,stage),stop_on_error=TRUE)
 stopifnot(validate_protvis_dataset(result),result$analysis_results$enrichment$status %in% c('success','skipped'))
 saved <- list.files(file.path(checkpoint,stage),pattern='checkpoint_network.*rds',full.names=TRUE)
 restored <- readRDS(tail(saved,1))
 stopifnot(validate_protvis_dataset(restored))
 summary[[stage]]<-data.frame(comparison=paste(g1,g2,sep='_vs_'),tested=nrow(limma),quantitative_significant=sum(limma$regulation!='Not significant'),detected_only=nrow(presence),pipeline_tested=nrow(result$analysis_results$differential_analysis$table),enrichment_terms=if(is.null(result$analysis_results$enrichment$table)) 0L else nrow(result$analysis_results$enrichment$table),max_direction_error=max(abs(tt$log2FC-limma$logFC[match(tt$protein_id,limma$ID)])))
 cat(stage,':',nrow(limma),'tested,',summary[[stage]]$quantitative_significant,'significant,',nrow(presence),'detected only\n')
}
summary<-do.call(rbind,summary)
write.csv(summary,file.path(checkpoint,'real_data_summary.csv'),row.names=FALSE)
print(summary)
