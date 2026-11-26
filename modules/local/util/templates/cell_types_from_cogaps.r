#!/usr/bin/env Rscript
library(CoGAPS)

cogaps_path <- '${cogaps_obj}'
process <- '${task.process}'

dir.create('${prefix}', showWarnings = FALSE, recursive = TRUE)
message("reading CoGAPS object from ", cogaps_path)
cogaps <- readRDS(cogaps_path)
cell_types <- cogaps@sampleFactors
write.csv(cell_types, file = '${prefix}/cogaps_cell_types.csv')
write.csv(cogaps@featureLoadings, file = '${prefix}/cogaps_patterns.csv')

#versions
message("reading session info")
sinfo <- sessionInfo()
versions <- lapply(sinfo[["otherPkgs"]], function(x) {sprintf("    %s: %s",x[["Package"]],x[["Version"]])})
versions[['R']] <- sprintf("    R: %s
",packageVersion("base"))
cat(paste0(process,":
"), file="versions.yml")
cat(unlist(versions), file="versions.yml", append=TRUE, sep="
")
