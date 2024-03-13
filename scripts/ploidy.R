
args <- commandArgs(trailingOnly = TRUE)

snakemake_dir<-args[1]
sample<-args[2]
h1.paf<-args[3]
h2.paf<-args[4]
out_pdf<-args[5]
out_summary<-args[6]

source(paste(snakemake_dir,'processAssemblyAlignments_functions.R', sep="/"))

cat(paste("\n##SAMPLE NAME:",sample,"\n"))

out <- asm2referenceCoverage(h1.paf = h1.paf, h2.paf = h2.paf, index = sample)

pdf (out_pdf)
plot(out$plot)
dev.off()
out$ploidy.ranges

cat(out$ploidy.summary, file=out_summary)

