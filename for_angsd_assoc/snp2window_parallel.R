## load libraries ----
require(optparse)
require(parallel)
require(zoo)
require(data.table)


## parameters ----
option_list <- list(
                    make_option(c('-i','--input'),  action='store', default=NULL, help='input file'),
                    make_option(c('-w','--snpwin'), action='store', type='numeric', default=100, help='fixed number of SNPs per window'),
                    make_option(c('-s','--snpstp'), action='store', type='numeric', default=10, help='fixed number of SNPs per step'),
                    make_option(c('-t','--quant'),  action='store', type='numeric', default=0.990, help='percentile to define outliers threshold'),
                    make_option(c('-m','--maxlength'),  action='store', type='numeric', default=1000000, help='maximum physical length of SNP window [default 1Mb]'),
                    make_option(c('-c','--cores'),  action='store', type='numeric', default=1, help='number of cores to use'))

## parameters ----
opt <- parse_args(OptionParser(option_list = option_list))
file = opt$input
snp.win = opt$snpwin
snp.stp = opt$snpstp
thr = opt$quant
maxlen = opt$maxlength
no_cores = opt$cores

## input data ----
# chose file and do necessary transformations
print("READING INPUT FILE ...")
#inp = read.table(gzfile(file), header=TRUE)
inp = fread(file)
names(inp) = c("SNP", "CHR", "BP", "Pvalue")
print("TAKING A LOOK TO INPUT FILE ... CHECK IF OK")
head(inp)

# tranform chrX name
#inp$CHR = ifelse(inp$CHR == "chrXX", "22", as.character(inp$CHR) )

## calculate mean Pvalue in snp windows ----
print("CALCULATING MEAN Pvalue IN SNP BLOCKS ...")
inp.df = data.frame(chr = as.character(), sta = as.numeric(), end = as.numeric(), Pvalue = as.numeric())

cl <- makeCluster(no_cores)
source("snp2window.function.R")
clusterExport(cl, "inp")
clusterExport(cl, "snp.win")
clusterExport(cl, "snp.stp")
clusterExport(cl, "inp.df")
ptm <- proc.time()
inp.df = as.data.frame(do.call(rbind, parLapply(cl, unique(inp$CHR), snp2window)))
#inp.df = as.data.frame(do.call(rbind, parLapply(cl, c(4,20,21), snp2window)))
proc.time() - ptm
stopCluster(cl)

# calculate physical length (bp) of SNP blocks
inp.df$len = inp.df$end - inp.df$sta + 1
write.table(inp.df, paste0("snp",snp.win,"_",snp.stp,".2.Pvalue"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

## check length of SNP blocks ----
print("SUMMARY OF SNP WINDOWS PHYSICAL LENGTHS DISTRIBUTION:")
summary(inp.df$len)

## filter windows that are too large ----
print(paste0("FILTERING SNP WINDOWS LARGER THAN ", maxlen, "KB"))
inp.flt = subset(inp.df, len < maxlen)
jpeg(paste0("Pvalue-distribution_",snp.win,"_",snp.stp,"2..jpeg"), units = "cm", width=30, height=15, res = 150)
mn = floor(min(inp.flt$Pvalue))
mx = ceiling(max(inp.flt$Pvalue))
hist(inp.flt$Pvalue, breaks=(seq(0,mx,(mx-mn)/100)), xlab="Mean -log(p-value)", main="Mean -log(p-value) distribution")
dev.off()
print(summary(inp.flt$len))


## define merging function ----
mrg <- function(x) {
    y = head(x, n=1)
    for (r in 2:nrow(x)) {
        # last position of mrg
        l = nrow(y)
        # same CHRosome?
        if (x$chr[r] == y$chr[l]) {
            # overlap
            if (x$sta[r] <= y$end[l]) {
                y$end[l] = x$end[r]
            }
            # do not overlap
            else {y = rbind(y, x[r,])}
        }
        # not same CHRomosome
        else {
            print(x$chr[r])
            y = rbind(y, x[r,])}
    }
    return(y)
}

## define outliers ----
outlier = subset(inp.flt, Pvalue >= thr)
jpeg(paste0("Pvalue_SNPwindow-lengths_",snp.win,"_",snp.stp,".2.jpeg"), units = "cm", width=30, height=15, res = 150)
boxplot(outlier$len/1000, inp.flt$len/1000, horizontal = TRUE, names=c("Outlier Windows","All Windows"), xlab="SNP window physical lengths (kb)")
dev.off()

## find chromosome middle points
middle_bp = as.numeric()
last_bp = as.numeric()
for (i in unique(inp.flt$chr)) {
  print(i)
  sub1 = subset(inp.flt, chr == i)
  if (length(middle_bp) == 0) {
    middle_bp = nrow(sub1)/2
    last_bp = nrow(sub1)
    }
  else {
    mbp = nrow(sub1)/2 + last_bp
    middle_bp = c(middle_bp, mbp)
    last_bp = last_bp + nrow(sub1)
    }
}

## plot Pvalue ----
if (mx < 15) {mx = 15}
jpeg(paste0("Pvalue_manhattan_plot_",snp.win,"_",snp.stp,".2.jpeg"), units = "cm", width=30, height=10, res = 300)
plot(inp.flt$Pvalue, cex=0.5, ylim = c(mn,10), ylab="mean -log10(P)", xlab="", col=ifelse(inp.flt$chr %in% c(seq(2,22,2)), "gray60", "gray10"), xaxt='n', pch=19, bty="l", axes=0)
abline(h=thr, col="black", lty=2)
#axis(1, at=middle_bp[1:22], labels=unique(inp.flt$chr)[1:22], las=2)
axis(1, at=middle_bp[1:22], pos=0, labels=FALSE)
axis(2, labels=FALSE)
#text(x=middle_bp[1:22], y=-2, labels=unique(inp.flt$chr)[1:22], srt=45, adj=1, xpd=FALSE)

dev.off()

## outlier coordinates ----
out.mrg = mrg(outlier)
write.table(out.mrg[,c(1:3)], paste0("outliers_snp",snp.win,"_",snp.stp,".2.coord_biomart.txt"), col.names=FALSE, row.names=FALSE, sep=":", quote=FALSE)

## background coordinates
bg.mrg = mrg(inp.flt)
write.table(bg.mrg[,c(1:3)], paste0("background_snp",snp.win,"_",snp.stp,".2.coord_biomart.txt"), col.names=FALSE, row.names=FALSE, sep=":", quote=FALSE)

