snp2window = function(x) {
  require(zoo)
  print(paste0("Analysing chromosome ", x))
  # subset snps from a certain chromosome
  sub = subset(inp, CHR == x)
  # if the number of SNPs is equal/greater then the SNP window size
  if (nrow(sub) >= snp.win) {
    # start position
    s = rollapply(data=sub$BP, width=snp.win, by=snp.stp, FUN=min)
    # end position
    e = rollapply(data=sub$BP, width=snp.win, by=snp.stp, FUN=max)
    # calculate mean Pvalue in chosen SNP windows
    f = rollapply(data=sub$Pvalue, width=snp.win, by=snp.stp, FUN=mean)
    # add calculations to table
    c = rep(x, length(s))
    df = data.frame(chr = c, sta = s, end = e, Pvalue = f)
    return(df)
  }
}
