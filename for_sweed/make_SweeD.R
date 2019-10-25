library(tidyverse)

inp <- read.table("angsd_mafs_file_pop_chr", header = T)
m <- mutate(inp, x2 = knownEM*nInd*2, n = nInd*2, folded = nInd*0)
input <- mutate(m, x = round(x2, digits = 0))
final <- select(input, position, x, n, folded)
write.table(final, "sweed_input_pop_chr.txt", col.names=T, row.names = F, quote = F, sep = "\t")
