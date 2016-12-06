library (ggplot2)

bed_file <- read.table("final_cov.bed", header = F, sep = '\t')

grepped_file <- grep (as.character("NM_000044:AR"), bed_file, fixed = TRUE)

plot (bed_file$V5, main="BRCA2 Coverage", xlab="Position", ylab="Number of Reads",  border = "gray")
abline(h=15, col = 'red')

ggplot(bed_file, aes(x=V2, y=V5)) +
  geom_point( colour="dodgerblue4", size=2)+
  ylab(substitute("Depth"))+
  xlab(substitute("Position"))+
  theme_bw()
  
