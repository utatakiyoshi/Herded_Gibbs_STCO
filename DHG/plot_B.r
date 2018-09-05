library(colorspace)
k_colors <- rainbow_hcl(6, c=200, l=40)[1:6];

argv <- commandArgs(TRUE)
i = 1;
fins = c();
while (i <= length(argv)){
	if(argv[i] == "-epsprefix"){ i=i+1; epsprefix=argv[i]; }
	if(argv[i] == "-fin"){ i=i+1; fins=c(fins, argv[i]); }
	i=i+1;
}

cor_index = 6;

cor_vec = c();
indexs_stat=c();
indexs_alg=c();
indexs_data=c();
data=data.frame();
for(i in 1:length(fins)){
		data <- read.table(fins[i]);
		cor = data[,cor_index]
		final_cor = cor[length(cor)];
		cor_vec = c(cor_vec, final_cor);
}

# Large fig #
epspostfix0 = "B";
epsname = paste(paste(epsprefix, epspostfix0, sep="_"), ".eps", sep="")
print(epsname);
postscript(file = epsname, width = 3.3, height = 3.3, family = "Helvetica",horizontal = FALSE, onefile = FALSE, paper = "special")
setEPS();
par(mar = c(2.1, 2.1, 1, 1.1))
par(mgp = c(1, 0.3, 0))
par(oma = c(0, 0, 0, 0))
par(ps=8)

plot(log2(4:16), cor_vec, type="o", log="y", ylab="Dcor", xlab="log_2(B)", cex=0.6, lwd=0.7)
dev.off()
