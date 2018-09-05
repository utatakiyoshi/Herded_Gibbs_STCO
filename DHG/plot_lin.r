library(colorspace)
k_colors <- rainbow_hcl(6, c=150, l=50, start=70, end=340)[c(6,1,2,3,4,5)]
ltys <- c(1,5,2,3)

argv <- commandArgs(TRUE)
stat_flags <- rep(TRUE,1);
dash_flag = 0;
dashprefix = "";
i = 1;
fins = c();
alg_names = c();
while (i <= length(argv)){
	if(argv[i] == "-epsprefix"){ i=i+1; epsprefix=argv[i]; }
	if(argv[i] == "-fin"){ i=i+1; fins=c(fins, argv[i]); }
	if(argv[i] == "-alg"){ i=i+1; alg_names=c(alg_names, argv[i]); }
	i=i+1;
}

indexs_single = 1;
labels_single = c(""); labels_single=paste("", labels_single,sep="");
indexs_data_single = 1;

indexs_data_single=indexs_data_single[stat_flags]
indexs_single=indexs_single[stat_flags]
labels_single=labels_single[stat_flags]

indexs_stat=c();
indexs_alg=c();
indexs_data=c();
labels=c();
data=data.frame();
for(i in 1:length(fins)){
	indexs_stat = c(indexs_stat, indexs_single);
	indexs_alg = c(indexs_alg, rep(i, length(indexs_single)));
	indexs_data = c(indexs_data, indexs_data_single + 1*(i-1))
	labels = c(labels, paste(alg_names[i], labels_single,sep=" "));
	print(fins[i]);
	if(i==1){
		data <- read.table(fins[i]);
	} else {
		data <- cbind(data, read.table(fins[i]));
	}
}

print(indexs_data);
#print(data[,indexs_data]);

maxi=max(data[,indexs_data])
mini=min(data[,indexs_data])
medi=sqrt(maxi*mini)
print(maxi)
maxi=10^(ceiling(log10(maxi)))
print(maxi)
mini=maxi/5000000

epspostfix0 = "";
epspostfix1 = paste(labels_single, collapse="_");
epspostfix2 = paste(alg_names, collapse="_");

# Large fig #
epsname = paste(paste(epsprefix,epspostfix0, epspostfix1, epspostfix2,sep="_"), ".eps", sep="")
epsname = gsub(" ","_",epsname)
print(epsname);
postscript(file = epsname, width = 3.3, height = 3.3, family = "Helvetica",horizontal = FALSE, onefile = FALSE, paper = "special")
setEPS();
par(mar = c(2.1, 2.1, 1, 1.1))
par(mgp = c(1, 0.3, 0))
par(oma = c(0, 0, 0, 0))
par(ps=8)
nonzeros = (1:length(indexs_data))[apply(data,2,max)[indexs_data]>0]
matplot(data[,indexs_data] ,type="o", ylab="estimation error", xlab="log_2(T)", log="y", ylim=c(mini,maxi), col=k_colors[indexs_alg], pch=indexs_alg, lty=1, cex=0.6, lwd=0.7)
par(xpd=NA);
legend("bottomleft", legend=labels[nonzeros], col=k_colors[indexs_alg[nonzeros]], pch=indexs_alg[nonzeros], lty=1, y.intersp=0.6, pt.cex=0.6);
dev.off()

