library(colorspace)
k_colors <- rainbow_hcl(4, c=150, l=50, start=70, end=340)[c(4,1,2,3)]
pchs <- c(19,1,2,3,4,5,6)
ltys <- c(1,5,2,3)

argv <- commandArgs(TRUE)
stat_flags <- rep(FALSE,6);
dash_flag = 0;
dashprefix = "";
i = 1;
fins = c();
alg_names = c();
while (i <= length(argv)){
	if(argv[i] == "-epsprefix"){ i=i+1; epsprefix=argv[i]; }
	if(argv[i] == "-dash"){ i=i+1; dash_flag=TRUE; dashprefix=argv[i];  }
	if(argv[i] == "-fin"){ i=i+1; fins=c(fins, argv[i]); }
	if(argv[i] == "-alg"){ i=i+1; alg_names=c(alg_names, argv[i]); }
	if(argv[i] == "-All"){ stat_flags = rep(TRUE,6); }
	if(argv[i] == "-overall"){ stat_flags[1] = TRUE; }
	if(argv[i] == "-prev"){ stat_flags[2] = TRUE; }
	if(argv[i] == "-herding"){ stat_flags[3] = TRUE; }
	if(argv[i] == "-approx"){ stat_flags[4] = TRUE; }
	if(argv[i] == "-cor"){ stat_flags[5] = TRUE; }
	if(argv[i] == "-bin"){ stat_flags[6] = TRUE; }
	i=i+1;
}

indexs_single = 1:6;
labels_single = c("all", "z", "herding", "approx", "cor", "bin"); labels_single=paste("D", labels_single,sep="");
indexs_data_single = 2:7;

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
	indexs_data = c(indexs_data, indexs_data_single + 7*(i-1))
	labels = c(labels, paste(alg_names[i], labels_single,sep=" "));
	if(i==1){
		data <- read.table(fins[i]);
	} else {
		data <- cbind(data, read.table(fins[i]));
	}
}

print(indexs_data);
#print(data[,indexs_data]);


epspostfix0 = "D";
if(dash_flag){
epspostfix0 = "D_dash";
}
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
matplot(data[,indexs_data] ,type="o", ylab="", xlab="log_2(T)", log="y", col=k_colors[indexs_alg], pch=pchs[indexs_stat], lty=ltys[indexs_alg], cex=0.6, lwd=0.7)

if(dash_flag){
print(dashprefix);
dashdata<-read.table(dashprefix);
print(dashdata);
dashpt=dashdata[1];
abline(h = dashpt, lty=2);
}

par(xpd=NA);
legend("bottomleft", legend=labels[nonzeros], col=k_colors[indexs_alg[nonzeros]], pch=pchs[indexs_stat[nonzeros]], lty=ltys[indexs_alg[nonzeros]], y.intersp=0.6, pt.cex=0.6);
dev.off()

