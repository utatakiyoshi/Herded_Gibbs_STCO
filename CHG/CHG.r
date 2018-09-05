library(colorspace)
y_legend = c("Gibbs", "Herding(2)", "Herding(4)", "Herding(6)", "Herding(8)");
k_colors <- rainbow_hcl(length(y_legend), c=200, l=40);
pchs <- (1: length(y_legend));

data<-read.table("out.txt"); tmp<-data[,1]
num_measure=22;
head=1;
x<-c();
for(i in 1:length(y_legend)) {
    x<-cbind(x, tmp[head:(head+num_measure-1)]);
    head<-head+num_measure;
}
measureT <- tmp[head:(head+num_measure-1)];
epsname = paste("CHG.eps",sep="")
postscript(file = epsname, width = 3.3, height = 3.3, family = "Helvetica",horizontal = FALSE, onefile = FALSE, paper = "special")
setEPS();
par(mar = c(2.1, 2.1, 1, 1.1))
par(mgp = c(1, 0.3, 0))
par(oma = c(0, 0, 0, 0))
par(ps=8)
x_label = "log_2(T)"
y_label = "median error"
matplot(x,type="o", log="y", lty=1, ylab=y_label, xlab=x_label,col=k_colors, pch=pchs, cex=0.6, lwd=0.7)
par(xpd=NA)
legend("bottomleft", legend=y_legend, pch=pchs, lty=1, col=k_colors, y.intersp=0.6, pt.cex=0.6);
dev.off()
