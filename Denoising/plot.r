library(dplyr)
library(ggplot2)
library(foreach)
library(grid)
library(gridExtra)

Nheader=2;
#settings<-read.table("out_settings.txt", nrows=1,stringsAsFactors=F);
#out_filename=(settings[1,1]);
settings<-read.table("out_settings.txt", skip=0,nrows=Nheader);
Tmax=(settings[1,1])
K=(settings[2,1])


settings<-read.table("out_settings.txt", skip=Nheader, nrows=1, colClasses=rep("character",K));
edge_names=settings[1,]
data<-as.data.frame(t(read.table("result1.csv",sep=","))[1:Tmax,]);
colnames(data)=settings;

data2 = data %>%
 tidyr::gather(key, value, -time) %>% 
 tidyr::separate(key, c("algorithm_ab","name"), sep='\\$') %>% 
 tidyr::separate(algorithm_ab, c("algorithm","algorithmb"), sep='-', remove=FALSE) %>% 
 tidyr::spread(key=name, value=value);

data2$"algorithmb"[is.na(data2$"algorithmb")] <- "mean";

data2 <- data2 %>% transform(algorithm= factor(algorithm, levels = c("Gibbs","LFSR","HG","HG_shared","HG_1Bin")))

errors <- aes(ymax = mean + std, ymin = mean - std)
gp = ggplot(data2, aes(x=time,y=mean,group=algorithm,colour=algorithm));
gp <- gp + geom_errorbar(errors, width = 0.5) + geom_path() 
gp <- gp + scale_x_continuous(breaks=c(1,3,7,15,23,31)) + xlab("T") + ylab("reconstruction error");
gp <- gp + theme_bw(base_family = "Helvetica") + theme(plot.margin=margin(0.2,0.2,0.1,0.1, "in"), axis.title=element_text(size=8),axis.text=element_text(size=8,colour="Black"),axis.text.y=element_text(angle = 90,hjust=0.5), panel.grid=element_line(linetype=0), panel.border=element_rect(size=0.8),legend.title=element_blank(), legend.key.size=unit(0.15, "in"), legend.justification=c(1,1), legend.position=c(0.995,0.995), legend.background = element_rect(fill="white"), legend.box.background = element_rect(fill = NULL, size=0.8), legend.box.margin = margin(.2,.2,.2,.2))
#gp <- gp + ylim(c(0,0.01));
ggsave(file = "resultMPEG.eps", plot = gp, width = 3.3, height=3.3)
