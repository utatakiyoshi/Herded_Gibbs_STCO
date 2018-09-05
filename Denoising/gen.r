library(magick)
library(foreach)

filenames = list.files("original/", pattern=".gif")

file.remove("file_list.txt");
write.table(as.matrix(c(length(filenames))),"file_list.txt", append=T, quote=F,col.names=F,row.names=F)
for(filename in filenames){
	data <- image_read(paste("original/",filename,sep=""));
	data2 = as.numeric(as.factor(as.raster(data)))
	H = image_info(data)$height;
	W = image_info(data)$width;
	dim(data2)=c(W,H);
	data2=t(data2);
	data2=data2-1
	answer_file_name = paste("inputs/", sub(".gif","_answer.txt",filename), sep="")
	input_file_name = paste("inputs/", sub(".gif",".txt",filename), sep="")
	write.table(t(as.matrix(c(H,W,input_file_name,answer_file_name))),"file_list.txt", append=T, quote=F,col.names=F,row.names=F)
	write.table(data2,answer_file_name,quote=F,col.names=F,row.names=F)
}