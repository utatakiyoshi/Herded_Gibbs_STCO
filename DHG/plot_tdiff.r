argv <- commandArgs(TRUE)
i=1
while (i <= length(argv)){
	if(argv[i] == "-epsprefix"){ i=i+1; epsprefix=argv[i]; }
	if(argv[i] == "-tmpprefix"){ i=i+1; tmpprefix=argv[i]; }
	if(argv[i] == "-fin"){ i=i+1; fin=argv[i]; }
	i=i+1;
}
N=8;
M=2^(N-3);

cnts=c();
names=c();
for(i in 1:M){
	cnt=0;
	it=i-1;
	name="";
	# it = edcba
	# cnt = abcde
	for(j in 1:(N-3)){
#		cnt = cnt + (it %% 2);
		cnt = cnt * 2 + (it %% 2);
		if(j==1){name=it%%2;}
		else{name = paste(name,",",it%%2,sep="");}
		it = floor(it/2);
	}
	names=c(names, paste("(",name,")",sep=""));
	cnts=c(cnts,cnt);
}
print(names);


indexs=c()
for(i in 0:M){
#indexs[abcde]= edcba + 1
	indexs=c(indexs, (1:M)[cnts==i])
}



index_y<-c(1,3,2,4)

txtname = paste(epsprefix, ".tex", sep="")
tmptxtname = paste(tmpprefix, "", sep="")
labels=c("$(0,0)$","$(0,1)$","$(1,0)$","$(1,1)$" )
out <- file(txtname, "w")
tmpout <- file(tmptxtname, "w")
data<-read.table(fin, nrows=4);
tmpsum=0;
for(i in index_y){
pyt=data[i,2]
t=data[i,3]
tdiff=data[i,4]
tmp = 2*pyt*t*t*tdiff;
tmpsum =tmpsum+tmp;
line <- paste(labels[i],
			formatC(pyt,digits=4, format="f"), # pyt
			formatC(t,digits=4, format="f"), # t
			formatC(tdiff,digits=2, format="e"), # Tdiff
			formatC(tmp,digits=2, format="e"), # Dcor
sep="&")
line <- paste(line, "\\\\", sep="");
writeLines(line,out);
}

line <- paste(
"\\multicolumn{4}{c}{}",
formatC(tmpsum,digits=2, format="e"),
sep="&")
line <- paste("\\cmidrule(){1-5}", line, "\\\\", sep="");
writeLines(line,out);

tmpline <- paste("",
formatC(tmpsum,digits=15, format="e"),
sep=" ")
writeLines(tmpline,tmpout);

close(out);
close(tmpout);
