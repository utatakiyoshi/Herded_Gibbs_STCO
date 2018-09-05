wget http://www.cis.temple.edu/~latecki/TestData/mpeg7shapeB.tar.gz 
tar -zxvf mpeg7shapeB.tar.gz 
R --vanilla --slave < gen.r
g++ gen_input.cpp -o gen_input.out
./gen_input.out -filelist file_list.txt -p 0.3
tail -n +2 file_list.txt | sort -R > input_file_shuffled.txt