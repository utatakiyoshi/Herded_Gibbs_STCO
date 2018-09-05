g++ -O3 image_denoising.cpp -o image_denoising.out
num_expe=5;
tail -n 1 observe_times.txt > result1.csv;
echo "" >> result1.csv
head -n 1 observe_times.txt > out_settings.txt;
K=$[1+${num_expe}*2]; 
echo ${K} >> out_settings.txt;
echo -n "time " >> out_settings.txt;
echo "0.42364893019" > input_file_list.txt
echo "1" >> input_file_list.txt
echo "10" >> input_file_list.txt
echo "100" >> input_file_list.txt
cat input_file_shuffled.txt >> input_file_list.txt
./image_denoising.out -expe_name Gibbs -Gibbs
./image_denoising.out -expe_name HG -HG
./image_denoising.out -expe_name HG_shared -sHG
./image_denoising.out -expe_name HG_1Bin -exsHG
./image_denoising.out -expe_name LFSR -qmcmc 
echo "" >> out_settings.txt

