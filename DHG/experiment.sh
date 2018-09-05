#!/bin/bash
rm result_expe_tdiff/*.txt
rm result_expe_disc/*.txt
rm result_expe_lin/*.txt

g++ main.cpp -O3 -o a.out

# ==== gen ====
./a.out -draw -seed 120122 -Hmax 1 -Hmin 0 -Kmax 1 -Kmin 0 -random -ring > graphs/G1_largec.txt
./a.out -draw -seed 120322 -Hmax 1 -Hmin -1 -Kmax 1 -Kmin -1 -random -ring > graphs/G2_largec.txt
./a.out -draw -seed 120522 -Hmax 0.5 -Hmin 0 -Kmax 0.5 -Kmin 0 -random -fullyconnected > graphs/G3_largec.txt
./a.out -draw -seed 120522 -Hmax 0.1 -Hmin 0 -Kmax 0.1 -Kmin 0 -random -fullyconnected > graphs/G3_smallc.txt
./a.out -draw -seed 120722 -Hmax 0.5 -Hmin -0.5 -Kmax 0.5 -Kmin -0.5 -random -fullyconnected > graphs/G4_largec.txt

# ==== expe ====
# tdiff
./a.out < graphs/G1_largec.txt -algorithm 1 -num_measure 30 -tdiff > result_expe_tdiff/out_G1_l.txt
./a.out < graphs/G2_largec.txt -algorithm 1 -num_measure 32 -tdiff > result_expe_tdiff/out_G2_l.txt

# HG
./a.out < graphs/G1_largec.txt -algorithm 3 > result_expe_disc/out_G1_l_HG.txt
./a.out < graphs/G2_largec.txt -algorithm 3 > result_expe_disc/out_G2_l_HG.txt

# disc
./a.out < graphs/G4_largec.txt -algorithm 3 -DHG_type 1 -range 0 -divide_B_inv 4 > result_expe_disc/out_G4_l_B4.txt
./a.out < graphs/G4_largec.txt -algorithm 3 -DHG_type 1 -range 0 -divide_B_inv 64 > result_expe_disc/out_G4_l_B64.txt

# B
./a.out < graphs/G4_largec.txt -algorithm 3 -DHG_type 1 -range 0 -divide_B_inv 4 -num_measure 24 > result_expe_disc/out_G4_long_l_B4.txt
./a.out < graphs/G4_largec.txt -algorithm 3 -DHG_type 1 -range 0 -divide_B_inv 5 -num_measure 24 > result_expe_disc/out_G4_long_l_B5.txt
./a.out < graphs/G4_largec.txt -algorithm 3 -DHG_type 1 -range 0 -divide_B_inv 6 -num_measure 24 > result_expe_disc/out_G4_long_l_B6.txt
./a.out < graphs/G4_largec.txt -algorithm 3 -DHG_type 1 -range 0 -divide_B_inv 7 -num_measure 24 > result_expe_disc/out_G4_long_l_B7.txt
./a.out < graphs/G4_largec.txt -algorithm 3 -DHG_type 1 -range 0 -divide_B_inv 8 -num_measure 24 > result_expe_disc/out_G4_long_l_B8.txt
./a.out < graphs/G4_largec.txt -algorithm 3 -DHG_type 1 -range 0 -divide_B_inv 9 -num_measure 24 > result_expe_disc/out_G4_long_l_B9.txt
./a.out < graphs/G4_largec.txt -algorithm 3 -DHG_type 1 -range 0 -divide_B_inv 10 -num_measure 24 > result_expe_disc/out_G4_long_l_B10.txt
./a.out < graphs/G4_largec.txt -algorithm 3 -DHG_type 1 -range 0 -divide_B_inv 11 -num_measure 24 > result_expe_disc/out_G4_long_l_B11.txt
./a.out < graphs/G4_largec.txt -algorithm 3 -DHG_type 1 -range 0 -divide_B_inv 12 -num_measure 24 > result_expe_disc/out_G4_long_l_B12.txt
./a.out < graphs/G4_largec.txt -algorithm 3 -DHG_type 1 -range 0 -divide_B_inv 13 -num_measure 26 > result_expe_disc/out_G4_long_l_B13.txt
./a.out < graphs/G4_largec.txt -algorithm 3 -DHG_type 1 -range 0 -divide_B_inv 14 -num_measure 26 > result_expe_disc/out_G4_long_l_B14.txt
./a.out < graphs/G4_largec.txt -algorithm 3 -DHG_type 1 -range 0 -divide_B_inv 15 -num_measure 26 > result_expe_disc/out_G4_long_l_B15.txt
./a.out < graphs/G4_largec.txt -algorithm 3 -DHG_type 1 -range 0 -divide_B_inv 16 -num_measure 26 > result_expe_disc/out_G4_long_l_B16.txt

#GDR (RDHG)
./a.out < graphs/G4_largec.txt -algorithm 1 > result_expe_disc/out_G4_l.txt
# ./a.out < graphs/G4_largec.txt -algorithm 3 -DHG_type 1 -range 0 -divide_B 0.015625 > result_expe_disc/out_G4_l_B64.txt
./a.out < graphs/G4_largec.txt -algorithm 3 -DHG_type 5 -range 0 -divide_B 0.015625 > result_expe_disc/out_RDHG_G4_l_B64.txt

#GBD (BER)
./a.out < graphs/G1_largec.txt -algorithm 4 -num_measure 24 -range 1 > result_expe_disc/out_G1_l_BER1.txt
./a.out < graphs/G1_largec.txt -algorithm 4 -num_measure 24 -range 10 > result_expe_disc/out_G1_l_BER10.txt
./a.out < graphs/G1_largec.txt -algorithm 1 -num_measure 24 > result_expe_disc/out_G1_l.txt

#lin
./a.out < graphs/G3_largec.txt -algorithm 1 -num_trial 100 -lin > result_expe_lin/out_G3_l_Gibbs.txt
./a.out < graphs/G3_largec.txt -algorithm 1 -num_trial 100 -lin -RB > result_expe_lin/out_G3_l_RB.txt
./a.out < graphs/G3_largec.txt -algorithm 3 -num_trial 100 -DHG_type 1 -range 0 -divide_B 0.015625 -lin  > result_expe_lin/out_G3_l_B64.txt
./a.out < graphs/G3_largec.txt -algorithm 3 -num_trial 100 -DHG_type 5 -range 0 -divide_B 0.015625 -lin  > result_expe_lin/out_RDHG_G3_l_B64.txt

./a.out < graphs/G3_smallc.txt -algorithm 1 -num_trial 100 -lin > result_expe_lin/out_G3_s_Gibbs.txt
./a.out < graphs/G3_smallc.txt -algorithm 1 -num_trial 100 -lin -RB > result_expe_lin/out_G3_s_RB.txt
./a.out < graphs/G3_smallc.txt -algorithm 3 -num_trial 100 -DHG_type 1 -range 0 -divide_B 0.015625 -lin  > result_expe_lin/out_G3_s_B64.txt
./a.out < graphs/G3_smallc.txt -algorithm 3 -num_trial 100 -DHG_type 5 -range 0 -divide_B 0.015625 -lin  > result_expe_lin/out_RDHG_G3_s_B64.txt

# additional lin
./a.out < graphs/G3_largec.txt -algorithm 3 -num_trial 100 -DHG_type 1 -range 0 -divide_B 0.015625 -lin -RB > result_expe_lin/out_G3_l_B64_RB.txt
./a.out < graphs/G3_largec.txt -algorithm 3 -num_trial 100 -DHG_type 5 -range 0 -divide_B 0.015625 -lin -RB > result_expe_lin/out_RDHG_G3_l_B64_RB.txt
./a.out < graphs/G3_smallc.txt -algorithm 3 -num_trial 100 -DHG_type 1 -range 0 -divide_B 0.015625 -lin -RB > result_expe_lin/out_G3_s_B64_RB.txt
./a.out < graphs/G3_smallc.txt -algorithm 3 -num_trial 100 -DHG_type 5 -range 0 -divide_B 0.015625 -lin -RB > result_expe_lin/out_RDHG_G3_s_B64_RB.txt

./a.out < graphs/G3_smallc.txt -algorithm 3 -num_trial 100 -DHG_type 1 -range 0 -divide_B 0.25 -lin  > result_expe_lin/out_G3_s_B4.txt
./a.out < graphs/G3_smallc.txt -algorithm 3 -num_trial 100 -DHG_type 1 -range 0 -divide_B 0.125 -lin  > result_expe_lin/out_G3_s_B8.txt
./a.out < graphs/G3_smallc.txt -algorithm 3 -num_trial 100 -DHG_type 1 -range 0 -divide_B 0.0625 -lin  > result_expe_lin/out_G3_s_B16.txt
./a.out < graphs/G3_smallc.txt -algorithm 3 -num_trial 100 -DHG_type 1 -range 0 -divide_B 0.03125 -lin  > result_expe_lin/out_G3_s_B32.txt

