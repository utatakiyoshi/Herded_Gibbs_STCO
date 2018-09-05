#!/bin/bash
# tdiff
r --vanilla --slave --args -fin result_expe_tdiff/out_G1_l.txt \
-epsprefix figs/tdiff_G1_l \
-tmpprefix result_expe_tdiff/tmp_tdiff_G1_l.txt \
< plot_tdiff.r

r --vanilla --slave --args -fin result_expe_tdiff/out_G2_l.txt \
-epsprefix figs/tdiff_G2_l \
-tmpprefix result_expe_tdiff/tmp_tdiff_G2_l.txt \
< plot_tdiff.r

# HG
# [HG_G1_l_D_dash_Dall_Dz_Dherding_Dcor_.eps]
r --vanilla --slave --args -D -overall -cor -prev -herding \
-fin result_expe_disc/out_G1_l_HG.txt -alg '' \
-dash result_expe_tdiff/tmp_tdiff_G1_l.txt \
-epsprefix figs/HG_G1_l < plot_stat.r

# [HG_G2_l_D_dash_Dall_Dz_Dherding_Dcor_.eps]
r --vanilla --slave --args -D -overall -cor -prev -herding \
-fin result_expe_disc/out_G2_l_HG.txt -alg '' \
-dash result_expe_tdiff/tmp_tdiff_G2_l.txt \
-epsprefix figs/HG_G2_l < plot_stat.r

# disc
# [disc_G4_l_B4_D_Dall_Dz_Dherding_Dapprox_Dcor_B=4_B=64.eps]
r --vanilla --slave --args -D -overall -cor -prev -approx -herding \
-fin result_expe_disc/out_G4_l_B4.txt -alg 'B=4' \
-fin result_expe_disc/out_G4_l_B64.txt -alg 'B=64' \
-epsprefix figs/disc_G4_l_B4 < plot_stat.r

# B
# [B_G4_l_B.eps]
r --vanilla --slave --args \
-fin result_expe_disc/out_G4_long_l_B4.txt \
-fin result_expe_disc/out_G4_long_l_B5.txt \
-fin result_expe_disc/out_G4_long_l_B6.txt \
-fin result_expe_disc/out_G4_long_l_B7.txt \
-fin result_expe_disc/out_G4_long_l_B8.txt \
-fin result_expe_disc/out_G4_long_l_B9.txt \
-fin result_expe_disc/out_G4_long_l_B10.txt \
-fin result_expe_disc/out_G4_long_l_B11.txt \
-fin result_expe_disc/out_G4_long_l_B12.txt \
-fin result_expe_disc/out_G4_long_l_B13.txt \
-fin result_expe_disc/out_G4_long_l_B14.txt \
-fin result_expe_disc/out_G4_long_l_B15.txt \
-fin result_expe_disc/out_G4_long_l_B16.txt \
-epsprefix figs/B_G4_l < plot_B.r

# GDR (RDHG)
# [GDR_G4_l_B64_D_Dall_Dapprox_Dbin_Gibbs_DHG_B=64_RDHG_B=64.eps]
r --vanilla --slave --args -D -overall -approx -bin \
-fin result_expe_disc/out_G4_l.txt -alg Gibbs \
-fin result_expe_disc/out_G4_l_B64.txt -alg 'DHG B=64' \
-fin result_expe_disc/out_RDHG_G4_l_B64.txt -alg 'RDHG B=64' \
-epsprefix figs/GDR_G4_l_B64 < plot_stat.r

# GBD (BEG)
# [GBD_G1_l_D_Dherding_Dcor_Gibbs_BEG_c=10_BEG_c=1_HG.eps]
r --vanilla --slave --args -D -herding -cor \
-fin result_expe_disc/out_G1_l.txt -alg Gibbs \
-fin result_expe_disc/out_G1_l_BER10.txt -alg 'BEG c=10' \
-fin result_expe_disc/out_G1_l_BER1.txt -alg 'BEG c=1' \
-fin result_expe_disc/out_G1_l_HG.txt -alg HG \
-epsprefix figs/GBD_G1_l < plot_stat.r

#lin
# [Lin_G3_l___Gibbs_Gibbs+RB_DHG_B=64_RDHG_B=64]
r --vanilla --slave --args \
-fin result_expe_lin/out_G3_l_Gibbs.txt -alg Gibbs \
-fin result_expe_lin/out_G3_l_RB.txt -alg Gibbs+RB \
-fin result_expe_lin/out_G3_l_B64.txt -alg 'DHG B=64' \
-fin result_expe_lin/out_RDHG_G3_l_B64.txt -alg 'RDHG B=64' \
-epsprefix figs/Lin_G3_l < plot_lin.r

# [Lin_G3_s___Gibbs_Gibbs+RB_DHG_B=64_RDHG_B=64]
r --vanilla --slave --args \
-fin result_expe_lin/out_G3_s_Gibbs.txt -alg Gibbs \
-fin result_expe_lin/out_G3_s_RB.txt -alg Gibbs+RB \
-fin result_expe_lin/out_G3_s_B64.txt -alg 'DHG B=64' \
-fin result_expe_lin/out_RDHG_G3_s_B64.txt -alg 'RDHG B=64' \
-epsprefix figs/Lin_G3_s < plot_lin.r

#additional lin
# [Lin_G3_l___Gibbs_Gibbs+RB_DHG_B=64_DHG_B=64+RB_RDHG_B=64_RDHG_B=64+RB.eps]
r --vanilla --slave --args \
-fin result_expe_lin/out_G3_l_Gibbs.txt -alg Gibbs \
-fin result_expe_lin/out_G3_l_RB.txt -alg Gibbs+RB \
-fin result_expe_lin/out_G3_l_B64.txt -alg 'DHG B=64' \
-fin result_expe_lin/out_G3_l_B64_RB.txt -alg 'DHG B=64+RB' \
-fin result_expe_lin/out_RDHG_G3_l_B64.txt -alg 'RDHG B=64' \
-fin result_expe_lin/out_RDHG_G3_l_B64_RB.txt -alg 'RDHG B=64+RB' \
-epsprefix figs/Lin_G3_l < plot_lin.r

# [Lin_G3_s___Gibbs_Gibbs+RB_DHG_B=64_DHG_B=64+RB_RDHG_B=64_RDHG_B=64+RB.eps]
r --vanilla --slave --args \
-fin result_expe_lin/out_G3_s_Gibbs.txt -alg Gibbs \
-fin result_expe_lin/out_G3_s_RB.txt -alg Gibbs+RB \
-fin result_expe_lin/out_G3_s_B64.txt -alg 'DHG B=64' \
-fin result_expe_lin/out_G3_s_B64_RB.txt -alg 'DHG B=64+RB' \
-fin result_expe_lin/out_RDHG_G3_s_B64.txt -alg 'RDHG B=64' \
-fin result_expe_lin/out_RDHG_G3_s_B64_RB.txt -alg 'RDHG B=64+RB' \
-epsprefix figs/Lin_G3_s < plot_lin.r

# [Lin_G3_s___Gibbs_Gibbs+RB_DHG_B=64_RDHG_B=64]
r --vanilla --slave --args \
-fin result_expe_lin/out_G3_s_Gibbs.txt -alg Gibbs \
-fin result_expe_lin/out_G3_s_B4.txt -alg 'DHG B=4' \
-fin result_expe_lin/out_G3_s_B16.txt -alg 'DHG B=16' \
-fin result_expe_lin/out_G3_s_B64.txt -alg 'DHG B=64' \
-epsprefix figs/Lin_G3_s < plot_lin.r


cp figs/HG_G1_l_D_dash_Dall_Dz_Dherding_Dcor_.eps figs/Fig3a.eps
cp figs/HG_G2_l_D_dash_Dall_Dz_Dherding_Dcor_.eps figs/Fig3b.eps
cp figs/disc_G4_l_B4_D_Dall_Dz_Dherding_Dapprox_Dcor_B=4_B=64.eps figs/Fig4.eps
cp figs/B_G4_l_B.eps figs/Fig6.eps
cp figs/GBD_G1_l_D_Dherding_Dcor_Gibbs_BEG_c=10_BEG_c=1_HG.eps figs/Fig7.eps
cp figs/GDR_G4_l_B64_D_Dall_Dapprox_Dbin_Gibbs_DHG_B=64_RDHG_B=64.eps figs/Fig8.eps
cp figs/Lin_G3_l___Gibbs_Gibbs+RB_DHG_B=64_RDHG_B=64.eps figs/Fig9a.eps
cp figs/Lin_G3_s___Gibbs_Gibbs+RB_DHG_B=64_RDHG_B=64.eps figs/Fig9b.eps
cp figs/Lin_G3_l___Gibbs_Gibbs+RB_DHG_B=64_DHG_B=64+RB_RDHG_B=64_RDHG_B=64+RB.eps figs/Fig9a_revised.eps 
cp figs/Lin_G3_s___Gibbs_Gibbs+RB_DHG_B=64_DHG_B=64+RB_RDHG_B=64_RDHG_B=64+RB.eps figs/Fig9b_revised.eps 
cp figs/Lin_G3_s___Gibbs_DHG_B=4_DHG_B=16_DHG_B=64.eps figs/Fig1c.eps
