!*************************************************************************
! This file is part of
!
!       HiggsBounds 1.2.0
!
! by Philip Bechtle, Oliver Brein, Sven Heinemyer, Georg Weiglein
!    and Karina E. Williams.
!
!  Journal Reference: e-Print: arXiv:0811.4169 [hep-ph], submitted to CPC.
!  Web site: http://www.ippp.dur.ac.uk/HiggsBounds
!
!10/09/2009
!*************************************************************************

        integer tlist_tev(TEV_NSCHN*NHMAX+TEV_NDCHN*NHMAX**2)
        common /key_tev/ tlist_tev

        double precision Mhmin_overall,Mhmax_overall
        common /overall_Mh_range/ Mhmin_overall,Mhmax_overall

        double precision CDF_ZH_llbb_pub_2_pred(1:11)
        double precision CDF_ZH_llbb_pub_2_obs(1:11)
        double precision Mhmin_CDF_ZH_llbb_pub_2,Mhmax_CDF_ZH_llbb_pub_2,sep_CDF_ZH_llbb_pub_2

        common /table_CDF_ZH_llbb_pub_2/ CDF_ZH_llbb_pub_2_pred,CDF_ZH_llbb_pub_2_obs
        common /table_CDF_ZH_llbb_pub_2_info/ Mhmin_CDF_ZH_llbb_pub_2,Mhmax_CDF_ZH_llbb_pub_2,
     &        sep_CDF_ZH_llbb_pub_2        
        

        double precision CDF_ZH_llbb_4_pred(1:11)
        double precision CDF_ZH_llbb_4_obs(1:11)
        double precision Mhmin_CDF_ZH_llbb_4,Mhmax_CDF_ZH_llbb_4,sep_CDF_ZH_llbb_4

        common /table_CDF_ZH_llbb_4/ CDF_ZH_llbb_4_pred,CDF_ZH_llbb_4_obs
        common /table_CDF_ZH_llbb_4_info/ Mhmin_CDF_ZH_llbb_4,Mhmax_CDF_ZH_llbb_4,
     &        sep_CDF_ZH_llbb_4
        
       
        double precision CDF_VH_bbmET_4_pred(1:11)
        double precision CDF_VH_bbmET_4_obs(1:11)
        double precision Mhmin_CDF_VH_bbmET_4,Mhmax_CDF_VH_bbmET_4,sep_CDF_VH_bbmET_4

        common /table_CDF_VH_bbmET_4/ CDF_VH_bbmET_4_pred,CDF_VH_bbmET_4_obs
        common /table_CDF_VH_bbmET_4_info/ Mhmin_CDF_VH_bbmET_4,Mhmax_CDF_VH_bbmET_4,
     &        sep_CDF_VH_bbmET_4            

        double precision CDF_VH_bbmET_pub_pred(1:7)
        double precision CDF_VH_bbmET_pub_obs(1:7)
        double precision Mhmin_CDF_VH_bbmET_pub,Mhmax_CDF_VH_bbmET_pub,sep_CDF_VH_bbmET_pub

        common /table_CDF_VH_bbmET_pub/ CDF_VH_bbmET_pub_pred,CDF_VH_bbmET_pub_obs
        common /table_CDF_VH_bbmET_pub_info/ Mhmin_CDF_VH_bbmET_pub,Mhmax_CDF_VH_bbmET_pub,
     &        sep_CDF_VH_bbmET_pub            
     
       
        
        double precision D0_ZH_llbb_3_pred(1:11)
        double precision D0_ZH_llbb_3_obs(1:11)
        double precision Mhmin_D0_ZH_llbb_3,Mhmax_D0_ZH_llbb_3,sep_D0_ZH_llbb_3

        common /table_D0_ZH_llbb_3/ D0_ZH_llbb_3_pred,D0_ZH_llbb_3_obs
        common /table_D0_ZH_llbb_3_info/ Mhmin_D0_ZH_llbb_3,Mhmax_D0_ZH_llbb_3,
     &        sep_D0_ZH_llbb_3


        double precision D0_VH_bbmET_pred(1:5)
        double precision D0_VH_bbmET_obs(1:5)
        double precision Mhmin_D0_VH_bbmET,Mhmax_D0_VH_bbmET,sep_D0_VH_bbmET

        common /table_D0_VH_bbmET/ D0_VH_bbmET_pred,D0_VH_bbmET_obs
        common /table_D0_VH_bbmET_info/ Mhmin_D0_VH_bbmET,Mhmax_D0_VH_bbmET,
     &        sep_D0_VH_bbmET
     

        double precision D0_VH_bbmET_pub_pred(1:4)
        double precision D0_VH_bbmET_pub_obs(1:4)
        double precision Mhmin_D0_VH_bbmET_pub,Mhmax_D0_VH_bbmET_pub,sep_D0_VH_bbmET_pub

        common /table_D0_VH_bbmET_pub/ D0_VH_bbmET_pub_pred,D0_VH_bbmET_pub_obs
        common /table_D0_VH_bbmET_pub_info/ Mhmin_D0_VH_bbmET_pub,Mhmax_D0_VH_bbmET_pub,
     &        sep_D0_VH_bbmET_pub


        double precision CDF_D0_SM_pred(1:19)
        double precision CDF_D0_SM_obs(1:19)
        double precision Mhmin_CDF_D0_SM,Mhmax_CDF_D0_SM,sep_CDF_D0_SM

        common /table_CDF_D0_SM/ CDF_D0_SM_pred,CDF_D0_SM_obs
        common /table_CDF_D0_SM_info/ Mhmin_CDF_D0_SM,Mhmax_CDF_D0_SM,
     &        sep_CDF_D0_SM
        
        double precision CDF_bH_bb_pred(1:13)
        double precision CDF_bH_bb_obs(1:13)
        double precision Mhmin_CDF_bH_bb,Mhmax_CDF_bH_bb,sep_CDF_bH_bb

        common /table_CDF_bH_bb/ CDF_bH_bb_pred,CDF_bH_bb_obs
        common /table_CDF_bH_bb_info/ Mhmin_CDF_bH_bb,Mhmax_CDF_bH_bb,
     &        sep_CDF_bH_bb

        double precision D0_WH_lnubb_3_pred(1:11)
        double precision D0_WH_lnubb_3_obs(1:11)
        double precision Mhmin_D0_WH_lnubb_3,Mhmax_D0_WH_lnubb_3,sep_D0_WH_lnubb_3

        common /table_D0_WH_lnubb_3/ D0_WH_lnubb_3_pred,D0_WH_lnubb_3_obs
        common /table_D0_WH_lnubb_3_info/ Mhmin_D0_WH_lnubb_3,Mhmax_D0_WH_lnubb_3,
     &        sep_D0_WH_lnubb_3

        double precision CDF_WH_lnubb_5_pred(1:11)
        double precision CDF_WH_lnubb_5_obs(1:11)
        double precision Mhmin_CDF_WH_lnubb_5,Mhmax_CDF_WH_lnubb_5,sep_CDF_WH_lnubb_5

        common /table_CDF_WH_lnubb_5/ CDF_WH_lnubb_5_pred,CDF_WH_lnubb_5_obs
        common /table_CDF_WH_lnubb_5_info/ Mhmin_CDF_WH_lnubb_5,Mhmax_CDF_WH_lnubb_5,
     &        sep_CDF_WH_lnubb_5
        


        double precision D0_bH_bb_tb_limit_2_pred(1:14)
        double precision D0_bH_bb_tb_limit_2_obs(1:14)
        double precision Mhmin_D0_bH_bb_tb_limit_2,Mhmax_D0_bH_bb_tb_limit_2,sep_D0_bH_bb_tb_limit_2

        common /table_D0_bH_bb_tb_limit_2/ D0_bH_bb_tb_limit_2_pred,D0_bH_bb_tb_limit_2_obs
        common /table_D0_bH_bb_tb_limit_2_info/ Mhmin_D0_bH_bb_tb_limit_2,Mhmax_D0_bH_bb_tb_limit_2,
     &        sep_D0_bH_bb_tb_limit_2
        
        double precision D0_bH_bb_2_pred(1:14)
        double precision D0_bH_bb_2_obs(1:14)
        double precision Mhmin_D0_bH_bb_2,Mhmax_D0_bH_bb_2,sep_D0_bH_bb_2

        common /table_D0_bH_bb_2/ D0_bH_bb_2_pred,D0_bH_bb_2_obs
        common /table_D0_bH_bb_2_info/ Mhmin_D0_bH_bb_2,Mhmax_D0_bH_bb_2,
     &        sep_D0_bH_bb_2
        
        double precision D0_ppH_tautau_pub_pred(1:22)
        double precision D0_ppH_tautau_pub_obs(1:22)
        double precision Mhmin_D0_ppH_tautau_pub,Mhmax_D0_ppH_tautau_pub,sep_D0_ppH_tautau_pub

        common /table_D0_ppH_tautau_pub/ D0_ppH_tautau_pub_pred,D0_ppH_tautau_pub_obs
        common /table_D0_ppH_tautau_pub_info/ Mhmin_D0_ppH_tautau_pub,Mhmax_D0_ppH_tautau_pub,
     &        sep_D0_ppH_tautau_pub
        
        double precision CDF_ppH_tautau_pub_pred(1:17)
        double precision CDF_ppH_tautau_pub_obs(1:17)
        double precision Mhmin_CDF_ppH_tautau_pub,Mhmax_CDF_ppH_tautau_pub,sep_CDF_ppH_tautau_pub

        common /table_CDF_ppH_tautau_pub/ CDF_ppH_tautau_pub_pred,CDF_ppH_tautau_pub_obs
        common /table_CDF_ppH_tautau_pub_info/ Mhmin_CDF_ppH_tautau_pub,Mhmax_CDF_ppH_tautau_pub,
     &        sep_CDF_ppH_tautau_pub
        
        double precision D0_WH_WWW_2_pred(1:5)
        double precision D0_WH_WWW_2_obs(1:5)
        double precision Mhmin_D0_WH_WWW_2,Mhmax_D0_WH_WWW_2,sep_D0_WH_WWW_2

        common /table_D0_WH_WWW_2/ D0_WH_WWW_2_pred,D0_WH_WWW_2_obs
        common /table_D0_WH_WWW_2_info/ Mhmin_D0_WH_WWW_2,Mhmax_D0_WH_WWW_2,
     &        sep_D0_WH_WWW_2

        double precision CDF_WH_WWW_2_pred(1:10)
        double precision CDF_WH_WWW_2_obs(1:10)
        double precision Mhmin_CDF_WH_WWW_2,Mhmax_CDF_WH_WWW_2,sep_CDF_WH_WWW_2

        common /table_CDF_WH_WWW_2/ CDF_WH_WWW_2_pred,CDF_WH_WWW_2_obs
        common /table_CDF_WH_WWW_2_info/ Mhmin_CDF_WH_WWW_2,Mhmax_CDF_WH_WWW_2,
     &        sep_CDF_WH_WWW_2

        double precision D0_ppH_gamgam_4_pred(1:21)
        double precision D0_ppH_gamgam_4_obs(1:21)
        double precision Mhmin_D0_ppH_gamgam_4,Mhmax_D0_ppH_gamgam_4,sep_D0_ppH_gamgam_4

        common /table_D0_ppH_gamgam_4/ D0_ppH_gamgam_4_pred,D0_ppH_gamgam_4_obs
        common /table_D0_ppH_gamgam_4_info/ Mhmin_D0_ppH_gamgam_4,Mhmax_D0_ppH_gamgam_4,
     &        sep_D0_ppH_gamgam_4

        double precision D0_ppH_gamgam_pub_pred(1:9)
        double precision D0_ppH_gamgam_pub_obs(1:9)
        double precision Mhmin_D0_ppH_gamgam_pub,Mhmax_D0_ppH_gamgam_pub,sep_D0_ppH_gamgam_pub

        common /table_D0_ppH_gamgam_pub/ D0_ppH_gamgam_pub_pred,D0_ppH_gamgam_pub_obs
        common /table_D0_ppH_gamgam_pub_info/ Mhmin_D0_ppH_gamgam_pub,Mhmax_D0_ppH_gamgam_pub,
     &        sep_D0_ppH_gamgam_pub
     
        double precision D0_ppH_gamgam_pub_2_pred(1:11)
        double precision D0_ppH_gamgam_pub_2_obs(1:11)
        double precision Mhmin_D0_ppH_gamgam_pub_2,Mhmax_D0_ppH_gamgam_pub_2,sep_D0_ppH_gamgam_pub_2

        common /table_D0_ppH_gamgam_pub_2/ D0_ppH_gamgam_pub_2_pred,D0_ppH_gamgam_pub_2_obs
        common /table_D0_ppH_gamgam_pub_2_info/ Mhmin_D0_ppH_gamgam_pub_2,Mhmax_D0_ppH_gamgam_pub_2,
     &        sep_D0_ppH_gamgam_pub_2
     
     

        
        double precision CDF_SM_tautau_pred(1:9)
        double precision CDF_SM_tautau_obs(1:9)
        double precision Mhmin_CDF_SM_tautau,Mhmax_CDF_SM_tautau,sep_CDF_SM_tautau

        common /table_CDF_SM_tautau/ CDF_SM_tautau_pred,CDF_SM_tautau_obs
        common /table_CDF_SM_tautau_info/ Mhmin_CDF_SM_tautau,Mhmax_CDF_SM_tautau,
     &        sep_CDF_SM_tautau

        double precision CDF_D0_SM_2_pred(1:19)
        double precision CDF_D0_SM_2_obs(1:19)
        double precision Mhmin_CDF_D0_SM_2,Mhmax_CDF_D0_SM_2,sep_CDF_D0_SM_2

        common /table_CDF_D0_SM_2/ CDF_D0_SM_2_pred,CDF_D0_SM_2_obs
        common /table_CDF_D0_SM_2_info/ Mhmin_CDF_D0_SM_2,Mhmax_CDF_D0_SM_2,
     &        sep_CDF_D0_SM_2

        double precision D0_ppH_WW_ll_2_pred(1:18)
        double precision D0_ppH_WW_ll_2_obs(1:18)
        double precision Mhmin_D0_ppH_WW_ll_2,Mhmax_D0_ppH_WW_ll_2,sep_D0_ppH_WW_ll_2

        common /table_D0_ppH_WW_ll_2/ D0_ppH_WW_ll_2_pred,D0_ppH_WW_ll_2_obs
        common /table_D0_ppH_WW_ll_2_info/ Mhmin_D0_ppH_WW_ll_2,Mhmax_D0_ppH_WW_ll_2,
     &        sep_D0_ppH_WW_ll_2
     
        double precision CDF_ppH_WW_4_pred(1:10)
        double precision CDF_ppH_WW_4_obs(1:10)
        double precision Mhmin_CDF_ppH_WW_4,Mhmax_CDF_ppH_WW_4,sep_CDF_ppH_WW_4

        common /table_CDF_ppH_WW_4/ CDF_ppH_WW_4_pred,CDF_ppH_WW_4_obs
        common /table_CDF_ppH_WW_4_info/ Mhmin_CDF_ppH_WW_4,Mhmax_CDF_ppH_WW_4,
     &        sep_CDF_ppH_WW_4     
     
     
        double precision CDF_D0_SM_3_pred(1:10)
        double precision CDF_D0_SM_3_obs(1:10)
        double precision Mhmin_CDF_D0_SM_3,Mhmax_CDF_D0_SM_3,sep_CDF_D0_SM_3

        common /table_CDF_D0_SM_3/ CDF_D0_SM_3_pred,CDF_D0_SM_3_obs
        common /table_CDF_D0_SM_3_info/ Mhmin_CDF_D0_SM_3,Mhmax_CDF_D0_SM_3,
     &        sep_CDF_D0_SM_3

        double precision CDF_D0_SM_4_pred(1:10)
        double precision CDF_D0_SM_4_obs(1:10)
        double precision Mhmin_CDF_D0_SM_4,Mhmax_CDF_D0_SM_4,sep_CDF_D0_SM_4

        common /table_CDF_D0_SM_4/ CDF_D0_SM_4_pred,CDF_D0_SM_4_obs
        common /table_CDF_D0_SM_4_info/ Mhmin_CDF_D0_SM_4,Mhmax_CDF_D0_SM_4,
     &        sep_CDF_D0_SM_4


        double precision CDF_SM_pred(1:21)
        double precision CDF_SM_obs(1:21)
        double precision Mhmin_CDF_SM,Mhmax_CDF_SM,sep_CDF_SM

        common /table_CDF_SM/ CDF_SM_pred,CDF_SM_obs
        common /table_CDF_SM_info/ Mhmin_CDF_SM,Mhmax_CDF_SM,
     &        sep_CDF_SM

        double precision CDF_SM_2_pred(1:21)
        double precision CDF_SM_2_obs(1:21)
        double precision Mhmin_CDF_SM_2,Mhmax_CDF_SM_2,sep_CDF_SM_2

        common /table_CDF_SM_2/ CDF_SM_2_pred,CDF_SM_2_obs
        common /table_CDF_SM_2_info/ Mhmin_CDF_SM_2,Mhmax_CDF_SM_2,
     &        sep_CDF_SM_2


        double precision D0_SM_pub_pred(1:21)
        double precision D0_SM_pub_obs(1:21)
        double precision Mhmin_D0_SM_pub,Mhmax_D0_SM_pub,sep_D0_SM_pub

        common /table_D0_SM_pub/ D0_SM_pub_pred,D0_SM_pub_obs
        common /table_D0_SM_pub_info/ Mhmin_D0_SM_pub,Mhmax_D0_SM_pub,
     &        sep_D0_SM_pub


        double precision CDF_ppH_WW_6_pred(1:19)
        double precision CDF_ppH_WW_6_obs(1:19)
        double precision Mhmin_CDF_ppH_WW_6,Mhmax_CDF_ppH_WW_6,sep_CDF_ppH_WW_6

        common /table_CDF_ppH_WW_6/ CDF_ppH_WW_6_pred,CDF_ppH_WW_6_obs
        common /table_CDF_ppH_WW_6_info/ Mhmin_CDF_ppH_WW_6,Mhmax_CDF_ppH_WW_6,
     &        sep_CDF_ppH_WW_6     


        double precision D0_ppH_WW_ll_3_pred(1:18)
        double precision D0_ppH_WW_ll_3_obs(1:18)
        double precision Mhmin_D0_ppH_WW_ll_3,Mhmax_D0_ppH_WW_ll_3,sep_D0_ppH_WW_ll_3

        common /table_D0_ppH_WW_ll_3/ D0_ppH_WW_ll_3_pred,D0_ppH_WW_ll_3_obs
        common /table_D0_ppH_WW_ll_3_info/ Mhmin_D0_ppH_WW_ll_3,Mhmax_D0_ppH_WW_ll_3,
     &        sep_D0_ppH_WW_ll_3

        double precision D0_WH_lnubb_pub_pred(1:11)
        double precision D0_WH_lnubb_pub_obs(1:11)
        double precision Mhmin_D0_WH_lnubb_pub,Mhmax_D0_WH_lnubb_pub,sep_D0_WH_lnubb_pub

        common /table_D0_WH_lnubb_pub/ D0_WH_lnubb_pub_pred,D0_WH_lnubb_pub_obs
        common /table_D0_WH_lnubb_pub_info/ Mhmin_D0_WH_lnubb_pub,Mhmax_D0_WH_lnubb_pub,
     &        sep_D0_WH_lnubb_pub

        double precision CDF_WH_lnubb_pub_2_pred(1:11)
        double precision CDF_WH_lnubb_pub_2_obs(1:11)
        double precision Mhmin_CDF_WH_lnubb_pub_2,Mhmax_CDF_WH_lnubb_pub_2,sep_CDF_WH_lnubb_pub_2

        common /table_CDF_WH_lnubb_pub_2/ CDF_WH_lnubb_pub_2_pred,CDF_WH_lnubb_pub_2_obs
        common /table_CDF_WH_lnubb_pub_2_info/ Mhmin_CDF_WH_lnubb_pub_2,Mhmax_CDF_WH_lnubb_pub_2,
     &        sep_CDF_WH_lnubb_pub_2


        double precision D0_bH_tautau_pub_pred(1:7)
        double precision D0_bH_tautau_pub_obs(1:7)
        double precision Mhmin_D0_bH_tautau_pub,Mhmax_D0_bH_tautau_pub,sep_D0_bH_tautau_pub

        common /table_D0_bH_tautau_pub/ D0_bH_tautau_pub_pred,D0_bH_tautau_pub_obs
        common /table_D0_bH_tautau_pub_info/ Mhmin_D0_bH_tautau_pub,Mhmax_D0_bH_tautau_pub,
     &        sep_D0_bH_tautau_pub

        double precision D0_bH_tautau_pred(1:24)
        double precision D0_bH_tautau_obs(1:24)
        double precision Mhmin_D0_bH_tautau,Mhmax_D0_bH_tautau,sep_D0_bH_tautau

        common /table_D0_bH_tautau/ D0_bH_tautau_pred,D0_bH_tautau_obs
        common /table_D0_bH_tautau_info/ Mhmin_D0_bH_tautau,Mhmax_D0_bH_tautau,
     &        sep_D0_bH_tautau

        double precision D0_ppH_tautau_pred(1:22)
        double precision D0_ppH_tautau_obs(1:22)
        double precision Mhmin_D0_ppH_tautau,Mhmax_D0_ppH_tautau,sep_D0_ppH_tautau

        common /table_D0_ppH_tautau/ D0_ppH_tautau_pred,D0_ppH_tautau_obs
        common /table_D0_ppH_tautau_info/ Mhmin_D0_ppH_tautau,Mhmax_D0_ppH_tautau,
     &        sep_D0_ppH_tautau

        double precision CDF_D0_ppH_tautau_pred(1:11)
        double precision CDF_D0_ppH_tautau_obs(1:11)
        double precision Mhmin_CDF_D0_ppH_tautau,Mhmax_CDF_D0_ppH_tautau,sep_CDF_D0_ppH_tautau

        common /table_CDF_D0_ppH_tautau/ CDF_D0_ppH_tautau_pred,CDF_D0_ppH_tautau_obs
        common /table_CDF_D0_ppH_tautau_info/ Mhmin_CDF_D0_ppH_tautau,Mhmax_CDF_D0_ppH_tautau,
     &        sep_CDF_D0_ppH_tautau
     
