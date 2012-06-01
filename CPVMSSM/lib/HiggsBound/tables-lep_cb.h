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

        integer tlist_lep(LEP_NSCHN*NHMAX+LEP_NDCHN*NHMAX**2)
        common /key_lep/ tlist_lep


        double precision S95_14b_pred(1:217)
        double precision S95_14b_obs(1:217)
        double precision S95_14c_pred(1:233)
        double precision S95_14c_obs(1:233)
        double precision S95_15_pred(0:119,0:60)
        double precision S95_15_obs(0:119,0:60)
        double precision S95_16_pred(0:119,0:60)
        double precision S95_16_obs(0:119,0:60)
        double precision S95_18_pred(0:180,0:180)
        double precision S95_18_obs(0:180,0:180)
        double precision S95_19_pred(0:180,0:180)
        double precision S95_19_obs(0:180,0:180)
        double precision S95_20_pred(0:179,0:90)
        double precision S95_20_obs(0:179,0:90)
        double precision S95_21_pred(0:179,0:90)
        double precision S95_21_obs(0:179,0:90)
        double precision S95_22_pred(0:119,0:60)
        double precision S95_22_obs(0:119,0:60)
        double precision S95_23_pred(0:180,0:180)
        double precision S95_23_obs(0:180,0:180)
        double precision S95_24_pred(0:180,0:180)
        double precision S95_24_obs(0:180,0:180)
        
        double precision Mhmin14,Mhmax14,sep14 
        double precision Mhmin14b,Mhmax14b,sep14b 
        double precision Mhmin14c,Mhmax14c,sep14c 
        double precision Mh1min15,Mh1max15,sep1_15,Mh2min15,Mh2max15,sep2_15
        double precision Mh1min18,Mh1max18,sep1_18,Mh2min18,Mh2max18,sep2_18
        double precision Mh1min20,Mh1max20,sep1_20,Mh2min20,Mh2max20,sep2_20

        common /table14b/ S95_14b_pred,S95_14b_obs
        common /table14c/ S95_14c_pred,S95_14c_obs
        common /table14info/ Mhmin14,Mhmax14,sep14 
        common /table14info/ Mhmin14b,Mhmax14b,sep14b 
        common /table14info/ Mhmin14c,Mhmax14c,sep14c 
        common /table15/ S95_15_pred,S95_15_obs
        common /table16/ S95_16_pred,S95_16_obs
        common /table22/ S95_22_pred,S95_22_obs
        common /table15info/ Mh1min15,Mh1max15,sep1_15,Mh2min15,Mh2max15,sep2_15
        common /table18/ S95_18_pred,S95_18_obs
        common /table18info/ Mh1min18,Mh1max18,sep1_18,Mh2min18,Mh2max18,sep2_18
        common /table19/ S95_19_pred,S95_19_obs
        common /table20/ S95_20_pred,S95_20_obs
        common /table20info/ Mh1min20,Mh1max20,sep1_20,Mh2min20,Mh2max20,sep2_20
        common /table21/ S95_21_pred,S95_21_obs
        
        common /table23/ S95_23_pred,S95_23_obs
        common /table24/ S95_24_pred,S95_24_obs

! new additions
        double precision S95_hZ_anyZ_pred(1:100),S95_hZ_anyZ_obs(1:100)
        double precision Mhmin_hZ_anyZ,Mhmax_hZ_anyZ,sep_hZ_anyZ 
        common /table_hZ_anyZ/ S95_hZ_anyZ_pred,S95_hZ_anyZ_obs
        common /table_hZ_anyZ_info/ Mhmin_hZ_anyZ,Mhmax_hZ_anyZ,sep_hZ_anyZ
         
        double precision S95_hZ_gamgamZ_pred(1:97),S95_hZ_gamgamZ_obs(1:97)
        double precision Mhmin_hZ_gamgamZ,Mhmax_hZ_gamgamZ,sep_hZ_gamgamZ 
        common /table_hZ_gamgamZ/ S95_hZ_gamgamZ_pred,S95_hZ_gamgamZ_obs
        common /table_hZ_gamgamZ_info/ Mhmin_hZ_gamgamZ,Mhmax_hZ_gamgamZ,sep_hZ_gamgamZ
       