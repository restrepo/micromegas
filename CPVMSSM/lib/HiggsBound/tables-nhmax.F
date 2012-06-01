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

!>>
c#define VERBOSE
c#define DEBUG
!<<


!maximum number of Higgs bosons
#define NHMAX 10

!maximum number of additional data columns in 'additional.dat'
#define MAX_ADDITIONAL_COLUMNS 20

#define TEV_NSCHN 36
#define TEV_NDCHN 0
#define delta_Mh_TEV 10d0

#define LEP_NSCHN 4
#define LEP_NDCHN 9
#define delta_Mh_LEP 2d0

! Value of variable eps in subroutine check_SM_likeness
! Suggested default value: 0.02
#define EPS_SM_LIKENESS 0.02d0

! LEP table identifiers: (first LEP_NSCHN numbers have to correspond to single Higgs analyses)
!			 (double Higgs analyses follow starting with number LEP_NSCHN+1)
#define LEPhZ_bbZ	1
#define LEPhZ_tautauZ	2
#define LEPhZ_anyZ	3
#define LEPhZ_gamgamZ	4
c#define LEPhZ_jjZ        
#define LEPh2Z_h1h1Z_4bZ	LEP_NSCHN+1
#define LEPh2Z_h1h1Z_4tauZ 	LEP_NSCHN+2
#define LEPh2h1_4b		LEP_NSCHN+3	
#define LEPh2h1_4tau		LEP_NSCHN+4
#define LEPh2h1_h1h1h1_6b	LEP_NSCHN+5
#define LEPh2h1_h1h1h1_6tau	LEP_NSCHN+6
#define LEPh2Z_h1h1Z_2b2tau	LEP_NSCHN+7
#define LEPh2h1_2b2tau		LEP_NSCHN+8
#define LEPh2h1_2tau2b		LEP_NSCHN+9




! Tevatron table identifiers:
#define CDF0908_3534  	1
#define CDF9889 	2
#define CDF9891 	3
#define D05876		4
#define D05586		5
#define CDF8961_D05536 	6
#define CDF9284		7
#define D05972		8
#define CDFwhlnubb_090814	9
#define D05726		10
#define D00805_3556	11
#define D00805_2491	12
#define CDF0906_1014	13
#define D05873		14
#define CDF7307v3	15
#define D05858		16
#define D00803_1514	101
#define CDF9248		17
#define CDF9290_D05645	18
#define D05757		19
#define CDF0809_3930 	20
#define CDF9465_D05754	21

c new additons (impl. starting 17.3.2009)
#define CDF9887		22
#define D05871		23
#define D00808_1970     24
#define CDF0906_5613	25
#define D00811_0024	26
#define D05740		27
#define D00808_1266	28
#define CDF0802_0432	29
#define D00901_1887	30
#define CDF9674		31
#define D00712_0598	32
#define CDF9888_D05980	33
#define CDF9713_D05889  34
#define D05985  	35
#define CDF9897         36

!If ONE_SELECTED is defined with an integer value between 1 and TEV_NSCHN, 
!only the Tevatron analysis with the selected number is considered for exclusion.
*#define ONE_SELECTED 36
