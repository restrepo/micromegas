# SUSY Les Houches Accord 2 - MSSM spectrum + Decays
# SPheno v3.2.3  
# W. Porod, Comput. Phys. Commun. 153 (2003) 275-315, hep-ph/0301101
# in case of problems send email to porod@physik.uni-wuerzburg.de
# Created: 16.03.2016,  10:46
Block SPINFO         # Program information
     1   SPheno      # spectrum calculator
     2   v3.2.3      # version number
#
Block SPhenoINFO     # SPheno specific information
    1      2         # using 2-loop RGEs
    2      1         # using running masses for boundary conditions at mZ
# Either the general MSSM or a model has been used
# which has not yet been implemented in the LesHouches standard
Block MINPAR  # Input parameters
    3    1.38621799E+01  # tanb at m_Z    
    4    1.00000000E+00  # Sign(mu)
Block EXTPAR  # non-universal input parameters
    1    7.31000000E+02  # M_1
    2    7.66000000E+02  # M_2
    3    1.90630000E+03  # M_3
   23    1.28630000E+03  # mu
   25    1.34000000E+01  # tan(beta)
   26    1.59290000E+03  # m_A, pole mass
Block SMINPUTS  # SM parameters
         1     1.27908953E+02  # alpha_em^-1(MZ)^MSbar
         2     1.16637000E-05  # G_mu [GeV^-2]
         3     1.17200000E-01  # alpha_s(MZ)^MSbar
         4     9.11884000E+01  # m_Z(pole)
         5     4.23000000E+00  # m_b(m_b), MSbar
         6     1.73300000E+02  # m_t(pole)
         7     1.77700000E+00  # m_tau(pole)
         8     0.00000000E+00  # m_nu_3
        11     5.10998910E-04  # m_e(pole)
        12     0.00000000E+00  # m_nu_1
        13     1.05658000E-01  # m_muon(pole)
        14     0.00000000E+00  # m_nu_2
        21     5.00000000E-03  # m_d(2 GeV), MSbar
        22     3.00000000E-03  # m_u(2 GeV), MSbar
        23     1.05000000E-01  # m_s(2 GeV), MSbar
        24     1.27000000E+00  # m_c(m_c), MSbar
Block gauge Q=  1.00000000E+03  # (SUSY scale)
   1    3.61188547E-01  # g'(Q)^DRbar
   2    6.35312107E-01  # g(Q)^DRbar
   3    1.03727668E+00  # g3(Q)^DRbar
Block Yu Q=  1.00000000E+03  # (SUSY scale)
  1  1     8.65073699E-06   # Y_u(Q)^DRbar
  2  2     3.66214424E-03   # Y_c(Q)^DRbar
  3  3     8.51065329E-01   # Y_t(Q)^DRbar
Block Yd Q=  1.00000000E+03  # (SUSY scale)
  1  1     1.86547821E-04   # Y_d(Q)^DRbar
  2  2     3.91750111E-03   # Y_s(Q)^DRbar
  3  3     1.87186139E-01   # Y_b(Q)^DRbar
Block Ye Q=  1.00000000E+03  # (SUSY scale)
  1  1     3.98960519E-05   # Y_e(Q)^DRbar
  2  2     8.24920001E-03   # Y_mu(Q)^DRbar
  3  3     1.38694002E-01   # Y_tau(Q)^DRbar
Block Au Q=  1.00000000E+03  # (SUSY scale)
  1  1     1.15597087E-05   # A_u(Q)^DRbar
  2  2     2.73064067E-08   # A_c(Q)^DRbar
  3  3    -3.28094672E+03   # A_t(Q)^DRbar
Block Ad Q=  1.00000000E+03  # (SUSY scale)
  1  1     5.36055579E-07   # A_d(Q)^DRbar
  2  2     2.55264765E-08   # A_s(Q)^DRbar
  3  3     5.34227589E-10   # A_b(Q)^DRbar
Block Ae Q=  1.00000000E+03  # (SUSY scale)
  1  1     2.50651369E-06   # A_e(Q)^DRbar
  2  2     1.21223876E-08   # A_mu(Q)^DRbar
  3  3     7.21011716E-10   # A_tau(Q)^DRbar
Block MSOFT Q=  1.00000000E+03  # soft SUSY breaking masses at Q
   1    7.31000000E+02  # M_1
   2    7.66000000E+02  # M_2
   3    1.90630000E+03  # M_3
  21    8.37502596E+05  # M^2_(H,d)
  22   -1.73717326E+06  # M^2_(H,u)
  31    3.58960001E+03  # M_(L,11)
  32    3.58960001E+03  # M_(L,22)
  33    3.58960001E+03  # M_(L,33)
  34    3.58960001E+03  # M_(E,11)
  35    3.58960001E+03  # M_(E,22)
  36    3.58960001E+03  # M_(E,33)
  41    3.25260001E+03  # M_(Q,11)
  42    3.25260001E+03  # M_(Q,22)
  43    1.63430000E+03  # M_(Q,33)
  44    3.25260001E+03  # M_(U,11)
  45    3.25260001E+03  # M_(U,22)
  46    1.05440000E+03  # M_(U,33)
  47    3.25260001E+03  # M_(D,11)
  48    3.25260001E+03  # M_(D,22)
  49    1.63430000E+03  # M_(D,33)
Block MASS  # Mass spectrum
#   PDG code      mass          particle
         6     1.73300000E+02  # m_t(pole)
        23     9.11884000E+01  # m_Z(pole)
        24     8.02829600E+01  # W+
        15     1.77700000E+00  # m_tau(pole)
        25     1.26264062E+02  # h0
        35     1.59266396E+03  # H0
        36     1.59290000E+03  # A0
        37     1.59510083E+03  # H+
   1000001     3.27161821E+03  # ~d_L
   2000001     3.26583699E+03  # ~d_R
   1000002     3.27086933E+03  # ~u_L
   2000002     3.26621255E+03  # ~u_R
   1000003     3.27162151E+03  # ~s_L
   2000003     3.26583473E+03  # ~s_R
   1000004     3.27086995E+03  # ~c_L
   2000004     3.26621309E+03  # ~c_R
   1000005     1.63100406E+03  # ~b_1
   2000005     1.66821717E+03  # ~b_2
   1000006     1.00903221E+03  # ~t_1
   2000006     1.68511088E+03  # ~t_2
   1000011     3.59585986E+03  # ~e_L-
   2000011     3.59263235E+03  # ~e_R-
   1000012     3.59463332E+03  # ~nu_eL
   1000013     3.59588129E+03  # ~mu_L-
   2000013     3.59261525E+03  # ~mu_R-
   1000014     3.59463478E+03  # ~nu_muL
   1000015     3.59035082E+03  # ~tau_1-
   2000015     3.59936745E+03  # ~tau_2-
   1000016     3.59504530E+03  # ~nu_tauL
   1000021     2.04994572E+03  # ~g
   1000022     7.38116743E+02  # ~chi_10
   1000023     8.02389530E+02  # ~chi_20
   1000025    -1.28843447E+03  # ~chi_30
   1000035     1.29450637E+03  # ~chi_40
   1000024     8.02336007E+02  # ~chi_1+
   1000037     1.29507533E+03  # ~chi_2+
# Higgs mixing
Block alpha # Effective Higgs mixing angle
          -7.53431907E-02   # alpha
Block Hmix Q=  1.00000000E+03  # Higgs mixing parameters
   1    1.28630000E+03  # mu
   2    1.34000000E+01  # tan[beta](Q)
   3    2.43737301E+02  # v(Q)
   4    2.51241342E+06  # m^2_A(Q)
Block stopmix # stop mixing matrix
   1  1     2.87029206E-01   # Re[R_st(1,1)]
   1  2     9.57921831E-01   # Re[R_st(1,2)]
   2  1    -9.57921831E-01   # Re[R_st(2,1)]
   2  2     2.87029206E-01   # Re[R_st(2,2)]
Block sbotmix # sbottom mixing matrix
   1  1     9.42847229E-01   # Re[R_sb(1,1)]
   1  2     3.33225303E-01   # Re[R_sb(1,2)]
   2  1    -3.33225303E-01   # Re[R_sb(2,1)]
   2  2     9.42847229E-01   # Re[R_sb(2,2)]
Block staumix # stau mixing matrix
   1  1     5.88915355E-01   # Re[R_sta(1,1)]
   1  2     8.08194720E-01   # Re[R_sta(1,2)]
   2  1    -8.08194720E-01   # Re[R_sta(2,1)]
   2  2     5.88915355E-01   # Re[R_sta(2,2)]
Block Nmix # neutralino mixing matrix
   1  1    -9.96329795E-01   # Re[N(1,1)]
   1  2     4.90757338E-02   # Re[N(1,2)]
   1  3    -5.94074774E-02   # Re[N(1,3)]
   1  4     3.72728268E-02   # Re[N(1,4)]
   2  1    -5.71952656E-02   # Re[N(2,1)]
   2  2    -9.91252111E-01   # Re[N(2,2)]
   2  3     9.90386268E-02   # Re[N(2,3)]
   2  4    -6.58733917E-02   # Re[N(2,4)]
   3  1    -1.44357834E-02   # Re[N(3,1)]
   3  2     2.43670896E-02   # Re[N(3,2)]
   3  3     7.06133489E-01   # Re[N(3,3)]
   3  4     7.07512083E-01   # Re[N(3,4)]
   4  1    -6.20261992E-02   # Re[N(4,1)]
   4  2     1.20071100E-01   # Re[N(4,2)]
   4  3     6.98596877E-01   # Re[N(4,3)]
   4  4    -7.02636524E-01   # Re[N(4,4)]
Block Umix # chargino mixing matrix
   1  1    -9.89437892E-01   # Re[U(1,1)]
   1  2     1.44957432E-01   # Re[U(1,2)]
   2  1     1.44957432E-01   # Re[U(2,1)]
   2  2     9.89437892E-01   # Re[U(2,2)]
Block Vmix # chargino mixing matrix
   1  1    -9.95315051E-01   # Re[V(1,1)]
   1  2     9.66847981E-02   # Re[V(1,2)]
   2  1     9.66847981E-02   # Re[V(2,1)]
   2  2     9.95315051E-01   # Re[V(2,2)]
Block SPhenoLowEnergy  # low energy observables
    1    2.84219019E-04   # BR(b -> s gamma)
    2    1.59158657E-06   # BR(b -> s mu+ mu-)
    3    3.54492855E-05   # BR(b -> s nu nu)
    4    5.19220889E-11   # BR(Bd -> mu+ mu-)
    5    3.39196504E-09   # BR(Bs -> mu+ mu-)
    6    1.10264264E-04   # BR(B_u -> tau nu)
    7    9.97058784E-01   # BR(B_u -> tau nu)/BR(B_u -> tau nu)_SM
    8    2.73462192E-01   # |Delta(M_Bd)| [ps^-1] 
    9    2.02923711E+01   # |Delta(M_Bs)| [ps^-1] 
   20    1.54035888E-14   # Delta(g-2)_electron/2
   21    3.33599139E-11   # Delta(g-2)_muon/2
   22   -2.08070902E-08   # Delta(g-2)_tau/2
   23    0.00000000E+00   # electric dipole moment of the electron
   24    0.00000000E+00   # electric dipole moment of the muon
   25    0.00000000E+00   # electric dipole moment of the tau
   26    0.00000000E+00   # Br(mu -> e gamma)
   27    0.00000000E+00   # Br(tau -> e gamma)
   28    0.00000000E+00   # Br(tau -> mu gamma)
   29    0.00000000E+00   # Br(mu -> 3 e)
   30    0.00000000E+00   # Br(tau -> 3 e)
   31    0.00000000E+00   # Br(tau -> 3 mu)
   39    1.15383112E-02   # Delta(rho_parameter)
   40    0.00000000E+00   # BR(Z -> e mu)
   41    0.00000000E+00   # BR(Z -> e tau)
   42    0.00000000E+00   # BR(Z -> mu tau)
Block FWCOEF Q=  9.11884000E+01  # Wilson coefficients at scale Q
#    id        order  M        value         comment
Block FWCOEF Q=  1.60000000E+02  # Wilson coefficients at scale Q
#    id        order  M        value         comment
     0305 4422   00   0    -1.87802455E-01   # C7
     0305 4422   00   1     1.88656214E-02   # C7
     0305 4322   00   1     3.84543422E-04   # C7'
     0305 6421   00   0    -9.49037728E-02   # C8
     0305 6421   00   1     1.82243664E-03   # C8
     0305 6321   00   1     2.98415449E-05   # C8'
 03051111 4133   00   0    -4.57850328E-01   # C9 e+e-
 03051111 4133   00   1     3.39121662E-05   # C9 e+e-
 03051111 4233   00   1    -4.00564905E-06   # C9' e+e-
 03051111 4137   00   0    -3.95218028E+00   # C10 e+e-
 03051111 4137   00   1     3.44351321E-04   # C10 e+e-
 03051111 4237   00   1     1.00382708E-04   # C10' e+e-
 03051313 4133   00   0    -4.57850328E-01   # C9 mu+mu-
 03051313 4133   00   1     3.39082987E-05   # C9 mu+mu-
 03051313 4233   00   1    -4.00569127E-06   # C9' mu+mu-
 03051313 4137   00   0    -3.95218028E+00   # C10 mu+mu-
 03051313 4137   00   1     3.44347376E-04   # C10 mu+mu-
 03051313 4237   00   1     1.00382750E-04   # C10' mu+mu-
 03051212 4237   00   0     1.49078738E+00   # C11 nu_1 nu_1
 03051212 4237   00   1    -8.64216215E-05   # C11 nu_1 nu_1
 03051212 4137   00   1    -2.41122281E-05   # C11' nu_1 nu_1
 03051414 4237   00   0     1.49078738E+00   # C11 nu_2 nu_2
 03051414 4237   00   1    -8.64216215E-05   # C11 nu_2 nu_2
 03051414 4137   00   1    -2.41122186E-05   # C11' nu_2 nu_2
 03051616 4237   00   0     1.49078738E+00   # C11 nu_3 nu_3
 03051616 4237   00   1    -8.64216204E-05   # C11 nu_3 nu_3
 03051616 4137   00   1    -2.41095381E-05   # C11' nu_3 nu_3
