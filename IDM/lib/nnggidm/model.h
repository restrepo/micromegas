*     LanHEP output produced at Mon Apr 13 18:09:41 2009
*     Model named 'MSSMrc cha1 nans'

      double precision Sqrt2, pi, degree, hbar_c2,bogus
      parameter (Sqrt2=1.41421356237309504880168872421D0)
      parameter (pi = 3.1415926535897932384626433832795029D0)
      parameter (degree = pi/180D0)
      parameter (hbar_c2 = 3.8937966D8)
      parameter (bogus = -1D123)
      double complex cI
      parameter (cI = (0D0, 1D0))

      double precision Divergence
      common /renorm/ Divergence

      double precision EE, GG, alfSMZ, SW, CW, MZ, MW, Q, Mtp
      double precision MbMb, McMc, Lqcd, Mb, Mt, Mc,  Mbp, Mcp
      double precision MHX, MH3, MHC, laL, la2, mu2, la3, la5
      double precision la4, wZ, wW, Mm, Ml, Mqu, Mqd, Ms, Mh, wh, Mqu2
      double precision Mqd2, Mc2, Ms2, Mt2, Mb2, Mm2, Ml2, Mh2
      double precision MH32, MHC2, MHX2, EE2, MW2, MZ2

      double precision AAABR(58)

      common /mdl_para/
     &    EE, GG, alfSMZ, SW, CW, MZ, MW, Q, Mtp, MbMb, McMc, Lqcd,
     &    Mb, Mt, Mc,  Mbp, Mcp, MHX, MH3, MHC, laL, la2, mu2,
     &    la3, la5, la4, wZ, wW, Mm, Ml, Mqu, Mqd, Ms, Mh, wh, Mqu2,
     &    Mqd2, Mc2, Ms2, Mt2, Mb2, Mm2, Ml2, Mh2, MH32, MHC2,
     &    MHX2, EE2, MW2, MZ2, AAABR

