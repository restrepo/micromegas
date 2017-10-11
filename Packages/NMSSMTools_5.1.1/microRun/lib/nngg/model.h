*     LanHEP output produced at Thu Oct 23 23:48:20 2014
*     Model named 'NMSSM(../spect.slha)+hgg'

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

      double precision EE, MZ, MW, MW2, GG, alfSMZ, MbMb, McMc, Mtp
      double precision Ml, Mq, CW, SW, C2W, LamQCD, rd, mu, Lambda
      double precision Kappa, aLambda, aKappa, muP, At, Ab, Al
      double precision Ml2, Ml3, Mr2, Mr3, Mq2, Mq3, Mu2, Mu3
      double precision Md2, Md3, tB, Mh1, Mh2, Mh3, Mha, Mhb, MHc
      double precision MNE1, MNE2, MNE3, MNE4, MNE5, MC1, MC2
      double precision MSG, MSuL, MSuR, MSdL, MSdR, MScL, MScR
      double precision MSsL, MSsR, MSt1, MSt2, MSb1, MSb2, QSUSY
      double precision MSeL, MSeR, MSmL, MSmR, MSl1, MSl2, MSne
      double precision MSnm, MSnl, Zl11, Zl12, Zl21, Zl22, Zh11
      double precision Zh12, Zh13, Zh21, Zh22, Zh23, Zh31, Zh32
      double precision Zh33, Za11, Za12, Za13, Za21, Za22, Za23
      double precision Zn11, Zn12, Zn13, Zn14, Zn15, Zn21, Zn22
      double precision Zn23, Zn24, Zn25, Zn31, Zn32, Zn33, Zn34
      double precision Zn35, Zn41, Zn42, Zn43, Zn44, Zn45, Zn51
      double precision Zn52, Zn53, Zn54, Zn55, Zu11, Zu12, Zu21
      double precision Zu22, Zv11, Zv12, Zv21, Zv22, Zt11, Zt12
      double precision Zt21, Zt22, Zb11, Zb12, Zb21, Zb22, la1
      double precision la2, la3, la4, la5, la6, la7, la1s, la2s
      double precision la3s, la4s, la5s, la6s, la7s, la8s, aa1
      double precision aa2, aa3, aa4, aa5, aa6, mB1, mB2, X, dMb
      double precision xif, xis, MM3, MSP, sb, cb, t2b, xvev, Pa12
      double precision Pa22, Pa11, Pa21, Td3, Q, Mt, Mb, Mc, xH2
      double precision yH2, dMd, Td2, Au, Ad, fiuu, fidd, ficc
      double precision fiss, Zuu11, Zuu12, Zuu21, Zuu22, Zdd11
      double precision Zdd12, Zdd21, Zdd22, Zcc11, Zcc12, Zcc21
      double precision Zcc22, Zss11, Zss12, Zss21, Zss22, Mele
      double precision Mmu, MbMM, MtMM, NMM11, NMM12, NMM13, NMM14
      double precision NMM15, NMM22, NMM23, NMM24, NMM25, NMM33
      double precision NMM34, NMM35, NMM44, NMM45, NMM55, MG1I
      double precision MG2I

      double precision AAABR(2062)

      double precision quuMass, qudMass, lpdMass, neuMass, chaMass
      double precision sluMass, sldMass, sleMass, squMass, sqvMass
      double precision sqdMass, sqeMass, hisMass, hiaMass, MTR001
      double precision hisW, hiaW
      double precision MTR002, MTR003, MTR004, MTR005, MTR006
      double precision MTR007, MTR008, MTR009, MTR010, MTR011
      double precision MTR012, MTR013, MTR014, MTR015, MTR016
      double precision MTR017, MTR018, MTR019, MTR020, MTR021
      double precision MTR022, MTR023, MTR024, MTR025, MTR026
      double precision MTR027, MTR028, MTR029, MTR030, MTR031
      double precision MTR032, MTR033, MTR034, MTR035, MTR036
      double precision MTR037, MTR038, MTR039, MTR040, MTR041
      double precision MTR042, MTR043, MTR044, MTR045, MTR046
      double precision MTR047, MTR048, MTR049, MTR050, MTR051
      double precision MTR052, MTR053, MTR054, MTR055, MTR056
      double precision MTR057, MTR058, MTR059, MTR060, MTR061
      double precision MTR062, MTR063, MTR064, MTR065, MTR066
      double precision MTR067, MTR068, MTR069, MTR070, MTR071
      double precision MTR072, MTR073, MTR074, MTR075, MTR076
      double precision MTR077, MTR078, MTR079, MTR080, MTR081
      double precision MTR082, MTR083, MTR084, MTR085, MTR086
      double precision MTR087, MTR088, MTR089, MTR090, MTR091
      double precision MTR092, MTR093, MTR094, MTR095, MTR096
      double precision MTR097, MTR098, MTR099, MTR100, MTR101
      double precision MTR102, MTR103, MTR104, MTR105, MTR106
      double precision MTR107, MTR108, MTR109, MTR110, MTR111
      double precision MTR112, MTR113, MTR114, MTR115, MTR116
      double precision MTR117, MTR118, MTR119, MTR120, MTR121
      double precision MTR122, MTR123, MTR124, MTR125, MTR126
      double precision MTR127, MTR128, MTR129, MTR130, MTR131
      double precision MTR132, MTR133, MTR134, MTR135, MTR136
      double precision MTR137, MTR138, MTR139, MTR140, MTR141
      double precision MTR142, MTR143, MTR144, MTR145, MTR146
      double precision MTR147, MTR148, MTR149, MTR150, MTR151
      double precision MTR152, MTR153, MTR154, MTR155, MTR156
      double precision MTR157, MTR158, MTR159, MTR160, MTR161
      double precision MTR162, MTR163, MTR164, MTR165, MTR166
      double precision MTR167, MTR168, MTR169, MTR170, MTR171
      double precision MTR172, MTR173, MTR174, MTR175, MTR176
      double precision MTR177, MTR178, MTR179, MTR180, MTR181
      double precision MTR182, MTR183, MTR184, MTR185, MTR186
      double precision MTR187, MTR188, MTR189, MTR190, MTR191
      double precision MTR192, MTR193, MTR194, MTR195, MTR196
      double precision MTR197, MTR198, MTR199, MTR200, MTR201
      double precision MTR202, MTR203, MTR204, MTR205, MTR206
      double precision MTR207, MTR208, MTR209, MTR210, MTR211
      double precision MTR212, MTR213, MTR214, MTR215, MTR216
      double precision MTR217, MTR218, MTR219, MTR220, MTR221
      double precision MTR222, MTR223, MTR224, MTR225, MTR226
      double precision MTR227, MTR228, MTR229, MTR230, MTR231
      double precision MTR232, MTR233, MTR234, MTR235, MTR236
      double precision MTR237, MTR238, MTR239, MTR240, MTR241
      double precision MTR242, MTR243, MTR244, MTR245, MTR246
      double precision MTR247, MTR248, MTR249, MTR250, MTR251
      double precision MTR252, MTR253, MTR254, MTR255, MTR256
      double precision MTR257, MTR258, MTR259, MTR260, MTR261
      double precision MTR262, MTR263, MTR264, MTR265, MTR266
      double precision MTR267, MTR268, MTR269, MTR270, MTR271
      double precision MTR272, MTR273, MTR274, MTR275, MTR276
      double precision MTR277, MTR278, MTR279, MTR280, MTR281
      double precision MTR282, MTR283, MTR284, MTR285, MTR286
      double precision MTR287, MTR288, MTR289, MTR290, MTR291
      double precision MTR292, MTR293, MTR294, MTR295, MTR296
      double precision MTR297, MTR298, MTR299, MTR300, MTR301
      double precision MTR302, MTR303, MTR304, MTR305, MTR306
      double precision MTR307, MTR308, MTR309, MTR310, MTR311
      double precision MTR312, MTR313, MTR314, MTR315, MTR316
      double precision MTR317, MTR318, MTR319, MTR320, MTR321
      double precision MTR322, MTR323, MTR324, MTR325, MTR326
      double precision MTR327, MTR328, MTR329, MTR330, MTR331
      double precision MTR332, MTR333, MTR334, MTR335, MTR336
      double precision MTR337, MTR338, MTR339, MTR340, MTR341
      double precision MTR342, MTR343, MTR344, MTR345, MTR346
      double precision MTR347, MTR348, MTR349, MTR350, MTR351
      double precision MTR352, MTR353, MTR354, MTR355, MTR356
      double precision MTR357, MTR358, MTR359, MTR360, MTR361
      double precision MTR362, MTR363, MTR364, MTR365, MTR366
      double precision MTR367, MTR368, MTR369, MTR370, MTR371
      double precision MTR372, MTR373, MTR374, MTR375, MTR376
      double precision MTR377, MTR378, MTR379, MTR380, MTR381
      double precision MTR382, MTR383, MTR384, MTR385, MTR386
      double precision MTR387, MTR388, MTR389, MTR390, MTR391
      double precision MTR392, MTR393, MTR394, MTR395, MTR396
      double precision MTR397, MTR398, MTR399, MTR400, MTR401


      dimension
     & quuMass(3), qudMass(3), lpdMass(3), neuMass(5), chaMass(2), 
     & sluMass(3), sldMass(3), sleMass(3), squMass(3), sqvMass(3), 
     & sqdMass(3), sqeMass(3), hisMass(3), hiaMass(2), MTR001(3),
     & hisW(3), hiaW(2),
     & MTR002(3), MTR003(3), MTR004(3), MTR005(3), MTR006(3), 
     & MTR007(3), MTR008(3), MTR009(3), MTR010(2), MTR011(3), 
     & MTR012(3), MTR013(2,3), MTR014(3,3), MTR015(3,3), 
     & MTR016(3,3), MTR017(3,3), MTR018(3), MTR019(3), MTR020(3), 
     & MTR021(2,3), MTR022(3,3), MTR023(3,3), MTR024(3), 
     & MTR025(3), MTR026(3,3), MTR027(3), MTR028(2,3), MTR029(3,3), 
     & MTR030(3,3), MTR031(3,3), MTR032(3), MTR033(2,3), 
     & MTR034(2,2,3), MTR035(3,3,3), MTR036(2), MTR037(3), 
     & MTR038(3), MTR039(3), MTR040(3), MTR041(3), MTR042(3), 
     & MTR043(3), MTR044(3), MTR045(3), MTR046(3), MTR047(3), 
     & MTR048(3), MTR049(3), MTR050(3), MTR051(3), MTR052(3), 
     & MTR053(2,3), MTR054(3), MTR055(2,2), MTR056(2,2), 
     & MTR057(2,2,2), MTR058(2,2,2), MTR059(2,2,3), MTR060(2,2,3), 
     & MTR061(2,3), MTR062(2,3), MTR063(2,5), MTR064(2,5), 
     & MTR065(2,5), MTR066(2,5), MTR067(2,3), MTR068(2,3), 
     & MTR069(2,3), MTR070(2,3), MTR071(3,2), MTR072(3,3), 
     & MTR073(5,3), MTR074(5,3), MTR075(5,3), MTR076(5,3), 
     & MTR077(2,3), MTR078(2,3), MTR079(5,3), MTR080(5,3), 
     & MTR081(5,3), MTR082(5,3), MTR083(5,3), MTR084(3), 
     & MTR085(3,2), MTR086(3,3), MTR087(3), MTR088(3), MTR089(3), 
     & MTR090(3), MTR091(2,3), MTR092(2,3), MTR093(2,3), 
     & MTR094(2,3), MTR095(3), MTR096(3), MTR097(3), MTR098(3), 
     & MTR099(5,3), MTR100(5,3), MTR101(5,3), MTR102(5,3), 
     & MTR103(3), MTR104(3,2), MTR105(3,3), MTR106(3), MTR107(3), 
     & MTR108(3), MTR109(3), MTR110(5,5), MTR111(5,5,2), 
     & MTR112(5,5,3), MTR113(2,2), MTR114(2,2), MTR115(2,5), 
     & MTR116(2,5), MTR117(5,5), MTR118(3), MTR119(3), MTR120(3), 
     & MTR121(3), MTR122(3), MTR123(3), MTR124(3), MTR125(3), 
     & MTR126(3), MTR127(2), MTR128(2,2), MTR129(3,3), MTR130(3), 
     & MTR131(3), MTR132(3), MTR133(3), MTR134(2,3), MTR135(2,3), 
     & MTR136(3,3), MTR137(3,3), MTR138(3), MTR139(3), MTR140(3), 
     & MTR141(3), MTR142(3), MTR143(3), MTR144(3), MTR145(2,3), 
     & MTR146(2,3), MTR147(3,3), MTR148(3,3), MTR149(3), 
     & MTR150(3), MTR151(3), MTR152(2,3), MTR153(2,3), MTR154(3,3), 
     & MTR155(3,3), MTR156(2), MTR157(3), MTR158(2,2), MTR159(2,3), 
     & MTR160(3,3), MTR161(3,3), MTR162(3,3), MTR163(3,3), 
     & MTR164(3,3), MTR165(3,3), MTR166(3,3), MTR167(3,3), 
     & MTR168(3,3), MTR169(3,3), MTR170(3,3), MTR171(3,3), 
     & MTR172(3,3), MTR173(3,3), MTR174(3,3), MTR175(3,3), 
     & MTR176(3,3), MTR177(3,3), MTR178(3,3), MTR179(3,3), 
     & MTR180(3,3), MTR181(3,3), MTR182(3,3), MTR183(3,3), 
     & MTR184(3,3), MTR185(3,3), MTR186(3,3), MTR187(3,3), 
     & MTR188(3,3), MTR189(3,3), MTR190(3), MTR191(3), MTR192(2,3), 
     & MTR193(3,3), MTR194(3), MTR195(2,3), MTR196(2,3), 
     & MTR197(3,3), MTR198(2,2,3), MTR199(2,2,3), MTR200(2,3,3), 
     & MTR201(3,3,3), MTR202(3,3,3), MTR203(3,3), MTR204(3,3), 
     & MTR205(3,3), MTR206(3,3), MTR207(3,3), MTR208(3,3), 
     & MTR209(3,3), MTR210(3,3), MTR211(3,3), MTR212(3,3), 
     & MTR213(3,3), MTR214(3,3), MTR215(3,3), MTR216(3,3), 
     & MTR217(3,3), MTR218(3), MTR219(2,3), MTR220(3,3), 
     & MTR221(3), MTR222(2,3), MTR223(2,2,3), MTR224(3,3,3), 
     & MTR225(3,3), MTR226(3,3), MTR227(3,3), MTR228(3,3), 
     & MTR229(3,3), MTR230(3,3), MTR231(3,3), MTR232(3,3), 
     & MTR233(3), MTR234(2,3), MTR235(2,2,3), MTR236(3,3,3), 
     & MTR237(3,3), MTR238(3,3), MTR239(3,3), MTR240(3,3), 
     & MTR241(3,3), MTR242(3,3), MTR243(3,3), MTR244(3,3), 
     & MTR245(3,3), MTR246(3,3), MTR247(3,3), MTR248(3,3), 
     & MTR249(3,3), MTR250(3,3), MTR251(3,3), MTR252(3,3), 
     & MTR253(3,3), MTR254(3,3), MTR255(3,3), MTR256(3,3), 
     & MTR257(3,3), MTR258(3,3), MTR259(3,3), MTR260(3,3), 
     & MTR261(3,3), MTR262(3,3), MTR263(3,3), MTR264(3,3), 
     & MTR265(3), MTR266(3), MTR267(3), MTR268(3), MTR269(2,3), 
     & MTR270(2,3), MTR271(3,3), MTR272(3,3), MTR273(3), 
     & MTR274(3), MTR275(2,3), MTR276(2,3), MTR277(3,3), 
     & MTR278(2,2,3), MTR279(2,2,3), MTR280(2,3,3), MTR281(3,3,3), 
     & MTR282(3,3,3), MTR283(3,3), MTR284(3,3), MTR285(3,3), 
     & MTR286(3,3), MTR287(3,3), MTR288(3,3), MTR289(3,3), 
     & MTR290(3,3), MTR291(3,3), MTR292(3,3), MTR293(3,3), 
     & MTR294(3,3), MTR295(3,3), MTR296(3), MTR297(3), MTR298(3), 
     & MTR299(2,3), MTR300(2,3), MTR301(3,3), MTR302(3,3), 
     & MTR303(3), MTR304(2,3), MTR305(2,2,3), MTR306(3,3,3), 
     & MTR307(3,3), MTR308(3,3), MTR309(3,3), MTR310(3,3), 
     & MTR311(3,3), MTR312(3,3), MTR313(3,3), MTR314(3,3), 
     & MTR315(3,3), MTR316(3,3), MTR317(3,3), MTR318(3,3), 
     & MTR319(3,3), MTR320(3), MTR321(3), MTR322(3), MTR323(3), 
     & MTR324(2,3), MTR325(2,3), MTR326(3,3), MTR327(2,2,3), 
     & MTR328(2,2,3), MTR329(2,3,3), MTR330(3,3,3), MTR331(3,3,3), 
     & MTR332(3,3), MTR333(3,3), MTR334(3), MTR335(3), MTR336(2,3), 
     & MTR337(2,2,3), MTR338(3,3,3), MTR339(2), MTR340(2,2), 
     & MTR341(3,3), MTR342(2,2), MTR343(3,3), MTR344(2,2,2), 
     & MTR345(2,3,3), MTR346(2,2,2,2), MTR347(2,2,3,3), 
     & MTR348(3,3,3,3), MTR349(2), MTR350(2), MTR351(3), 
     & MTR352(3), MTR353(3), MTR354(3), MTR355(3), MTR356(3), 
     & MTR357(3), MTR358(3), MTR359(3), MTR360(3), MTR361(3), 
     & MTR362(3), MTR363(3), MTR364(3), MTR365(3), MTR366(3), 
     & MTR367(3), MTR368(3), MTR369(3), MTR370(3), MTR371(3), 
     & MTR372(3), MTR373(3), MTR374(3), MTR375(3), MTR376(3), 
     & MTR377(3), MTR378(3), MTR379(3), MTR380(3), MTR381(3), 
     & MTR382(3), MTR383(3), MTR384(3), MTR385(3), MTR386(3), 
     & MTR387(3), MTR388(3), MTR389(3), MTR390(3), MTR391(3), 
     & MTR392(3), MTR393(3), MTR394(3), MTR395(3), MTR396(3), 
     & MTR397(3), MTR398(2,2), MTR399(2,2), MTR400(3,3), 
     & MTR401(3,3)

      common /mdl_mtrces/
     & quuMass, qudMass, lpdMass, neuMass, chaMass, sluMass, 
     & sldMass, sleMass, squMass, sqvMass, sqdMass, sqeMass, 
     & hisMass, hiaMass, MTR001, MTR002, MTR003, MTR004, 
     & MTR005, MTR006, MTR007, MTR008, MTR009, MTR010, MTR011, 
     & MTR012, MTR013, MTR014, MTR015, MTR016, MTR017, MTR018, 
     & MTR019, MTR020, MTR021, MTR022, MTR023, MTR024, MTR025, 
     & MTR026, MTR027, MTR028, MTR029, MTR030, MTR031, MTR032, 
     & MTR033, MTR034, MTR035, MTR036, MTR037, MTR038, MTR039, 
     & MTR040, MTR041, MTR042, MTR043, MTR044, MTR045, MTR046, 
     & MTR047, MTR048, MTR049, MTR050, MTR051, MTR052, MTR053, 
     & MTR054, MTR055, MTR056, MTR057, MTR058, MTR059, MTR060, 
     & MTR061, MTR062, MTR063, MTR064, MTR065, MTR066, MTR067, 
     & MTR068, MTR069, MTR070, MTR071, MTR072, MTR073, MTR074, 
     & MTR075, MTR076, MTR077, MTR078, MTR079, MTR080, MTR081, 
     & MTR082, MTR083, MTR084, MTR085, MTR086, MTR087, MTR088, 
     & MTR089, MTR090, MTR091, MTR092, MTR093, MTR094, MTR095, 
     & MTR096, MTR097, MTR098, MTR099, MTR100, MTR101, MTR102, 
     & MTR103, MTR104, MTR105, MTR106, MTR107, MTR108, MTR109, 
     & MTR110, MTR111, MTR112, MTR113, MTR114, MTR115, MTR116, 
     & MTR117, MTR118, MTR119, MTR120, MTR121, MTR122, MTR123, 
     & MTR124, MTR125, MTR126, MTR127, MTR128, MTR129, MTR130, 
     & MTR131, MTR132, MTR133, MTR134, MTR135, MTR136, MTR137, 
     & MTR138, MTR139, MTR140, MTR141, MTR142, MTR143, MTR144, 
     & MTR145, MTR146, MTR147, MTR148, MTR149, MTR150, MTR151, 
     & MTR152, MTR153, MTR154, MTR155, MTR156, MTR157, MTR158, 
     & MTR159, MTR160, MTR161, MTR162, MTR163, MTR164, MTR165, 
     & MTR166, MTR167, MTR168, MTR169, MTR170, MTR171, MTR172, 
     & MTR173, MTR174, MTR175, MTR176, MTR177, MTR178, MTR179, 
     & MTR180, MTR181, MTR182, MTR183, MTR184, MTR185, MTR186, 
     & MTR187, MTR188, MTR189, MTR190, MTR191, MTR192, MTR193, 
     & MTR194, MTR195, MTR196, MTR197, MTR198, MTR199, MTR200, 
     & MTR201, MTR202, MTR203, MTR204, MTR205, MTR206, MTR207, 
     & MTR208, MTR209, MTR210, MTR211, MTR212, MTR213, MTR214, 
     & MTR215, MTR216, MTR217, MTR218, MTR219, MTR220, MTR221, 
     & MTR222, MTR223, MTR224, MTR225, MTR226, MTR227, MTR228, 
     & MTR229, MTR230, MTR231, MTR232, MTR233, MTR234, MTR235, 
     & MTR236, MTR237, MTR238, MTR239, MTR240, MTR241, MTR242, 
     & MTR243, MTR244, MTR245, MTR246, MTR247, MTR248, MTR249, 
     & MTR250, MTR251, MTR252, MTR253, MTR254, MTR255, MTR256, 
     & MTR257, MTR258, MTR259, MTR260, MTR261, MTR262, MTR263, 
     & MTR264, MTR265, MTR266, MTR267, MTR268, MTR269, MTR270, 
     & MTR271, MTR272, MTR273, MTR274, MTR275, MTR276, MTR277, 
     & MTR278, MTR279, MTR280, MTR281, MTR282, MTR283, MTR284, 
     & MTR285, MTR286, MTR287, MTR288, MTR289, MTR290, MTR291, 
     & MTR292, MTR293, MTR294, MTR295, MTR296, MTR297, MTR298, 
     & MTR299, MTR300, MTR301, MTR302, MTR303, MTR304, MTR305, 
     & MTR306, MTR307, MTR308, MTR309, MTR310, MTR311, MTR312, 
     & MTR313, MTR314, MTR315, MTR316, MTR317, MTR318, MTR319, 
     & MTR320, MTR321, MTR322, MTR323, MTR324, MTR325, MTR326, 
     & MTR327, MTR328, MTR329, MTR330, MTR331, MTR332, MTR333, 
     & MTR334, MTR335, MTR336, MTR337, MTR338, MTR339, MTR340, 
     & MTR341, MTR342, MTR343, MTR344, MTR345, MTR346, MTR347, 
     & MTR348, MTR349, MTR350, MTR351, MTR352, MTR353, MTR354, 
     & MTR355, MTR356, MTR357, MTR358, MTR359, MTR360, MTR361, 
     & MTR362, MTR363, MTR364, MTR365, MTR366, MTR367, MTR368, 
     & MTR369, MTR370, MTR371, MTR372, MTR373, MTR374, MTR375, 
     & MTR376, MTR377, MTR378, MTR379, MTR380, MTR381, MTR382, 
     & MTR383, MTR384, MTR385, MTR386, MTR387, MTR388, MTR389, 
     & MTR390, MTR391, MTR392, MTR393, MTR394, MTR395, MTR396, 
     & MTR397, MTR398, MTR399, MTR400, MTR401

      common /mdl_para/
     &    EE, MZ, MW, MW2, GG, alfSMZ, MbMb, McMc, Mtp, Ml, Mq, CW,
     &    SW, C2W, LamQCD, rd, mu, Lambda, Kappa, aLambda, aKappa,
     &    muP, At, Ab, Al, Ml2, Ml3, Mr2, Mr3, Mq2, Mq3, Mu2, Mu3,
     &    Md2, Md3, tB, Mh1, Mh2, Mh3, Mha, Mhb, MHc, MNE1, MNE2,
     &    MNE3, MNE4, MNE5, MC1, MC2, MSG, MSuL, MSuR, MSdL, MSdR,
     &    MScL, MScR, MSsL, MSsR, MSt1, MSt2, MSb1, MSb2, QSUSY,
     &    MSeL, MSeR, MSmL, MSmR, MSl1, MSl2, MSne, MSnm, MSnl,
     &    Zl11, Zl12, Zl21, Zl22, Zh11, Zh12, Zh13, Zh21, Zh22,
     &    Zh23, Zh31, Zh32, Zh33, Za11, Za12, Za13, Za21, Za22,
     &    Za23, Zn11, Zn12, Zn13, Zn14, Zn15, Zn21, Zn22, Zn23,
     &    Zn24, Zn25, Zn31, Zn32, Zn33, Zn34, Zn35, Zn41, Zn42,
     &    Zn43, Zn44, Zn45, Zn51, Zn52, Zn53, Zn54, Zn55, Zu11,
     &    Zu12, Zu21, Zu22, Zv11, Zv12, Zv21, Zv22, Zt11, Zt12,
     &    Zt21, Zt22, Zb11, Zb12, Zb21, Zb22, la1, la2, la3, la4,
     &    la5, la6, la7, la1s, la2s, la3s, la4s, la5s, la6s, la7s,
     &    la8s, aa1, aa2, aa3, aa4, aa5, aa6, mB1, mB2, X, dMb, xif,
     &    xis, MM3, MSP, sb, cb, t2b, xvev, Pa12, Pa22, Pa11, Pa21,
     &    Td3, Q, Mt, Mb, Mc, xH2, yH2, dMd, Td2, Au, Ad, fiuu,
     &    fidd, ficc, fiss, Zuu11, Zuu12, Zuu21, Zuu22, Zdd11,
     &    Zdd12, Zdd21, Zdd22, Zcc11, Zcc12, Zcc21, Zcc22, Zss11,
     &    Zss12, Zss21, Zss22, Mele, Mmu, MbMM, MtMM, NMM11, NMM12,
     &    NMM13, NMM14, NMM15, NMM22, NMM23, NMM24, NMM25, NMM33,
     &    NMM34, NMM35, NMM44, NMM45, NMM55, MG1I, MG2I, AAABR

