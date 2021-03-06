* vert_sld5.F
* form factors for {MNE1, MNE1} -> {0, 0}
* generated by FormCalc 12 Apr 2013 15:13

	subroutine vert_sld5
	implicit none

#include "vars.h"

        Cloop(1) = Cloop(1) + 
     -   (Cval(cc00,iint14(sld5))*
     -       (1/(4.D0*Pi**2)*
     -          (Abb5*AAABR(1072)*MTR019(sld5)*MTR181(1,1)) + 
     -         1/(4.D0*Pi**2)*
     -          (Abb7*AAABR(1072)*MTR019(sld5)*MTR182(1,1))) + 
     -      Cval(cc12,iint14(sld5))*
     -       (-(1/(4.D0*Pi**2)*
     -            (AbbSum30*AAABR(1072)*MTR019(sld5)*MTR181(1,1)))-
     -           1/(4.D0*Pi**2)*
     -          (AbbSum31*AAABR(1072)*MTR019(sld5)*MTR182(1,1))) + 
     -      1/(16.D0*Pi**2)*(AbbSum68*MNE1*AAABR(1072)*
     -         B1(S,sldMass(sld5)**2,sldMass(sld5)**2)*MTR060(sld5)*
     -         (MTR191(1,1) + MTR192(1,1))) + 
     -      B0(S,sldMass(sld5)**2,sldMass(sld5)**2)*
     -       (-(1/(16.D0*Pi**2)*
     -            (Abb5*AAABR(1072)*MTR019(sld5)*MTR181(1,1))) + 
     -         AAABR(1072)*
     -          (-(1/(16.D0*Pi**2)*
     -               (Abb7*MTR019(sld5)*MTR182(1,1))) + 
     -            AbbSum68*MNE1*
     -             (1/(32.D0*Pi**2)*MTR191(1,1) + 
     -               1/(32.D0*Pi**2)*MTR192(1,1))*MTR060(sld5))))/
     -    (-MZ2 + S)

#ifdef DEBUG
	DEB 'vert_sld5: Cloop(1) =', Cloop(1)
#endif

	end

