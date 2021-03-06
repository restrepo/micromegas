#if 0
* gzvert_sqd5.F
* form factors for {F[6, {1}], F[6, {1}]} -> {V[1], V[2]}
* generated by FormCalc 8.0 on 6-Apr-2015 0:18
#endif

#include "gzvars.h"


	subroutine gzvert_sqd5
	implicit none

#include "gzvars.h"


#include "inline.h"

        Cloop(1) = Cloop(1) + 
     &    1/Pi**2*(MNE1*Pair3*
     &        (F1*(-(1/576.D0*pave19(bb0,sqd5)) - 
     &             1/288.D0*pave19(bb1,sqd5)) + 
     &          F2*(1/576.D0*pave19(bb0,sqd5) + 
     &             1/288.D0*pave19(bb1,sqd5)))*MTR041(sqd5)*
     &        (MTR139(1,1) + MTR140(1,1))*MTR152(sqd5))/(-MZ2 + S)

#ifdef DEBUG
	DEB("vert_sqd5: Cloop(1) =", Cloop(1))
#endif

	end
