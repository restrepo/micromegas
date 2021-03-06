#if 0
* vert_sqd5.F
* form factors for {F[6, {1}], F[6, {1}]} -> {V[1], V[1]}
* generated by FormCalc 8.0 on 6-Apr-2015 1:04
#endif

#include "vars.h"


	subroutine vert_sqd5
	implicit none

#include "vars.h"


#include "inline.h"

        Cloop(1) = Cloop(1) + 
     &    1/Pi**2*(MNE1*Pair3*AAABR(373)*
     &        (F1*(-(1/288.D0*pave13(bb0,sqd5)) - 
     &             1/144.D0*pave13(bb1,sqd5)) + 
     &          F2*(1/288.D0*pave13(bb0,sqd5) + 
     &             1/144.D0*pave13(bb1,sqd5)))*MTR041(sqd5)*
     &        (MTR139(1,1) + MTR140(1,1)))/(-MZ2 + S)

#ifdef DEBUG
	DEB("vert_sqd5: Cloop(1) =", Cloop(1))
#endif

	end
