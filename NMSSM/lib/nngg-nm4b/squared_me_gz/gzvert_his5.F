#if 0
* gzvert_his5.F
* form factors for {F[6, {1}], F[6, {1}]} -> {V[1], V[2]}
* generated by FormCalc 8.0 on 6-Apr-2015 0:18
#endif

#include "gzvars.h"


	subroutine gzvert_his5
	implicit none

#include "gzvars.h"


#include "inline.h"

        Cloop(1) = Cloop(1) + 
     &    (1/Pi**2*(1/32.D0*(Pair4*Pair5*pave34(cc12)*Sub417(his5) - 
     &             pave34(cc1)*Sub414(his5)*Sub594(his5) - 
     &             pave34(cc2)*Sub414(his5)*Sub597(his5) + 
     &             pave34(cc0)*Sub414(his5)*Sub604(his5)) + 
     &          EE*AAABR(175)*
     &           (-(1/8.D0*(Pair3*pave31(cc00))) + 
     &             1/8.D0*(Pair4*Pair5*pave31(cc12)))*MTR003(his5)*
     &           Sub414(his5)) + 
     &       Pair3*(1/Pi**2*
     &           (-(1/32.D0*(pave34(cc00)*Sub417(his5))) + 
     &             1/32.D0*((AAABR(372)*MTR003(his5)*pave30(bb0) + 
     &                   AAABR(232)*MTR188(his5)*pave32(bb0) + 
     &                   AAABR(5)*MTR189(his5)*pave33(bb0))*
     &                 Sub414(his5) + pave35(bb0)*Sub613(his5))) - 
     &          1/8.D0*(Sub414(his5)*Sub614(his5))))/
     &     (S - hisMass(his5)**2)

#ifdef DEBUG
	DEB("vert_his5: Cloop(1) =", Cloop(1))
#endif

	end
