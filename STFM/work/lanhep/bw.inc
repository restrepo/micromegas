
% Dark Sector   1910.01549


parameter MG1=100,MG2=102, LAM=10000, kappa=0.1, kappa2=0.2,lmbd=0.3.

parameter MM11 = (+1*MG1+2/EE**2/LAM*MW**2*SW**2*kappa2).
parameter MM12 = (-1/EE**2/LAM*MW**2*SW**2*lmbd)      .
parameter MM22 = (+1*MG2+2/EE**2/LAM*MW**2*SW**2*kappa) .

parameter hDiag= rDiagonal(2,MM11,MM12,MM22).
parameter MNE1= MassArray(hDiag,1), MNE2 = MassArray(hDiag,2). 
parameter MC1=MG2+2/EE**2/LAM*MW**2*SW**2*kappa.

_j=1-2 in parameter Z1_j=MixMatrix(hDiag,_j,1).
_j=1-2 in parameter Z2_j=MixMatrix(hDiag,_j,2).


OrthMatrix({  {Z11, Z12}, {Z21, Z22}  }).

%let Z11=1,Z12=0,Z21=0,Z22=0.


spinor    ~o1 / ~o1:('neutralino 1', pdg 1000022,mass MNE1),
          ~o2 / ~o2:('neutralino 2', pdg 1000023,mass MNE2, width wNE2=auto),
         '~1+'/'~1-':('chargino 1', pdg 1000024,mass MC1,  width wC1=auto).


let f_W   =  { ('~1+'+cc('~1-'))/Sqrt2, -('~1+'-cc('~1-'))/Sqrt2/i,   Z21*   '~o1'  + Z22*'~o2'     } .
let f_W_T =  { ('~1-'+cc('~1+'))/Sqrt2,  ('~1-'-cc('~1+'))/Sqrt2/i,   Z21*cc('~o1') + Z22*cc('~o2') } .

lterm   i/2*f_W_T*gamma*(deriv +g2*eps*W)*f_W.
lterm -MG2/2*f_W_T *f_W.

let Bino = Z11*'~o1'+Z12*'~o2', BinoT = Z11*cc('~o1') + Z12*cc('~o2').
lterm   i/2* BinoT*gamma*deriv*Bino.
lterm   -MG1/2*BinoT*Bino.


lterm  - kappa/2/LAM*f_W_T*f_W*shD*shd .
lterm  - kappa2/2/LAM*BinoT*Bino*shD*shd  .
lterm  - lmbd/LAM*f_W_T*Bino * shD*tau*shd/2.

