%    ==============================
%    *  CalcHEP version 2.5b     *
%    ==============================
inParticles:={"q1","~f"}$
outParticles:={"q1","~f"}$
%
parameters:={
 Mq=>1.000000E-02
,MX=>1.000000E+02
,Seven=>1.000000E-02
,MSQ=>1.000000E+00
,A=>1.000000E+00
,B=>1.000000E+00
}$
%
substitutions:={
}$


vector p1,p2,p3,p4,p5,p6$

let p4 = +p1+p2-p3;
mass p1  = Mq; Mshell p1;
mass p2  = MX; Mshell p2;
mass p3  = Mq; Mshell p3;
let p2.p3 = -1*(MX^2-Mq^2-MX^2-Mq^2-2*p1.p2+2*p1.p3)/2;
let p1.p3=Mq^2;

vector !=p_,!=q_$
operator propDen$
for all p_,q_,m,w let propDen(0*p_+q_,m,w)=propDen(q_,m,w)$
for all p_,m,w such that ordp(p_,-p_) let propDen(p_,m,w)=propDen(-p_,m,w);$

initSum();

DiagrNumber:="1_2"$

%                  q1    q1   !  q1          q1                      
%                ==>==@==>====!==>==\     /==>==                     
%                  P1 |  P3   !  P3 |     |  P1                      
%                   s0|       !     |     |                          
%                  ~f |  ~f   !  ~f | ~Sq |  ~f                      
%                =====@=======!=====@-->--@=====                     
%                  P2    P4   !  P4          P2                      
totS:=(8*Seven)/(1)$
numS:=p1.p3*p1.p2*B^2-p1.p3*p1.p2*A^2+p1.p3*B^2*MX^2-p1.p3*B^2*MX*Mq+
 p1.p3*B^2*Mq^2-p1.p3*A^2*MX^2-p1.p3*A^2*MX*Mq-p1.p3*A^2*Mq^2+2*p1.p2*B^2*
 MX*Mq-p1.p2*B^2*Mq^2+2*p1.p2*A^2*MX*Mq+p1.p2*A^2*Mq^2+B^2*MX^2*Mq^2+B^2*MX*
 Mq^3-B^2*Mq^4-A^2*MX^2*Mq^2+A^2*MX*Mq^3+A^2*Mq^4$
denS:=MSQ^2-(p1+p2).(p1+p2);


addToSum()$

DiagrNumber:="1_3"$

%                                !  ~f    q1                         
%                               /!=====@==>==                        
%                               |!  P4 |  P1                         
%                               |!  ~Sq|                             
%                     q1    q1  |!  q1 |  ~f                         
%                   ==>==@==>===+!==>==@=====                        
%                     P1 |  P3  |!  P3    P2                         
%                      s0|      |!                                   
%                     ~f |  ~f  |!                                   
%                   =====@======/!                                   
%                     P2    P4   !                                   
totU:=(8*Seven)/(1)$
numU:=p1.p3^2*B^2-p1.p3^2*A^2-p1.p3*p1.p2*B^2+p1.p3*p1.p2*A^2+p1.p3*B^
 2*MX^2+p1.p3*B^2*MX*Mq-p1.p3*B^2*Mq^2-p1.p3*A^2*MX^2+p1.p3*A^2*MX*Mq+p1.p3*
 A^2*Mq^2-2*p1.p2*B^2*MX*Mq+p1.p2*B^2*Mq^2-2*p1.p2*A^2*MX*Mq-p1.p2*A^2*Mq^2+
 B^2*MX^2*Mq^2-B^2*MX*Mq^3-A^2*MX^2*Mq^2-A^2*MX*Mq^3$
denU:=MSQ^2-(p1-p4).(p1-p4);


addToSum()$
finishSum();
End$
