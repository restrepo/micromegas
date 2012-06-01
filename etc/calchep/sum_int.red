
procedure initSum();
begin
 
   sum:=0;
let p1.p2= Eg*Mh3;
let p1.p3=Mh3*(1/2*(Mh3-Eg) + Eg*kappa*cosfi/2);

%let kappa^2= 1- 4*mf^2/Mh/(Mh-2*Eg);


   for all p, m,w let propDen(p,m,w)= -1/(p.p-m^2);

end;

procedure s_subst(s,x) ;
begin
   return (s where cosfi => x);
end; 
              

procedure addToSum();
begin

  int_fun:=int(numerator*denominator ,cosfi);

  int_fun:=int_fun;
  sum:=sum+  totFactor*(s_subst(int_fun,1)- s_subst(int_fun,-1));
end$

procedure finishSum();
begin

  on combinelogs;
  sum:=sum;
  off combinelogs;
   let kappa^2= 1- 4*MW^2/Mh3/(Mh3-2*Eg);
let Eg=Mh3/2*x;
let MW=eps*Mh3/2;
  sum:=sum*kappa;
off nat;
  on factor;

%  if not (sub_S=0) then  sum:=(sum where log(~x)=>log(sub({s=sub_S},x)));
%  sum:= (sum where log(~x)=>-log(1/x) when  ordp(x,1/x));
%  sum:=(sum/((s- p1.p1-p2.p2)^2 -4*p1.p1*p2.p2)/(16*pi) where substitutions);
end;

end;

