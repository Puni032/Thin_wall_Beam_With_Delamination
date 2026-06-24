function [u1,t1,u2,u3]=Analytical1(Ksectionh,Ksection,Dcoord,Lcoord)

ldoc=length(Dcoord);
Fc=zeros(28*(ldoc+1));
uc=zeros(28*(ldoc+1));
Fq=zeros(28*(ldoc+1));
uq=zeros(28*(ldoc+1));

for i=1:ldoc+1
    
    [Fce,Fqe,uce,uqe]=Anasing(S,L);
    
    
end
    

   plot(xp/L,t1)
     
       
   