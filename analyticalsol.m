function [Nsol]=analyticalsol(S)
S11=S(1,1);
S12=S(1,2);
S13=S(1,3);
S14=S(1,4);
S15=S(1,5);
S16=S(1,6);
S22=S(2,2);
S23=S(2,3);
S24=S(2,4);
S25=S(2,5);
S26=S(2,6);
S33=S(3,3);
S34=S(3,4);
S35=S(3,5);
S36=S(3,6);
S44=S(4,4);
S45=S(4,5);
S46=S(4,6);
S55=S(5,5);
S56=S(5,6);
S66=S(6,6);

l = 0.762;
Q3 = 0;
Q2 = 0;
M2 = 0;  M3 = 0;  Na = 0; 
M1 = 0.113;
x=0:0.0001:0.762;
kappa1=zeros(length(x),1);

for n=1:20
 for i=1:length(x)
F = [S11 -S14 -S13;-S14 S44 S34;-S13 S34 S33];
em = F\[S12 + S15*kappa1(i); -S24 - S45*kappa1(i);-S34 - S35*kappa1(i)];
fm = F\[S16;-S46;-S36];
gm = F\[0; Q3; -Q2];
hm = F\[Na; M2-l*Q3; M3 + l*Q2];
e = S22 + 2*S25*kappa1(i) +  S55*kappa1(i)^2 - em(1)*(S12 + S15*kappa1(i)) + em(2)*(S24 + S45*kappa1(i)) + em(3)*(S23 + S35*kappa1(i));
g = S66-fm(1)*S16 + fm(2)*S46 + fm(3)*S36; 
q = gm(1)*(S12 + S15*kappa1(i)) - gm(2)*(S24 + S45*kappa1(i)) - gm(3)*(S23 + S35*kappa1(i));
h = M1 + gm(1)*S16 - gm(2)*S46 - gm(3)*S36 - hm(1)*(S12 + S15*kappa1(i)) + hm(2)*(S24 + S45*kappa1(i)) + hm(3)*(S23 + S35*kappa1(i));
lambda = sqrt(e/g);
p1 = -q/(2*e);
p2 = h/e;
A2 = (p2*(exp(lambda*l) - 1)-2*p1*l)/(lambda*(exp(lambda*l)-exp(-lambda*l)));
A1 = A2 - p2/lambda;
p3 = -A1 - A2;
theta(i) = A1*exp(lambda*x(i)) + A2*exp(-lambda*x(i))+p1*x(i)^2+ p2*x(i)+p3;
kappa1(i)= A1*lambda*exp(lambda*x(i)) - A2*lambda*exp(-lambda*x(i))+p1*x(i)^2+p2*x(i)+p3;
 end
end
Nsol=plot(x/0.762,theta);
end