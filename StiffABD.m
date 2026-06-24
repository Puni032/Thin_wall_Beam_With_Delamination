function [A,B,D]=StiffABD(t,all_ang,C)
t= thickness(p)/tn(p);
all_ang = l(p,:);
N=length(all_ang);
h=zeros(N+1);
h(N+1)=N.*t/2;
theta=unique(all_ang);
nl=length(theta);
Qb=cell(nl,1);
lstf=zeros(3*N,3);
qq=zeros(3,3);
   for k = 1 : N
    
        %%Compliance Matrix elements
        S = [1/E1     -n21/E1    -n31/E1    0         0       0;
            -n21/E1   1/E2       -n23/E2    0         0       0;
            -n31/E1   -n23/E2     1/E3      0         0       0;
            0           0           0       1/G23     0       0;
            0           0           0        0       1/G13    0;
            0            0          0        0       0       1/G12];

        %%Stiffness Matrix elements
        C = inv(S);
        
%         [Stiffnew1]=Stiffnessnew(C,b(p),t,Stripcentroid(p,1),Stripcentroid(p,2),alpha(p));
%         Stiffnew=Stiffnew+Stiffnew1;
        %Reduced Stiffness
        Q = zeros(3,3);
        Q(1,1) = C(1,1) - (C(1,3))^2 / C(3,3);
        Q(2,2) = C(2,2) - (C(2,3))^2 / C(3,3);
        Q(1,2) = C(1,2) - C(1,3)*C(2,3) / C(3,3);
        Q(2,1) = Q(1,2);
        Q(3,3) = C(6,6);

%Transformatioin Matrix        

m = cosd (all_ang(k));
n = sind (all_ang(k));
       
        T = zeros(3,3);
        T(1,1)=m^2;
        T(1,2)=n^2;
        T(1,3)=2*m*n;
        T(2,1)=n^2;
        T(2,2)=m^2;
        T(2,3)=-2*m*n;
        T(3,1)=-m*n;
        T(3,2)=m*n;
        T(3,3)=m^2-n^2;
     
%Transformed Stiffness

q=zeros(3,3);
q([1:2],[1:2])=Q([1:2],[1:2]);
q(3,3)=2*Q(3,3);


q=T\q*T;


    for i = 1:3
        q(i,3)=q(i,3)/2;
    end
    qq(1,1,k)=q(1,1);
    qq(1,2,k)=q(1,2);
    qq(1,3,k)=q(1,3);
    qq(2,2,k)=q(2,2);
    qq(2,3,k)=q(2,3);
    qq(3,3,k)=q(3,3);
    qq(2,1,k)=qq(1,2,k);
    qq(3,1,k)=qq(1,3,k);
    qq(3,2,k)=qq(2,3,k);

%lstf([3*k-2:3*k],[1:3])=q([1:3],[1:3]);

h(k)=(k-N/2-1)*t;
  end
Qb{k,1}=q;

%lstf([3*k-2:3*k],[1:3])=q([1:3],[1:3]);



%To calculate ABD for laminate
% [A,B,D]=Stiffness_ABD(N,nl,l,theta,Qb,thickness(p))
A=zeros(3,3);
B=zeros(3,3);
D=zeros(3,3);
    
    for i=1:3
        for j=1:3
          for k = 1 : N
            A(i,j) = qq(i,j,k) * (t) + A(i,j);
            B(i,j) = 1/2*(qq(i,j,k) * (h(k+1)^2 - h(k)^2)) + B(i,j);
            D(i,j) = 1/3*(qq(i,j,k) * (h(k+1)^3 - h(k)^3)) + D(i,j);    
          end
       end
    end