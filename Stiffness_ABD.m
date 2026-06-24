function [A,B,D]=Stiffness_ABD(N,nl,l,theta,Qb,tk)
Qbar=cell(N,1);
thick=tk/N;
for i=1:N
    for j=1:nl
        if l(i)==theta(j);
            Qbar{i,1}=Qb{j,1};
            break
        end
    end
end
z=zeros(N+1,1);
z(1)=(-tk/2);
for k=2:N+1
    z(k)=z(k-1)+thick;
end
A=zeros(3,3);
B=zeros(3,3);
D=zeros(3,3);
for i=1:N
     k=i+1;
        Aply=vpa(Qbar{i,1}*(z(k)-z(k-1)));
        Bply=vpa(Qbar{i,1}*((z(k))^2-(z(k-1))^2)*0.5);
        Dply=vpa((1/3)*Qbar{i,1}*((z(k))^3-(z(k-1))^3));
        A=A+Aply;
        B=B+Bply;
        D=D+Dply;
end
A=round((A*10^5))/10^5;
B=round(B*10^5)/10^5;
D=round(D*10^5)/10^5;
end
