function [Ksection]=Sectional_stiffness(Properties,MId,th,ns,Xcg,alpha,conn,node,nd,thickness,k1,tn,Lconf,itnn,tlam,tdof,ln,lnm,dnum)

Mat1=zeros(tdof,6);
Mat2=zeros(tdof,4);
Mat3=zeros(tdof);

CowNm=zeros(4,tdof);
CoLNm=zeros(4,4);
CosNm=zeros(4,6);
A=cell(tlam,1);
B=cell(tlam,1);
D=cell(tlam,1);
cdof=4*ln;
lamnum=ns;

for p=1:ns
    n1=conn(p,1);
    n2=conn(p,2);
%     t= thickness(p)/tn(p);
Lns=Lconf{lnm(p,2)};
Mns=MId{(lnm(p,2))};
tnl=th{lnm(p,2)};
    if ismember(p,nd)
        
        for nsub=1:dnum(nd==p)
            if sum (bsxfun (@eq,conn(p,1), conn(:)))==1
                if nsub==1
                    edof=[4*n1-3 4*n1-2 4*n1-1 4*n1 4*n2-3 4*n2-2 4*n2-1 4*n2];
                    lnum=p;
                else
                    edof=[cdof+1 cdof+2 cdof+3 cdof+4 4*n2-3 4*n2-2 4*n2-1 4*n2];
                    cdof=cdof+4;
                    lnum=lamnum+1;
                    lamnum=lamnum+1;
                end
            end
            if sum (bsxfun (@eq,conn(p,2), conn(:)))==1
                if nsub==1
                    edof=[4*n1-3 4*n1-2 4*n1-1 4*n1 4*n2-3 4*n2-2 4*n2-1 4*n2];
                    lnum=p;
                else
                    edof=[4*n1-3 4*n1-2 4*n1-1 4*n1 cdof+1 cdof+2 cdof+3 cdof+4];
                    cdof=cdof+4;
                    lnum=lamnum+1;
                    lamnum=lamnum+1;
                end
            end
             if sum(bsxfun (@eq,conn(p,1), conn(:)))~=1 && sum(bsxfun (@eq,conn(p,2), conn(:)))~=1
                 if nsub==1
                    edof=[4*n1-3 4*n1-2 4*n1-1 4*n1 4*n2-3 4*n2-2 4*n2-1 4*n2];
                    lnum=p;
                else
                    edof=[4*n1-3 4*n1-2 4*n1-1 4*n1 4*n2-3 4*n2-2 4*n2-1 4*n2];
                    lnum=lamnum+1;
                    lamnum=lamnum+1;
                 end
             end
            itn=itnn{nd==p};
 
            tsub=tnl(itn(nsub)+1:itn(nsub+1));
            all_ang=Lns(1+itn(nsub):itn(nsub+1));
            Mnsub=Mns(1+itn(nsub):itn(nsub+1));
            Nsub=length(all_ang);
            N=tn(p,1);
            
            [As,Bs,Ds]=ABD_MATRIX(Properties,Mnsub,Nsub,all_ang,tsub);
            A{lnum,1}=As;
            B{lnum,1}=Bs;
            D{lnum,1}=Ds;
            hnd=-sign(itn(nsub)+itn(nsub+1)-N)*sum(tsub)/2;
            alp=alpha(p);
            Xn2=Xcg(p,1)-hnd*sin(alp);
            Xn3=Xcg(p,2)+hnd*cos(alp);

            a=(node(conn(p,1),1)-Xn2)*cos(alp)+(node(conn(p,1),2)-Xn3)*sin(alp);
            b=(node(conn(p,2),1)-Xn2)*cos(alp)+(node(conn(p,2),2)-Xn3)*sin(alp);
            
            hn=0;
            
           [wsvec,wcvec,wLvec,Consvec,Conscvec,ConsLvec,BCsvec,BCcvec,BCLvec]=WBCCEQ(As,Bs,Ds,hn,a,b,Xn2,Xn3,k1,alp);

           MATE1=BCcvec*(wcvec\wsvec)-BCsvec;
           MATE2=BCcvec*(wcvec\wLvec)-BCLvec;
           MATE3=BCcvec/(wcvec);
           
           CosNem=Consvec-Conscvec*(wcvec\wsvec);
           CowNem=Conscvec/wcvec;
           CoLNem=Conscvec*(wcvec\wLvec)-ConsLvec;
           Mat1(edof,:)=Mat1(edof,:)+MATE1;
           Mat2(edof,:)=Mat2(edof,:)+MATE2;
           Mat3(edof,edof)=Mat3(edof,edof)+MATE3;
           CosNm=CosNm+CosNem;
           CowNm(:,edof)=CowNm(:,edof)+CowNem;
           CoLNm=CoLNm+CoLNem;
        end
        
    else
        edof=[n1*4-3:n1*4,n2*4-3:n2*4];
        all_ang = Lns;
        N=length(all_ang);
       [As,Bs,Ds]=ABD_MATRIX(Properties,Mns,N,all_ang,tnl);
        A{p,1}=As;
        B{p,1}=Bs;
        D{p,1}=Ds;
        %[A,B,D]=StiffABD(t,all_ang,C);
        
        Xn2=Xcg(p,1);
        Xn3=Xcg(p,2);
        alp=alpha(p);
        a=(node(conn(p,1),1)-Xn2)*cos(alp)+(node(conn(p,1),2)-Xn3)*sin(alp);
        b=(node(conn(p,2),1)-Xn2)*cos(alp)+(node(conn(p,2),2)-Xn3)*sin(alp);
        hn=0;
        [wsvec,wcvec,wLvec,Consvec,Conscvec,ConsLvec,BCsvec,BCcvec,BCLvec]=WBCCEQ(As,Bs,Ds,hn,a,b,Xn2,Xn3,k1,alp);
        MATE1=BCcvec*((wcvec)\wsvec)-BCsvec;
        MATE2=BCcvec*((wcvec)\wLvec)-BCLvec;
        MATE3=BCcvec/(wcvec);
        CosNem=Consvec-Conscvec*((wcvec)\wsvec);
        CowNem=Conscvec/(wcvec);
        CoLNem=Conscvec*((wcvec)\wLvec)-ConsLvec;
        
        
        Mat1(edof,:)=Mat1(edof,:)+MATE1;
        Mat2(edof,:)=Mat2(edof,:)+MATE2;
        Mat3(edof,edof)=Mat3(edof,edof)+MATE3;
        CosNm=CosNm+CosNem;
        CowNm(:,edof)=CowNm(:,edof)+CowNem;
        CoLNm=CoLNm+CoLNem;
    end
end

Coglm=CoLNm\CosNm;
Conwlm=CoLNm\CowNm;

Kw=Mat3-Mat2*Conwlm;
Kr=Mat1+Mat2*Coglm;

we=(Kw)\Kr;
Lg=Coglm+Conwlm*we;
Ksection=zeros(6,6);
lamnum=ns;
cdof=4*ln;
for i=1:ns
    n1=conn(i,1);
    n2=conn(i,2);
    if ismember(i,nd)
        Ksd=zeros(6,6);
        for nsub=1:dnum(nd==i)
            if sum (bsxfun (@eq,conn(i,1), conn(:)))==1
                if nsub==1
                    edof=[4*n1-3 4*n1-2 4*n1-1 4*n1 4*n2-3 4*n2-2 4*n2-1 4*n2];
                    lnum=i;
                else
                    edof=[cdof+1 cdof+2 cdof+3 cdof+4 4*n2-3 4*n2-2 4*n2-1 4*n2];
                    cdof=cdof+4;
                    lnum=lamnum+1;
                    lamnum=lamnum+1;
                end
            end
            if sum (bsxfun (@eq,conn(i,2), conn(:)))==1
                if nsub==1
                    edof=[4*n1-3 4*n1-2 4*n1-1 4*n1 4*n2-3 4*n2-2 4*n2-1 4*n2];
                    lnum=i;
                
                else
                    edof=[4*n1-3 4*n1-2 4*n1-1 4*n1 cdof+1 cdof+2 cdof+3 cdof+4];
                    cdof=cdof+4;
                    lnum=lamnum+1;
                    lamnum=lamnum+1;
                end
            end
            if sum (bsxfun (@eq,conn(i,1), conn(:)))~=1 && sum(bsxfun (@eq,conn(i,2), conn(:)))~=1
                 if nsub==1
                    edof=[4*n1-3 4*n1-2 4*n1-1 4*n1 4*n2-3 4*n2-2 4*n2-1 4*n2];
                    lnum=i;
                else
                    edof=[4*n1-3 4*n1-2 4*n1-1 4*n1 4*n2-3 4*n2-2 4*n2-1 4*n2];
                    lnum=lamnum+1;
                    lamnum=lamnum+1;
                 end
            end
        hnd=-sign((itn(nsub)+itn(nsub+1)-N))*sum(tsub)/2;
        alp=alpha(i);
        Xn2=Xcg(i,1)-hnd*sin(alp); Xn3=Xcg(i,2)+hnd*cos(alp);
%         hn=-(itn(nsub)+itn(nsub+1)-N)/2*t*0;
        As=A{lnum,1};
        Bs=B{lnum,1};
        Ds=D{lnum,1};
        a=(node(conn(i,1),1)-Xn2)*cos(alp)+(node(conn(i,1),2)-Xn3)*sin(alp);
        b=(node(conn(i,2),1)-Xn2)*cos(alp)+(node(conn(i,2),2)-Xn3)*sin(alp);
        
        [wsvec,wcvec,wLvec,~,~,~,~,~,~]=WBCCEQ(As,Bs,Ds,hn,a,b,Xn2,Xn3,k1,alp);

        constant=(wcvec)\we(edof,:)-((wcvec)\wsvec)-((wcvec)\wLvec)*Lg;
        [Mgg,Mcg,Mlg,Mgc,Mcc,Mlc,Mgl,Mcl,Mll]=Energyterms(As,Bs,Ds,hn,a,b,Xn2,Xn3,k1,alp);
        
        Kss=Ksd+(Mgg+constant'*Mcg+Lg'*Mlg+Mgc*constant+constant'*Mcc*constant+Lg'*Mlc*constant+Mgl*Lg+constant'*Mcl*Lg+Lg'*Mll*Lg);
        Ksd=Ksd+(Mgg+constant'*Mcg+Lg'*Mlg+Mgc*constant+constant'*Mcc*constant+Lg'*Mlc*constant+Mgl*Lg+constant'*Mcl*Lg+Lg'*Mll*Lg);
        end
    else
        n1=conn(i,1);
        n2=conn(i,2);
        edof=[4*n1-3 4*n1-2 4*n1-1 4*n1 4*n2-3 4*n2-2 4*n2-1 4*n2];
        Xn2=Xcg(i,1); Xn3=Xcg(i,2);
        alp=alpha(i);
        hn=0;
        As=A{i,1};
        Bs=B{i,1};
        Ds=D{i,1};
        a=(node(conn(i,1),1)-Xn2)*cos(alp)+(node(conn(i,1),2)-Xn3)*sin(alp);
        b=(node(conn(i,2),1)-Xn2)*cos(alp)+(node(conn(i,2),2)-Xn3)*sin(alp);
        [wsvec,wcvec,wLvec,~,~,~,~,~,~]=WBCCEQ(As,Bs,Ds,hn,a,b,Xn2,Xn3,k1,alp);
        
        constant=wcvec\we(edof,:)-(wcvec\wsvec)-(wcvec\wLvec)*Lg;
        [Mgg,Mcg,Mlg,Mgc,Mcc,Mlc,Mgl,Mcl,Mll]=Energyterms(As,Bs,Ds,hn,a,b,Xn2,Xn3,k1,alp);
        
        Kss=(Mgg+constant'*Mcg+Lg'*Mlg+Mgc*constant+constant'*Mcc*constant+Lg'*Mlc*constant+Mgl*Lg+constant'*Mcl*Lg+Lg'*Mll*Lg);

                
    end
        Ksection=Ksection+Kss;
end