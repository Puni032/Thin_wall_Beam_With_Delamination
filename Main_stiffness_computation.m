%%Program for modal and transient analysis of pretwisted anisotropic strips based on spectral methods%%%
clc;
clear vars;
close all;
digits(32);
global Length breadth thickness rho
format short e;
disp('Welcome!')
%%Program for modal and transient analysis of pretwisted anisotropic strips based on spectral methods%%%
Ksections=zeros(6,6);
Inputfile=input('Data file name in txt format:','s');
fileID = fopen(Inputfile);
% read the whole file, interpret each line as a string
MyText = textscan(fileID, '%s%[^\n\r]', 'Delimiter', '', 'WhiteSpace', '',  'ReturnOnError', false);
Mytext=MyText{1};
Length=str2double(Mytext{3,1});

Nmater=str2num(Mytext{6,1});

% first the start locations
StartLineIdx = cellfun(@(x) contains(x, "Begin"),Mytext(:,1),'uni',0);
StartLineIdx = find([StartLineIdx{:}]==1) + 1;
% now find the stop locations of each data block
StopLineIdx = cellfun(@(x) contains(x, "End"),Mytext(:,1),'uni',0);
StopLineIdx = find([StopLineIdx{:}]) - 1;
% Allocate the cell array
InData = cell(numel(StartLineIdx),1);
% Get the data
for i = 1:numel(StartLineIdx)
    InData{i} = cellfun(@(x) textscan(x,'%f','Delimiter','\t','ReturnOnError',false),Mytext(StartLineIdx(i):StopLineIdx(i)));
    %     InData{i} = cellfun(@transpose,MyData{i},'UniformOutput',false);
    %     MyData{i} = cell2mat(MyData{i});
end

Properties=cell2mat(cellfun(@transpose,InData{1},'UniformOutput',false));

Matdis=InData{3};

th=InData{4};
tni=cellfun(@length,InData{4});
Lconf=InData{2};
lnm=cell2mat(cellfun(@transpose,InData{5},'UniformOutput',false));
thickn=cellfun(@sum,th);
thickness=thickn(lnm(:,2));
tn=tni(lnm(:,2));

node=cell2mat(cellfun(@transpose,InData{6},'UniformOutput',false));
% node(3,1)=(1-delp)*0.0127;
node=node(:,2:3);
xpoints=node(:,1);
ypoints=node(:,2);

conn=cell2mat(cellfun(@transpose,InData{7},'UniformOutput',false));
conn=conn(:,2:3);

nd=cell2mat(cellfun(@transpose,InData{8},'UniformOutput',false));
dnum=cell2mat(cellfun(@transpose,InData{9},'UniformOutput',false));

itnn=cellfun(@transpose,InData{10},'UniformOutput',false);

ns=length(conn(:,1));
alpha=zeros(ns,1);


ln=length(node(:,1));
% k1=3.94*0;
% tdof=4*ln;

Coeffmat=zeros(ns);
cns=1;
bcoeff=zeros(ns,1);

lnse=zeros(ns,1);
Xcg=zeros(ns,2);



for i=1:ns
    Xn2=(xpoints(conn(i,1))+xpoints(conn(i,2)))/2;
    Xn3=(ypoints(conn(i,1))+ypoints(conn(i,2)))/2;
    Xcg(i,1)=Xn2;
    Xcg(i,2)=Xn3;
    xc=(xpoints(conn(i,2))-xpoints(conn(i,1)));
    yc=(ypoints(conn(i,2))-ypoints(conn(i,1)));
    alphai=pi-pi/2*(1+sign(xc))*(1-sign(yc^2))-pi/4*(2+sign(xc))*sign(yc)-sign(xc*yc)*atan((abs(xc)-abs(yc))/(abs(xc)+abs(yc)));
    alpha(i)=alphai;
    lnse(i)=sqrt((ypoints(conn(i,2))-ypoints(conn(i,1)))^2+(xpoints(conn(i,2))-xpoints(conn(i,1)))^2);
    breadth(i)=lnse(i);
    Coeffmat(ns,i)=lnse(i);
end
% alpha=[0 pi/2 pi 3*pi/2];
alp=alpha;
for i=1:ln
    fln=find(conn==i);
    [r,~]=find(conn==i);
    if length(r)>1
        for j=1:length(r)-1
            e1=r(j);
            e2=r(j+1);
            Coeffmat(cns,e1)=1;
            Coeffmat(cns,e2)=-1;
            ye1=(node(i,1)-Xcg(e1,1))*cos(alp(e1))+(node(i,2)-Xcg(e1,2))*sin(alp(e1));
            bterm1=cos(alp(e1))*sin(alp(e1))*Xcg(e1,1)^2-cos(alp(e1))*sin(alp(e1))*Xcg(e1,2)^2-cos(alp(e1))^2*Xcg(e1,1)*Xcg(e1,2)+sin(alp(e1))^2*Xcg(e1,1)*Xcg(e1,2)+sin(alp(e1))*Xcg(e1,1)*ye1-cos(alp(e1))*Xcg(e1,2)*ye1;
            ye2=(node(i,1)-Xcg(e2,1))*cos(alp(e2))+(node(i,2)-Xcg(e2,2))*sin(alp(e2));
            bterm2=cos(alp(e2))*sin(alp(e2))*Xcg(e2,1)^2-cos(alp(e2))*sin(alp(e2))*Xcg(e2,2)^2-cos(alp(e2))^2*Xcg(e2,1)*Xcg(e2,2)+sin(alp(e2))^2*Xcg(e2,1)*Xcg(e2,2)+sin(alp(e2))*Xcg(e2,1)*ye2-cos(alp(e2))*Xcg(e2,2)*ye2;
            bcoeff(cns,1)=bterm2-bterm1;
            cns=cns+1;
        end
    end
end

k1=0;
tdof=4*ln;
tlam=ns;
tns=ns;

if nd(1)~=0
    for bn=1:length(nd)
        if sum (bsxfun (@eq,conn(nd(bn),1), conn(:)))==1
            tdof=tdof+4*(dnum(bn)-1);
        end
        if sum (bsxfun (@eq,conn(nd(bn),2), conn(:)))==1
            tdof=tdof+4*(dnum(bn)-1);
        end
        tlam=tlam+(dnum(bn)-1);
    end
end

ndh=0;
[Ksectionh]=Sectional_stiffness(Properties,Matdis,th,ns,Xcg,alpha,conn,node,ndh,thickness,k1,tn,Lconf,itnn,ns,ln*4,ln,lnm,dnum);

[Ksection]=Sectional_stiffness(Properties,Matdis,th,ns,Xcg,alpha,conn,node,nd,thickness,k1,tn,Lconf,itnn,tlam,tdof,ln,lnm,dnum);

%Delamination coordinates along the length
Dcoord=[Length Length/2];
Lcoord=unique([0 Dcoord Length]);

[Nsol]=Anasing(Ksectionh,Length);
% [Nsold]=Anasing(Ksection,Length);


disp('----------------------------------------------------------')
disp(' What is the pretwist value?');
twist=0; %input('Enter pretwist (rad/meter)>> ');
% STIFFFIN=TotalStiffnessnew(StiffIb,breadth,alpha,Xcg,ns,twist);
% STIFFIN=[5.133*10^4 0 0 15.072 190.92 0; 0 0.178 -1.278 0 0.0165/2 1.145*10^(-5); 0 -1.278 154.54 0 0 -5.144*10^(-4); 15.072 0 0 227.146 8.65*10^(-4) 0; 190.92 0.0165/2 0 8.65*10^(-4) 0.0107 0; 0 1.145*10^(-5) -5.144*10^(-4) 0 0 7.752*10^(-3)];
% [u]=static(Length,STIFFFIN);
% [Nsol]=analyticalsol(STIFFFIN);
% disp('----------------------------------------------------------')
% disp('What analysis do you want to perform?');
% disp('1 MODAL');
% disp('2 TRANSIENT');
% result=input('Enter the number correspond to your choice>> ');

% STIFFLIN=STIFFFIN(1:4,1:4);
% STIFFd=STIFFLIN;