
gam   =      1; 
beta2=0.01;
delta2=1e-4;
N     =      1e5;
Avg   =     10;  % Number of Realizations
 %x1    = stblrnd(1.85,beta,gam,delta,N,Avg); 
 load("input_185.mat")
 Ref_Signal=x1';
mew=[1e-1 1e-2 1e-3 1e-4 1e-5 1e-6 1e-7];
lambda=0.9999;
delta=0.3;
filename='P1Z.BIN'; fid=fopen(filename,'r'); pz1_num=fread(fid,'float'); fclose(fid);
% denominator of P1(z)
filename='P1P.BIN'; fid=fopen(filename,'r'); pz1_den=fread(fid,'float'); fclose(fid);
filename='P2Z.BIN'; fid=fopen(filename,'r'); pz2_num=fread(fid,'float'); fclose(fid);
% denominator of P1(z)
filename='P2P.BIN'; fid=fopen(filename,'r'); pz2_den=fread(fid,'float'); fclose(fid);
% numerator of S11(z)
filename='S11Z.BIN'; fid=fopen(filename,'r'); sz11_num=fread(fid,'float'); fclose(fid);
% denominator of S11(z)
filename='S11P.BIN'; fid=fopen(filename,'r'); sz11_den=fread(fid,'float'); fclose(fid);
filename='S12Z.BIN'; fid=fopen(filename,'r'); sz12_num=fread(fid,'float'); fclose(fid);
% denominator of S11(z)
filename='S12P.BIN'; fid=fopen(filename,'r'); sz12_den=fread(fid,'float'); fclose(fid);
filename='S21Z.BIN'; fid=fopen(filename,'r'); sz21_num=fread(fid,'float'); fclose(fid);
% denominator of S11(z)
filename='S21P.BIN'; fid=fopen(filename,'r'); sz21_den=fread(fid,'float'); fclose(fid);
filename='S22Z.BIN'; fid=fopen(filename,'r'); sz22_num=fread(fid,'float'); fclose(fid);
% denominator of S11(z)
filename='S22P.BIN'; fid=fopen(filename,'r'); sz22_den=fread(fid,'float'); fclose(fid);
IMP_PZ1  =   impz(pz1_num,pz1_den);
IMP_PZ2  =   impz(pz2_num,pz2_den);
IMP_SZ11 =   impz(sz11_num,sz11_den);
IMP_SZ12 =   impz(sz12_num,sz12_den);
IMP_SZ21 =   impz(sz21_num,sz21_den);
IMP_SZ22 =   impz(sz22_num,sz22_den);
L     =   256;%256;    % The primapry path
M     =   128;% 32;    % The secondary Path
Lw    =   192; %64;    % The ANC filter W(z)
pz1 =   IMP_PZ1(1:L);
pz2  =  IMP_PZ2(1:L);
sz11  = IMP_SZ11(1:M);
sz12 =  IMP_SZ12(1:M);
sz21  = IMP_SZ21(1:M);
sz22=   IMP_SZ22(1:M);
sh11=sz11;
sh12=sz12;
sh21=sz21;
sh22=sz22;
for b=1:length(mew)
s = num2str(mew);
EE1h    =       double(zeros(Avg,N));    % for all realizations
%EE1   =       double(zeros(N,DR)); 
PEE1     =       double(zeros(Avg,N));   % power of e(n) for all realizations
D1       =       double(zeros(1,N));          % averaged
DD1      =       double(zeros(Avg,N));    % for all realizations
PDD1     =       double(zeros(Avg,N));   % power of e(n) for all realizations
NRR1     =       double(zeros(Avg,N));  % Noise reduction for all realizations
EE2h    =       double(zeros(Avg,N));    % for all realizations
%EE2   =       double(zeros(N,DR));    % for all realizations
PEE1     =       double(zeros(Avg,N));   % power of e(n) for all realizations
D2       =       double(zeros(1,N));          % averaged
DD2      =       double(zeros(Avg,N));    % for all realizations
PDD2     =       double(zeros(Avg,N));   % power of e(n) for all realizations
NRR2     =       double(zeros(Avg,N));  % Noise reduction for all realizations
E1      =       double(zeros(1,N)); 
E2=     double(zeros(1,N)); 
E1h      =       double(zeros(1,N)); 
E2h=     double(zeros(1,N)); 
for avg=1:Avg  
x= Ref_Signal(avg,:);
xh11=zeros(1,L);
xh12=zeros(1,L);
xh21=zeros(1,L);
xh22=zeros(1,L);
xx=zeros(1,L);
y11=zeros(1,M);
y12=zeros(1,M);
y21=zeros(1,M);
y22=zeros(1,M);
y1=zeros(1,M);
y2=zeros(1,M);
y1a=zeros(1,M);
y2a=zeros(1,M);
w1=zeros(Lw,1);
d1=zeros(1,N);
d2=zeros(1,N);
w2=zeros(Lw,1);
e1a=zeros(1,N);
e2a=zeros(1,N);

for n=1:N
  
    xx=[x(n) xx(1:L-1)]; 
    d1(n)=pz1'*xx';
    d2(n)=pz2'*xx';
    y1 = [w1'*xx(1:Lw)'    y1(1:M-1)];
    y2 = [w2'*xx(1:Lw)'    y2(1:M-1)];
        y11(n)=y1*sz11;
       y12(n)=y2*sz12;
       y21(n)=y1*sz21;
       y22(n)=y2*sz22;
        e1a(n) = d1(n)-y11(n)-y12(n);
       e2a(n)= d2(n)-y22(n)-y21(n);
       xh11=[sh11'*xx(1:M)' xh11(1:L-1)];
       xh12=[sh21'*xx(1:M)' xh12(1:L-1)];
       xh21=[sh12'*xx(1:M)' xh21(1:L-1)];
       xh22=[sh22'*xx(1:M)' xh22(1:L-1)];
        E1(n+1)=(lambda*E1(n))+(1-lambda)*(abs(e1a(n))^2);
         E2(n+1)=(lambda*E2(n))+(1-lambda)*(abs(e2a(n))^2);
         mew1 = mew(b)  / ((norm(xh11+xh12)^2)+delta2+E1(n+1));
          mew2 = mew(b)  / ((norm(xh21+xh22)^2)+delta2+E2(n+1));
        mu1= mew1/(1 + beta2* (abs(e1a(n))^3));
        mu2= mew2/(1 + beta2* (abs(e2a(n))^3));
        w1  = w1 + ( mu1 * (e1a(n)^2)*sign(e1a(n)) ) * xh11(1:Lw)'+ ( mu2 * (e2a(n)^2)*sign(e2a(n)) ) * xh12(1:Lw)';
         w2  = w2 + ( mu1 * (e1a(n)^2)*sign(e1a(n)) ) * xh21(1:Lw)'+ ( mu2 * (e2a(n)^2)*sign(e2a(n)) ) * xh22(1:Lw)';
end
   E1h = E1h + (abs(e1a));
    EE1h(avg,:) = e1a;
    D1 = D1 + (abs(d1));
    DD1(avg,:) = d1;
     E2h = E2h + (abs(e2a));
    EE2h(avg,:) = e2a;
    D2 = D2 + (abs(d2));
    DD2(avg,:) = d2;
end
%EW1(b,:) = EE1;
 %EW2(b,:) = EE2;
%     D1 = D1 + (abs(d1));
    %DW1(b,:) = DD1;
%     D2 = D2 + (abs(d2));%DW2(b,:) = DD2;
% Ew1 = (1/Avg) * EE1;
% D1 = (1/Avg) * DD1;
% Ew2 = (1/Avg) * EE2;f

% +-------------------------------------------------+
for avg = 1 : Avg
    Temp1 = double(zeros(1,N));
    P=1;
    for n = 1:N
        P=0.99*P+(1-0.99)*abs(DD2(avg,n));
        Temp1(n)=P; %Power estimate of averaged absolute error
    end
    PDD(avg,:) = Temp1;
    
    Temp2 = double(zeros(1,N));
    P=1;
    for n = 1:N
        P=0.99*P+(1-0.99)*abs(EE2h(avg,n));
        Temp2(n)=P; %Power estimate of averaged absolute error
    end
    PEE(avg,:) = Temp2;
    
    NRR1(avg,:) = Temp2./Temp1;
    
end

%plot(e1a, 'linewidth', 2)
for avg = 1 : Avg
    Temp1 = double(zeros(1,N));
    P=1;
    for n = 1:N
        P=0.99*P+(1-0.99)*abs(DD1(avg,n));
        Temp1(n)=P; %Power estimate of averaged absolute error
    end
    PDD(avg,:) = Temp1;
    
    Temp2 = double(zeros(1,N));
    P=1;
    for n = 1:N
        P=0.99*P+(1-0.99)*abs(EE1h(avg,n));
        Temp2(n)=P; %Power estimate of averaged absolute error
    end
    PEE(avg,:) = Temp2;
    
    NRR2(avg,:) = Temp2./Temp1;
    
end

MNR2 = 1/Avg * sum(NRR1,1);
MNR1 = 1/Avg * sum(NRR2,1);
subplot(2,1,1)
plot(MNR1, 'linewidth', 2);
title(' E1')
hold on
ylim([0 1.1])
legend('µ=1e-1','µ=1e-2','µ=1e-3','µ=1e-4','µ=1e-5','µ=1e-6','µ=1e-7')
subplot(2,1,2)
plot(MNR2, 'linewidth', 2)
title(' E2')
hold on
ylim([0 1.1])
legend('µ=1e-1','µ=1e-2','µ=1e-3','µ=1e-4','µ=1e-5','µ=1e-6','µ=1e-7')
%title('error')
end