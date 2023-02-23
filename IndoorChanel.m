clear all;
O=[0,0,0];
NO=[10,10,3];
OX=[1,0,0];
OY=[0,1,0];
OZ=[0,0,1];
NXOYD=OZ;
NXOYU=-OZ;
NXOZR=OY;
NXOZL=-OY;
NYOZB=OX;
NYOZF=-OX;
M=[OX;OY;OZ];
NM=[NXOYD;NXOZR;NYOZB;NXOYU;NXOZL;NYOZF];
L=10;W=10;H=3;
T=[5,5,3];
nt=[0,0,-1];
R=[2.5,2.5,0];
nr=[0,0,1];
n=1;
N=500;%光子数
FOV=pi;%单位弧度
AR=1*10^(-4);%单位cm^2,注意单位！
G=1;
rho=0.8;
K=4;
nu=2*pi*rand(N,K);
mu=rand(N,K);
nu0=2*pi*rand(N,1);
mu0=rand(N,1);
P=zeros(N,K);%表示概率
Pik=zeros(3,N,K);%表示反射面上的光子的位置
Lik=zeros(N,K);%记录每个光子所走的总光程
e=zeros(N,K);
d=zeros(N,K);
c=3*10^8;%
num=900;
bi=10;
h=zeros(1,num+1);
% nik=zeros(N,K);

%***求Plos***%
Plos=((n+1)/(2*pi))*((R-T)*nt'/norm(R-T))*(AR*nr*(T-R)'/(norm(R-T))^3)*rect(1/FOV*acos((T-R)*nr'/norm(T-R)));

for i=1:N
    Lt=zeros(1,K);%记录一个光子每次反射所走的光程
    for k=1:K
        %***求P(i,k)***%
        lt=zeros(1,6);
        if k==1
            %第1面的光的方向矢量
            E0=[(1-mu0(i,1))^(0.5)*cos(nu0(i,1)),(1-mu0(i,1))^(0.5)*sin(nu0(i,1)),mu0(i,1)^(0.5)];
            e0=E0*M;
            for tt=1:6
                if(tt<=3)
                    if(O*NM(tt,:)'-T*NM(tt,:)'~=0&&NM(tt,:)*e0'~=0)
                        lt(tt)=(O*NM(tt,:)'-T*NM(tt,:)')/(NM(tt,:)*e0');
                    else
                        lt(tt)=9999;
                    end;
                end;
                if(tt>=4)
                    if(NO*NM(tt,:)'-T*NM(tt,:)'~=0&&NM(tt,:)*e0'~=0)
                        lt(tt)=(NO*NM(tt,:)'-T*NM(tt,:)')/(NM(tt,:)*e0');
                    else
                        lt(tt)=9999;
                    end;
                end;
            end;
            [tm,tn]=min(abs(lt));
            Pik(:,i,k)=T'+lt(tn)*e0';
            Lt(k)=abs(lt(tn));
            d=R-Pik(:,i,k)';
            nik=NM(tn,:);
            Lik(i,k)=sum(Lt(1,1:k))+norm(d);
        else
            %第k-1个反射面反射光的方向矢量
            Ek1=[(1-mu(i,k-1))^(0.5)*cos(nu(i,k-1)),(1-mu(i,k-1))^(0.5)*sin(nu(i,k-1)),mu(i,k-1)^(0.5)];
            ek1=Ek1*M;
            for tt=1:6
                if(tt<=3)
                    if(O*NM(tt,:)'-Pik(:,i,k-1)'*NM(tt,:)'~=0&&NM(tt,:)*ek1'~=0)
                        lt(tt)=(O*NM(tt,:)'-Pik(:,i,k-1)'*NM(tt,:)')/(NM(tt,:)*ek1');
                    else
                        lt(tt)=9999;
                    end;
                end;
                if(tt>=4)
                    if(NO*NM(tt,:)'-Pik(:,i,k-1)'*NM(tt,:)'~=0&&NM(tt,:)*ek1'~=0)
                        lt(tt)=(NO*NM(tt,:)'-Pik(:,i,k-1)'*NM(tt,:)')/(NM(tt,:)*ek1');
                    else
                        lt(tt)=9999;
                    end;
                end;
            end;
            [tm,tn]=min(abs(lt));
            Pik(:,i,k)=Pik(:,i,k-1)+lt(tn)*ek1';
            Lt(k)=abs(lt(tn));
            d=R-Pik(:,i,k)';
            nik=NM(tn,:);
            Lik(i,k)=sum(Lt(1,1:k))+norm(d);
        end;
        P(i,k)=((n+1)/(2*pi))*(d*nik'/norm(d))*AR*(-d*nr'/(norm(d))^3)*rect(1/FOV*acos(-d*nr'/(norm(d))));
        Gik=G/N*rho^k;
        for time=0:num
            timet=time/bi;
            h(time+1)=h(time+1)+Gik*P(i,k)*delta(timet,Lik(i,k)/c*10^(9));%换算成纳秒
        end;
    end;
end;

for time=0:num
    timet=time/bi;
    h(time+1)=h(time+1)+G*Plos*delta(timet,norm(R-T)/c*10^(9));
end;
timeLik=Lik/c*10^(9);
timet=0:1/bi:90;
plot(timet,h);
max=1.3*10^(-8);
axis([0 90 0 max]);
xlabel('t/ns');
ylabel('h(t)/W');
title('500光子数时多路径脉冲响应函数');
timemax=norm(R-T)/c*10^(9);
text(timemax,max*39/40,' \leftarrow LOS(6.1975e-7W)');