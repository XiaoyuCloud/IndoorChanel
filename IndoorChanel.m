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
N=500;%������
FOV=pi;%��λ����
AR=1*10^(-4);%��λcm^2,ע�ⵥλ��
G=1;
rho=0.8;
K=4;
nu=2*pi*rand(N,K);
mu=rand(N,K);
nu0=2*pi*rand(N,1);
mu0=rand(N,1);
P=zeros(N,K);%��ʾ����
Pik=zeros(3,N,K);%��ʾ�������ϵĹ��ӵ�λ��
Lik=zeros(N,K);%��¼ÿ���������ߵ��ܹ��
e=zeros(N,K);
d=zeros(N,K);
c=3*10^8;%
num=900;
bi=10;
h=zeros(1,num+1);
% nik=zeros(N,K);

%***��Plos***%
Plos=((n+1)/(2*pi))*((R-T)*nt'/norm(R-T))*(AR*nr*(T-R)'/(norm(R-T))^3)*rect(1/FOV*acos((T-R)*nr'/norm(T-R)));

for i=1:N
    Lt=zeros(1,K);%��¼һ������ÿ�η������ߵĹ��
    for k=1:K
        %***��P(i,k)***%
        lt=zeros(1,6);
        if k==1
            %��1��Ĺ�ķ���ʸ��
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
            %��k-1�������淴���ķ���ʸ��
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
            h(time+1)=h(time+1)+Gik*P(i,k)*delta(timet,Lik(i,k)/c*10^(9));%���������
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
title('500������ʱ��·��������Ӧ����');
timemax=norm(R-T)/c*10^(9);
text(timemax,max*39/40,' \leftarrow LOS(6.1975e-7W)');