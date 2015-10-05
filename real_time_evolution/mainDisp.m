clear all,
clc

load SystemParameters;
load Measurment;

%======================================================================
%======================================================================
%======================================================================
figure;plot(space,randV,'-or')

figure;plot(time,HalfNum,'-');

figure;semilogx(time,Mdepart,'-r');
figure;semilogx(time,Mentropy,'-b');

% figure;plot(space,Mdensity1f(:,Nt),'-sr');hold on;plot(space,Mdensity2f(:,Nt),'-ob');
% xlim([1,length(space)]);ylim([0,1.1]);



% figure;surf(Mdensity2f');shading interp
% view([0,90]);
% xlim([1,L]);
% ylim([0,50]);
% yTicks=[0,10,20,30,40,50];
% yTicks2=5+yTicks*Nt;
% % yTicks2=yTicks;
% set(gca,'ytick',yTicks,'yticklabel',yTicks2,'fontsize',14);
% title('W=10,U=5 (\chi=30)')
% xlabel('site','fontsize',14);
% ylabel('time (units of J^{-1})','fontsize',14);

% a=Mdensity2f;
% M=moviein(Nt);
% figure;
% for j=1:Nt,
% plot(space,a(:,j),'-o');
% M(:,j)=getframe;
% end
% figure;movie(M,1);
% % movie2avi(M,'showL6beta');

a0=Mdensity2f0;b0=Mdensity1f0;
a1=Mdensity2f;b1=Mdensity1f;
M=moviein(Nt+Nt0);
figure;
for j=1:Nt0,
    plot(space,a0(:,j),'-o',space,b0(:,j),'-vr');ylim([0,1]);
    M(:,j)=getframe;
end
for j=1:Nt,
plot(space,a1(:,j),'-o',space,b1(:,j),'-vr');ylim([0,1]);
M(:,j+Nt0)=getframe;
end
% figure;movie(M,1);
movie2avi(M,['showU',num2str(U),'W',num2str(W),'.avi']);

densteady=zeros(L,1);
for j=0:9,
    densteady=densteady+Mdensity2f(:,Nt-j);
end
densteady=densteady/10;
% figure;
% plot(space,densteady,'-ob');ylim([0,1]);xlim([1,L]);
hold on;
plot(space,densteady,'-');ylim([0,1]);xlim([1,L]);
xlabel('site');ylabel('\rho_{steady}');title(['W=',num2str(W)]);
set(gca,'fontsize',12);

figure;
subplot(3,1,1);
plot(time,real(EimpKin),'-r');hold on;plot(time,real(EimpKin(1))+time*0,'-.r');
hold on;plot(time,real(EimpPot),'-b');hold on;plot(time,real(EimpPot(1))+time*0,'-.b');
subplot(3,1,2);
plot(time,real(EmajKin),'-r');hold on;plot(time,real(EmajKin(1))+time*0,'-.r');
hold on;plot(time,real(EmajPot),'-b');hold on;plot(time,real(EmajPot(1))+time*0,'-.b');
subplot(3,1,3);
plot(time,real(Eint),'-k');hold on;plot(time,real(Eint(1))+time*0,'-.k');
hold on;plot(time,real(EimpKin+EimpPot+EmajKin+EmajPot+Eint),'-g');

