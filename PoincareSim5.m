clc;
clear;
close all;
%Poincaré diagrams!
dth=0;
tf=3800; dt=0.01; T=0.84;
t=[0:dt:tf];a0=30.20;ang0=pi/4;th=0; B=0.05; w=2*pi/T;w0=sqrt(w^2-B^2);w1=w0/1.73; A=0.009*(2*rand(1,length(t))-1)
N=length(t);
for i=1:N-1;
k1t=dth(i);
    k1dt=proy(w,w0,w1,a0,ang0,th(i),t(i),dth(i),B);
    %
    k2t=dth(i)+k1dt*dt/2;
      k2dt=proy(w,w0,w1,a0,ang0,th(i)+k1t*dt/2,t(i)+dt/2,dth(i)+dt/2*k1dt,B);
    %
    k3t=dth(i)+k2dt*dt/2;
      k3dt=proy(w,w0,w1,a0,ang0,th(i)+k2t*dt/2,t(i)+dt/2,dth(i)+dt/2*k2dt,B);
    %
    k4t=dth(i)+k3dt*dt;
   k4dt=proy(w,w0,w1,a0,ang0,th(i)+k3t*dt,t(i)+dt,dth(i)+dt*k3dt,B);
    %
    th(i+1)=th(i)+(dt/6)*(k1t+2*k2t+2*k3t+k4t);
    dth(i+1)=dth(i)+(dt/6)*(k1dt+2*k2dt+2*k3dt+k4dt);
end
%%
th=round(th*180/pi*100)/100*pi/180;dth=round(dth*180/pi*100)/100*pi/180;
th=th+0.0035*pi*(2*rand(1,length(th))-1);dth=dth+0.0035*pi*(2*rand(1,length(th))-1);th(1)=0;
%%
Y=abs(fft(th))/N; f=100*(0:length(Y)-1)/length(Y);
figure(333)
plot(f,Y)
%%
th=th-pi/3.6;
th=th/1.2;
dth=dth/2;
figure(1)
plot(th,dth,'-b','LineWidth',1.1)
axis([min(th)-0.12 max(th)+0.08 -max(dth)-0.4 max(dth)+0.4])
xlabel('Θ (rad)')
ylabel('ω (rad/s)')
%%
inc=84; delta=1;
figure(2)
pause(1)
for i=1:150;
plot(th(delta:inc:end),dth(delta:inc:end),'.b','MarkerSize',10)
axis([min(th)-0.12 max(th)+0.08 -max(dth)-1.39 max(dth)+0.4])
xlabel('Θ (rad)')
ylabel('ω (rad/s)')
delta=delta+1;
title(['δ= ',num2str(2*pi*delta/inc),' rad'])
Fa(i)=getframe(gcf);
hold off
end
video=VideoWriter('test1');
video.Quality=100;
video.FrameRate=20;
open(video);
writeVideo(video,Fa);
close(video);
%%
figure(3);
plot(t,th,'-b')
hold on
plot(t(1:inc:end),th(1:inc:end),'.r','MarkerSize',14)
ylabel('Θ (rad)')
xlabel('t (s)')
hold off
figure(4)
D=abs(fft(th)/N);
plot(1./t,D)
figure(5);
plot(t,dth,'-b')
hold on
plot(t(1:inc:end),dth(1:inc:end),'.r','MarkerSize',14)
ylabel('ω (rad)')
xlabel('t (s)')
hold off