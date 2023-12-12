clc;
clear;
close all;
%Bifurcation Diagrams... Beware this proyect is slow: takes time to record,
%and all that
dth=0;
tf=1800; dt=0.001; T=1.04;
t=[0:dt:tf];a0=0;ang0=pi/4;th=0; B=0.45; w=2*pi/T;w0=sqrt(w^2-B^2);w1=w0/1.6382; %A=0.009*(2*rand(1,length(t))-1)
N=length(t); D=86.1/0.01;
for j=1:D;
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
a0=a0+0.01;
G=th(350*T/dt:T/dt:N);Ga=dth(350/dt:T/dt:N);
plot(a0*ones(1,length(G)),G,'.b','MarkerSize',4)
axis([0 86 -5.1 1.6])
drawnow
hold on
 th=0;dth=0;
end
ylabel('Θ(rad)')
xlabel('a0(rad/s^2)')
saveas(gcf,'Bifur.png')
%%
a0=0;
for j=1:D;
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
a0=a0-0.05;
G=th(200/dt:T/dt:N);Ga=dth(200/dt:T/dt:N);
plot(abs(a0)*ones(1,length(G)),G,'.b','MarkerSize',4.2)
hold on
drawnow
end