function acel=proy(w,w0,w1,a0,ang0,ang,t,dang,B);
acel=-w0^2*sin(ang)+a0*cos(w*t)-w1^2*(ang-ang0)-dang*2*B;
end