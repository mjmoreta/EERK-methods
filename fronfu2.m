function valz=fronfu2(tt,kk,xx,yy)

% Function that gives the boundary value 
% f(tn+k2,u(tn)+k2 ut(tn)+k2^2/4 Aut(tn))

aux1=(cos(tt+xx+yy)-kk*sin(tt+xx+yy)+(kk^2)*sin(tt+xx+yy))^2;
aux2=2*cos(tt+kk+xx+yy)-sin(tt+kk+xx+yy)-cos(tt+kk+xx+yy)^2;

valz=aux1+aux2;
