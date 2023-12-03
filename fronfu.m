function valz=fronfu(tt,kk,xx,yy)

% It is used to calculate the boundary value f(tn+k/2,u(tn)+k/2 ut(tn))

aux1=(cos(tt+xx+yy)-kk*sin(tt+xx+yy))^2;
aux2=2*cos(tt+kk+xx+yy)-sin(tt+kk+xx+yy)-cos(tt+kk+xx+yy)^2;

valz=aux1+aux2;