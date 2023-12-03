function valz=fronf(tt,kk,xx,yy)

% Function that gives the boundary value 
% f(tn+k,u(tn+k ut(tn)+k^2/2 utt(tn))

aux1=(cos(tt+xx+yy)-kk*sin(tt+xx+yy)-(kk^2/2)*cos(tt+xx+yy))^2;
aux2=2*cos(tt+kk+xx+yy)-sin(tt+kk+xx+yy)-cos(tt+kk+xx+yy)^2;

valz=aux1+aux2;