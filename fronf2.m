function valz=fronf2(tt,kk,xx,yy)

% Function that gives the boundary value 
% f(tn+k2,u(tn)+k2 ut(tn)+k2^2/4 (utt(tn)+h_t(tn)+f_u(tn,u(tn))u_t(tn)))

aux1=(cos(tt+xx+yy)-kk*sin(tt+xx+yy)+(kk^2)*(-cos(tt+xx+yy)-sin(tt+xx+yy)))^2;
aux2=2*cos(tt+kk+xx+yy)-sin(tt+kk+xx+yy)-cos(tt+kk+xx+yy)^2;

valz=aux1+aux2;