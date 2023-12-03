function valz=fronAfu(tt,kk,xx,yy)

% Auxiliary function to calculate the boundary value
% Af(tn+k,u(tn)+kut(tn))

term2=-4*(cos(tt+xx+yy)-kk*sin(tt+xx+yy))^2;
term3=-4*cos(tt+kk+xx+yy)+2*sin(tt+kk+xx+yy)+4-8*sin(tt+kk+xx+yy)^2;

valz=term2+term3;