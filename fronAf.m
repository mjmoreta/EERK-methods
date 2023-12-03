function valz=fronAf(tt,xx,yy)

% Auxiliary function to calculate the boundary value Af(tn,u(tn)) 

valz=-4*sin(tt+xx+yy)^2-4*cos(tt+xx+yy)+2*sin(tt+xx+yy);

