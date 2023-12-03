function valz=frona2u(tt,xx,yy)

% Auxiliary function to calculate the boundary value A^2u(tn,u(tn))

valz=4*cos(tt+xx+yy)+4*sin(tt+xx+yy)^2;