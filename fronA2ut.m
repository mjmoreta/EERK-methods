function valz=fronA2ut(tt,xx,yy)

% Auxiliary function to calculate the boundary value A^2u_t(tn,u(tn))

valz=8*sin(tt+xx+yy)*cos(tt+xx+yy)-4*sin(tt+xx+yy);
