% Program that implements Krogstad method, given by tableau
%
% 0    |   
% 1/2  |  phi_{1,2}/2
% 1/2  |  phi_{1,3}/2-phi_{2,3}   phi_{2,3}
%  1   |  phi_{1,4}-2phi_{2,4}       0           2phi_{2,4}   
% -------------------------------------------------------------------------
%      |  phi_1-3phi_2+4phi_3   2phi_2-4phi_3   2phi_2-4phi_3   4phi_3-phi2  
%
% avoiding the order reduction with p=3 for problem 
% u_t=u_xx+u_yy+g(t,u), (x,y)\in [0,1]x[0,1] with boundary conditions
% u(t,0,y)=gx0(t,y), u(t,1,y)=gx1(t,y), u(t,x,0)=gy0(t,x), and u(t,x,1)=gy1(t,x) 
% with g(t,u)=u^2+h(t,u) and with exact solution u(t,x,y)=cos(t+x+y)
% For the spatial discretization we have used the 9 point formula. 
% Krogstad method appears in 
% S. Krogstad, Generalized integrated factor methods for stiff PDEs.
% J. Comput. Phys. 203 (2005) 72-88.

% Some initial values, such as h and some more that are used several times 
% along the program
% JJ is such that h=1/JJ is the grid diameter in [0,1] for x and y

JJ=80;
dim=JJ-1;
dim2=(JJ-1)*(JJ-1);
h=1/JJ;

Chtilde=zeros(dim2,1);
Dhtilde=zeros(dim2,1);
F1=zeros(dim2,1);
F2=zeros(dim2,1);
F3=zeros(dim2,1);
F4=zeros(dim2,1);
K2=zeros(dim2,1);
K3=zeros(dim2,1);
K4=zeros(dim2,1);
vech=zeros(dim2,1);
vechk2=zeros(dim2,1);
vechk=zeros(dim2,1);
uxx0=zeros(dim,1);
uxx1=zeros(dim,1);
uy0y=zeros(dim,1);
uy1y=zeros(dim,1);
ux0y=zeros(dim,1);
ux1y=zeros(dim,1);
uyx0=zeros(dim,1);
uyx1=zeros(dim,1);
uxtx0=zeros(dim,1);
uxtx1=zeros(dim,1);
uyt0y=zeros(dim,1);
uyt1y=zeros(dim,1);
uxt0y=zeros(dim,1);
uxt1y=zeros(dim,1);
uytx0=zeros(dim,1);
uytx1=zeros(dim,1);
ERROR1=zeros(dim,dim);
x=zeros(dim,1);
y=zeros(dim,1);
x=[h:h:1]';
y=x;

% Matrix A_{h,0} is given by M^{-1}A. M^{-1}, is not computed, when
% necessary, a system is solved or phipmM9points is used. Matrices A and M 
% are not computed because they are very large. We use their values different 
% from 0 when necessary

% n is such that the time step size is k=1/n.
n=2;
k=1/n;
Tf=n*k;
h2=h^2;
mult=12/h2;
mult2=h2/12;

% The program runs for 9 different values of k, from k=1/2, in order to calculate
% the error and the order of the method
for ll=1:9     
    k2=k/2;
    ksq=k^2;
    kd=2*k;
    ksq2=ksq/2;
    ksqd=2*ksq;
    
    t=0;
     
    % U is the initial value of u(0,x,y)
    U=zeros(dim2,1);
    for ii=1:dim
        for jj=1:dim
            U(ii+(jj-1)*dim,1)=cos(x(ii)+y(jj));
        end
    end
    
   % No se calculan las funciones phi. Al terminar se usa la subrutina
   % phipmM4points.

   % For the local error r=1. For the global one, r=n         
    for kk=1:1

        % Ch and Dh are calculated by solving a system M x=u. Here, u is
        % calculated u and then, system M x=u is solved. 
        % There are two types of boundary values, Ch g and D_h g. 
        % Boundaries of the form C_h g are of the form (12/h^2) M^{-1} bound Ag 
        % and the ones of the form D_h g are M^{-1} fron M g. 
        
        % All the vector are (JJ-1)*(JJ-1). The result corresponding to
        % (x(i),y(JJ)) is at position ii+JJ(jj-1)
        
        % Some numerical differentiation is needed. We refer to the
        % article that is indicated at readme.txt to the explanation        
        
        bb=2/3;
        cc=1/6;   
                
        if kk==1
            U0=U;
        end        
        if kk==2
            U1=U;
        end        
        if kk==3
            U2=U;
        end
               
        % Calculus of the numerical derivatives that are needed
        for ii=1:dim
            % ux0y approximates u_x(tn,0,y), ux1y approximates u_x(tn,1,y),
            % uxx0 approximates u_x(tn,x,0), uxx1 approximates u_x(tn,x,1),
            % uy0y approximates u_y(tn,0,y), uy1y approximates u_y(tn,1,y),
            % uyx0 approximates u_y(tn,x,0), uyx1 approximates u_y(tn,x,1)            
            
            % When t=0, which corresponds to kk=1, all the values of
            % ux0y, ux1y, uyx0 and uyx1 are known as u(0,x,y) is known.
            % The values that are known for every t are those of
            % ux(t,x,0), ux(t,x,1), uy(t,0,y) and uy(t,1,y)
            
            uxx0(ii)=-sin(t+x(ii));
            uxx1(ii)=-sin(t+x(ii)+1);            
            uy0y(ii)=-sin(t+y(ii));
            uy1y(ii)=-sin(t+1+y(ii));
                        
            if kk==1 
                ux0y(ii)=-sin(y(ii));
                ux1y(ii)=-sin(t+1+y(ii));
                uyx0(ii)=-sin(t+x(ii));
                uyx1(ii)=-sin(t+x(ii)+1); 
            end
            
            % For t>0, from the second step, the approximated values of 
            % ux0y, ux1y, uyx0 and uyx1 are calculated by using a 4th order 
            % BDF fórmula. It is progressive when calculating at x=0 or y=0 
            % and regressive if x=1 or y=1
            % o y=1
            
            if kk>1
                sum1=-25*cos(t+y(ii))/12+4*U(1+(ii-1)*dim,1)-3*U(2+(ii-1)*dim,1);
                sum2=4*U(3+(ii-1)*dim,1)/3-U(4+(ii-1)*dim,1)/4;
                ux0y(ii)=(sum1+sum2)/h; 
                
                sum1=25*cos(1+t+y(ii))/12-4*U(dim+(ii-1)*dim,1)+3*U(dim-1+(ii-1)*dim,1);
                sum2=-4*U(dim-2+(ii-1)*dim,1)/3+U(dim-3+(ii-1)*dim,1)/4;
                ux1y(ii)=(sum1+sum2)/h; 
                
                sum1=-25*cos(t+x(ii))/12+4*U(ii,1)-3*U(ii+dim,1);
                sum2=4*U(ii+2*dim,1)/3-U(ii+3*dim,1)/4;
                uyx0(ii)=(sum1+sum2)/h; 
                
                sum1=25*cos(t+x(ii)+1)/12-4*U(ii+(dim-1)*dim,1)+3*U(ii+(dim-2)*dim,1);
                sum2=-4*U(ii+(dim-3)*dim,1)/3+U(ii+(dim-4)*dim,1)/4;
                uyx1(ii)=(sum1+sum2)/h;    
            end          
            
            
            % uxt0y approximates u_xt(tn,0,y), uxt1y approximates u_xt(tn,1,y)  
            % uxtx0 approximates u_xt(tn,x,0), uxtx1 approximates u_xt(tn,x,1)
            % uyt0y approximates u_yt(tn,0,y), uyt1y approximates u_yt(tn,1,y)
            % uytx0 approximates u_yt(tn,x,0), uytx1 approximates u_yt(tn,x,1)            

            % The values of uxt(t,x,0), uxt(t,x,1), uyt(t,0,y) and
            % uyt(t,1,y) are known for all t>0
            
            uxtx0(ii)=-cos(t+x(ii));
            uxtx1(ii)=-cos(t+x(ii)+1);            
            uyt0y(ii)=-cos(t+y(ii));
            uyt1y(ii)=-cos(t+1+y(ii));
            
            % To calculate the rest of the values, a 4th order BDF formula
            % is used. The values of u_t(t,x,y) are ne3eded. 

            % When t=0 (kk=1), those values are known as u_t(0,x,y) is
            % known. For k=2 and k=3, a second order Taylor polynomial is used                  
            
            if kk==1                             
                ut0=-sin(y(ii));
                ut1=-sin(x(1)+y(ii));               
                ut2=-sin(x(2)+y(ii));
                ut3=-sin(x(3)+y(ii));
                ut4=-sin(x(4)+y(ii));        
                uxt0y(ii)=(-25*ut0/12+4*ut1-3*ut2+4*ut3/3-ut4/4)/h;
                  
                ut0=-sin(1+y(ii));
                ut1=-sin(x(dim)+y(ii));               
                ut2=-sin(x(dim-1)+y(ii));
                ut3=-sin(x(dim-2)+y(ii));
                ut4=-sin(x(dim-3)+y(ii));                
                uxt1y(ii)=(25*ut0/12-4*ut1+3*ut2-4*ut3/3+ut4/4)/h;
                
                ut0=-sin(x(ii));
                ut1=-sin(x(ii)+y(1));               
                ut2=-sin(x(ii)+y(2));
                ut3=-sin(x(ii)+y(3));
                ut4=-sin(x(ii)+y(4));                
                uytx0(ii)=(-25*ut0/12+4*ut1-3*ut2+4*ut3/3-ut4/4)/h;                
                          
                ut0=-sin(1+x(ii));
                ut1=-sin(x(ii)+y(dim));               
                ut2=-sin(x(ii)+y(dim-1));
                ut3=-sin(x(ii)+y(dim-2));
                ut4=-sin(x(ii)+y(dim-3));                
                uytx1(ii)=(25*ut0/12-4*ut1+3*ut2-4*ut3/3+ut4/4)/h;                
                             
            end
                
            if kk==2  
                
                % Firstly, the values of u_t(t,x,y) are approximated. Then, 
                % these valuesa are used to approximated the needed ones
                ut0=-sin(t+y(ii));
                ut1=-sin(x(1)+y(ii))-k*cos(x(1)+y(ii))+ksq2*sin(x(1)+y(ii));
                ut2=-sin(x(2)+y(ii))-k*cos(x(2)+y(ii))+ksq2*sin(x(2)+y(ii));
                ut3=-sin(x(3)+y(ii))-k*cos(x(3)+y(ii))+ksq2*sin(x(3)+y(ii));
                ut4=-sin(x(4)+y(ii))-k*cos(x(4)+y(ii))+ksq2*sin(x(4)+y(ii));        
                uxt0y(ii)=(-25*ut0/12+4*ut1-3*ut2+4*ut3/3-ut4/4)/h;
                  
                ut0=-sin(t+1+y(ii));
                ut1=-sin(x(dim)+y(ii))-k*cos(x(dim)+y(ii))+ksq2*sin(x(dim)+y(ii));
                ut2=-sin(x(dim-1)+y(ii))-k*cos(x(dim-1)+y(ii))+ksq2*sin(x(dim-1)+y(ii));
                ut3=-sin(x(dim-2)+y(ii))-k*cos(x(dim-2)+y(ii))+ksq2*sin(x(dim-2)+y(ii));
                ut4=-sin(x(dim-3)+y(ii))-k*cos(x(dim-3)+y(ii))+ksq2*sin(x(dim-3)+y(ii));            
                uxt1y(ii)=(25*ut0/12-4*ut1+3*ut2-4*ut3/3+ut4/4)/h;
                
                ut0=-sin(t+x(ii));
                ut1=-sin(y(1)+x(ii))-k*cos(y(1)+x(ii))+ksq2*sin(y(1)+x(ii));
                ut2=-sin(y(2)+x(ii))-k*cos(y(2)+x(ii))+ksq2*sin(y(2)+x(ii));
                ut3=-sin(y(3)+x(ii))-k*cos(y(3)+x(ii))+ksq2*sin(y(3)+x(ii));
                ut4=-sin(y(4)+x(ii))-k*cos(y(4)+x(ii))+ksq2*sin(y(4)+x(ii));                
                uytx0(ii)=(-25*ut0/12+4*ut1-3*ut2+4*ut3/3-ut4/4)/h;                
                          
                ut0=-sin(t+1+x(ii));
                ut1=-sin(y(dim)+x(ii))-k*cos(y(dim)+x(ii))+ksq2*sin(y(dim)+x(ii));
                ut2=-sin(y(dim-1)+x(ii))-k*cos(y(dim-1)+x(ii))+ksq2*sin(y(dim-1)+x(ii));
                ut3=-sin(y(dim-2)+x(ii))-k*cos(y(dim-2)+x(ii))+ksq2*sin(y(dim-2)+x(ii));
                ut4=-sin(y(dim-3)+x(ii))-k*cos(y(dim-3)+x(ii))+ksq2*sin(y(dim-3)+x(ii));            
                uytx1(ii)=(25*ut0/12-4*ut1+3*ut2-4*ut3/3+ut4/4)/h;                 
                
            end
            
            if kk==3      
               
                % Firstly, the values of u_t(t,x,y) are approximated. Then, 
                % these valuesa are used to approximated the needed ones                
                ut0=-sin(t+y(ii));
                ut1=-sin(x(1)+y(ii))-kd*cos(x(1)+y(ii))+ksqd*sin(x(1)+y(ii));
                ut2=-sin(x(2)+y(ii))-kd*cos(x(2)+y(ii))+ksqd*sin(x(2)+y(ii));
                ut3=-sin(x(3)+y(ii))-kd*cos(x(3)+y(ii))+ksqd*sin(x(3)+y(ii));
                ut4=-sin(x(4)+y(ii))-kd*cos(x(4)+y(ii))+ksqd*sin(x(4)+y(ii));        
                uxt0y(ii)=(-25*ut0/12+4*ut1-3*ut2+4*ut3/3-ut4/4)/h;
                  
                ut0=-sin(t+1+y(ii));
                ut1=-sin(x(dim)+y(ii))-kd*cos(x(dim)+y(ii))+ksqd*sin(x(dim)+y(ii));
                ut2=-sin(x(dim-1)+y(ii))-kd*cos(x(dim-1)+y(ii))+ksqd*sin(x(dim-1)+y(ii));
                ut3=-sin(x(dim-2)+y(ii))-kd*cos(x(dim-2)+y(ii))+ksqd*sin(x(dim-2)+y(ii));
                ut4=-sin(x(dim-3)+y(ii))-kd*cos(x(dim-3)+y(ii))+ksqd*sin(x(dim-3)+y(ii));            
                uxt1y(ii)=(25*ut0/12-4*ut1+3*ut2-4*ut3/3+ut4/4)/h;
                
                ut0=-sin(t+x(ii));
                ut1=-sin(y(1)+x(ii))-kd*cos(y(1)+x(ii))+ksqd*sin(y(1)+x(ii));
                ut2=-sin(y(2)+x(ii))-kd*cos(y(2)+x(ii))+ksqd*sin(y(2)+x(ii));
                ut3=-sin(y(3)+x(ii))-kd*cos(y(3)+x(ii))+ksqd*sin(y(3)+x(ii));
                ut4=-sin(y(4)+x(ii))-kd*cos(y(4)+x(ii))+ksqd*sin(y(4)+x(ii));                
                uytx0(ii)=(-25*ut0/12+4*ut1-3*ut2+4*ut3/3-ut4/4)/h;                
                          
                ut0=-sin(t+1+x(ii));
                ut1=-sin(y(dim)+x(ii))-kd*cos(y(dim)+x(ii))+ksqd*sin(y(dim)+x(ii));
                ut2=-sin(y(dim-1)+x(ii))-kd*cos(y(dim-1)+x(ii))+ksqd*sin(y(dim-1)+x(ii));
                ut3=-sin(y(dim-2)+x(ii))-kd*cos(y(dim-2)+x(ii))+ksqd*sin(y(dim-2)+x(ii));
                ut4=-sin(y(dim-3)+x(ii))-kd*cos(y(dim-3)+x(ii))+ksqd*sin(y(dim-3)+x(ii));            
                uytx1(ii)=(25*ut0/12-4*ut1+3*ut2-4*ut3/3+ut4/4)/h;          
                
               
            end      
            
            if kk>3   
                
                ut0=-sin(t+y(ii));
                % uti is the aproximation to u_t(xi,yj). 
                ut1=(11*U(1+(ii-1)*dim,1)/6-3*U2(1+(ii-1)*dim,1)+3*U1(1+(ii-1)*dim,1)/2-U0(1+(ii-1)*dim,1)/3)/k;
                ut2=(11*U(2+(ii-1)*dim,1)/6-3*U2(2+(ii-1)*dim,1)+3*U1(2+(ii-1)*dim,1)/2-U0(2+(ii-1)*dim,1)/3)/k;
                ut3=(11*U(3+(ii-1)*dim,1)/6-3*U2(3+(ii-1)*dim,1)+3*U1(3+(ii-1)*dim,1)/2-U0(3+(ii-1)*dim,1)/3)/k;
                ut4=(11*U(4+(ii-1)*dim,1)/6-3*U2(4+(ii-1)*dim,1)+3*U1(4+(ii-1)*dim,1)/2-U0(4+(ii-1)*dim,1)/3)/k;                
                uxt0y(ii)=(-25*ut0/12+4*ut1-3*ut2+4*ut3/3-ut4/4)/h;           
     
                ut0=-sin(1+t+y(ii)); 
                ut1=(11*U(dim+(ii-1)*dim,1)/6-3*U2(dim+(ii-1)*dim,1)+3*U1(dim+(ii-1)*dim,1)/2-U0(dim+(ii-1)*dim,1)/3)/k;
                ut2=(11*U(dim-1+(ii-1)*dim,1)/6-3*U2(dim-1+(ii-1)*dim,1)+3*U1(dim-1+(ii-1)*dim,1)/2-U0(dim-1+(ii-1)*dim,1)/3)/k;
                ut3=(11*U(dim-2+(ii-1)*dim,1)/6-3*U2(dim-2+(ii-1)*dim,1)+3*U1(dim-2+(ii-1)*dim,1)/2-U0(dim-2+(ii-1)*dim,1)/3)/k;
                ut4=(11*U(dim-3+(ii-1)*dim,1)/6-3*U2(dim-3+(ii-1)*dim,1)+3*U1(dim-3+(ii-1)*dim,1)/2-U0(dim-3+(ii-1)*dim,1)/3)/k;                
                uxt1y(ii)=(25*ut0/12-4*ut1+3*ut2-4*ut3/3+ut4/4)/h;               
                
                
                ut0=-sin(t+x(ii));
                % uti is the aproximation to u_t(xi,yj).  
                ut1=(11*U(ii,1)/6-3*U2(ii,1)+3*U1(ii,1)/2-U0(ii,1)/3)/k;
                ut2=(11*U(ii+dim,1)/6-3*U2(ii+dim,1)+3*U1(ii+dim,1)/2-U0(ii+dim,1)/3)/k;
                ut3=(11*U(ii+2*dim,1)/6-3*U2(ii+2*dim,1)+3*U1(ii+2*dim,1)/2-U0(ii+2*dim,1)/3)/k;
                ut4=(11*U(ii+3*dim,1)/6-3*U2(ii+3*dim,1)+3*U1(ii+3*dim,1)/2-U0(ii+3*dim,1)/3)/k;                
                uytx0(ii)=(-25*ut0/12+4*ut1-3*ut2+4*ut3/3-ut4/4)/h;               
                
                ut0=-sin(1+t+x(ii));
                % uti is the aproximation to u_t(xi,yj).  
                ut1=(11*U(ii+(dim-1)*dim,1)/6-3*U2(ii+(dim-1)*dim,1)+3*U1(ii+(dim-1)*dim,1)/2-U0(ii+(dim-1)*dim,1)/3)/k;
                ut2=(11*U(ii+(dim-2)*dim,1)/6-3*U2(ii+(dim-2)*dim,1)+3*U1(ii+(dim-2)*dim,1)/2-U0(ii+(dim-2)*dim,1)/3)/k;
                ut3=(11*U(ii+(dim-3)*dim,1)/6-3*U2(ii+(dim-3)*dim,1)+3*U1(ii+(dim-3)*dim,1)/2-U0(ii+(dim-3)*dim,1)/3)/k;
                ut4=(11*U(ii+(dim-4)*dim,1)/6-3*U2(ii+(dim-4)*dim,1)+3*U1(ii+(dim-4)*dim,1)/2-U0(ii+(dim-4)*dim,1)/3)/k;                
                uytx1(ii)=(-25*ut0/12+4*ut1-3*ut2+4*ut3/3-ut4/4)/h;       
                                
            end    
        end
               
        % Once the values of U0, U1 and U2 have been used, they are updated
        if kk>3
            U0=U1;
            U1=U2;
            U2=U;
        end
        
        % Calculus of all the boundary values                 
        MM1=cos(t+y(1))+cos(t+x(1));
        MM2=cos(t+x(2))+cos(t+y(2))+cos(t);
        Chtilde(1,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=cos(t+x(ii));
            MM2=cos(t+x(ii+1))+cos(t+x(ii-1));
            Chtilde(ii,1)=bb*MM1+cc*MM2;
        end
        MM1=cos(t+1+y(1))+cos(t+x(dim));
        MM2=cos(t+1+y(2))+cos(t+1)+cos(t+x(dim-1));
        Chtilde(dim,1)=bb*MM1+cc*MM2;
           
        for m=2:dim-1
            MM1=cos(t+y(m));
            MM2=cos(t+y(m+1))+cos(t+y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;
            
            MM1=cos(t+1+y(m));
            MM2=cos(t+1+y(m+1))+cos(t+1+y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
        end
           
        MM1=cos(t+x(1)+1)+cos(t+y(dim));
        MM2=cos(t+x(2)+1)+cos(t+1)+cos(t+y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=cos(t+x(ii)+1);
            MM2=cos(t+x(ii+1)+1)+cos(t+x(ii-1)+1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        end
        MM1=cos(t+1+y(dim))+cos(t+x(dim)+1);
        MM2=cos(t+2)+cos(t+1+y(dim-1))+cos(t+x(dim-1)+1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        
        Chu=mult*multMm1(dim,Chtilde);       
        
        % Boundaries Chut and D_h Au.          
        MM1=-sin(t+y(1))-sin(t+x(1));
        MM2=-sin(t+x(2))-sin(t+y(2))-sin(t); 
        Dhtilde(1,1)=-2*cos(t+y(1))-2*cos(t+x(1));        
        Chtilde(1,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=-sin(t+x(ii));
            MM2=-sin(t+x(ii+1))-sin(t+x(ii-1));
            Chtilde(ii,1)=bb*MM1+cc*MM2;
            Dhtilde(ii,1)=-2*cos(t+x(ii));
        end
        MM1=-sin(t+1+y(1))-sin(t+x(dim));
        MM2=-sin(t+1+y(2))-sin(t+1)-sin(t+x(dim-1));
        Chtilde(dim,1)=bb*MM1+cc*MM2;
        Dhtilde(dim,1)=-2*cos(t+1+y(1))-2*cos(t+x(dim));              
        for m=2:dim-1
            MM1=-sin(t+y(m));
            MM2=-sin(t+y(m+1))-sin(t+y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;
            Dhtilde(1+(m-1)*dim,1)=-2*cos(t+y(m));            
            MM1=-sin(t+1+y(m));
            MM2=-sin(t+1+y(m+1))-sin(t+1+y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;                        
            Dhtilde(dim+(m-1)*dim,1)=-2*cos(t+1+y(m)); 
        end           
        MM1=-sin(t+x(1)+1)-sin(t+y(dim));
        MM2=-sin(t+x(2)+1)-sin(t+1)-sin(t+y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;        
        Dhtilde(1+(dim-1)*dim,1)=-2*cos(t+y(dim))-2*cos(t+x(1)+1); 
        for ii=2:dim-1
            MM1=-sin(t+x(ii)+1);
            MM2=-sin(t+x(ii+1)+1)-sin(t+x(ii-1)+1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
            Dhtilde(ii+(dim-1)*dim,1)=-2*cos(t+x(ii)+1);  
        end
        MM1=-sin(t+1+y(dim))-sin(t+x(dim)+1);
        MM2=-sin(t+2)-sin(t+1+y(dim-1))-sin(t+x(dim-1)+1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        Dhtilde(dim2,1)=-2*cos(t+1+y(dim))-2*cos(t+x(dim)+1);
         
        Chut=mult*multMm1(dim,Chtilde);  
        DhAu=multMm1(dim,Dhtilde);          
           
        
        % Boundary Chf. It is used to calculate ChAU=Chut-Chf
        MM1=2*(cos(t+y(1))+cos(t+x(1)))-sin(t+y(1))-sin(t+x(1));
        MM2=2*(cos(t+x(2))+cos(t+y(2))+cos(t))-sin(t+x(2))-sin(t+y(2))-sin(t);
        Chtilde(1,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=2*cos(t+x(ii))-sin(t+x(ii));
            MM2=2*(cos(t+x(ii+1))+cos(t+x(ii-1)))-sin(t+x(ii+1))-sin(t+x(ii-1));
            Chtilde(ii,1)=bb*MM1+cc*MM2;
        end
        MM1=2*(cos(t+1+y(1))+cos(t+x(dim)))-sin(t+1+y(1))-sin(t+x(dim));
        MM2=2*(cos(t+1+y(2))+cos(t+1)+cos(t+x(dim-1)))-sin(t+1+y(2))-sin(t+1)-sin(t+x(dim-1));
        Chtilde(dim,1)=bb*MM1+cc*MM2;           
        for m=2:dim-1
            MM1=2*cos(t+y(m))-sin(t+y(m));
            MM2=2*(cos(t+y(m+1))+cos(t+y(m-1)))-sin(t+y(m+1))-sin(t+y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;            
            MM1=2*cos(t+1+y(m))-sin(t+1+y(m));
            MM2=2*(cos(t+1+y(m+1))+cos(t+1+y(m-1)))-sin(t+1+y(m+1))-sin(t+1+y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
        end
           
        MM1=2*(cos(t+x(1)+1)+cos(t+y(dim)))-sin(t+x(1)+1)-sin(t+y(dim));        
        MM2=2*(cos(t+x(2)+1)+cos(t+1)+cos(t+y(dim-1)))-sin(t+x(2)+1)-sin(t+1)-sin(t+y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=2*cos(t+x(ii)+1)-sin(t+x(ii)+1);
            MM2=2*(cos(t+x(ii+1)+1)+cos(t+x(ii-1)+1))-sin(t+x(ii+1)+1)-sin(t+x(ii-1)+1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        end
        MM1=2*(cos(t+1+y(dim))+cos(t+x(dim)+1))-sin(t+1+y(dim))-sin(t+x(dim)+1);
        MM2=2*(cos(t+2)+cos(t+1+y(dim-1))+cos(t+x(dim-1)+1))-sin(t+2)-sin(t+1+y(dim-1))-sin(t+x(dim-1)+1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        
        Chf=mult*multMm1(dim,Chtilde);
          
        
        
        % Boundaries DhAut and ChAut
        MM1=2*sin(t+y(1))+2*sin(t+x(1));
        MM2=2*sin(t+x(2))+2*sin(t+y(2))+2*sin(t);
        Chtilde(1,1)=bb*MM1+cc*MM2;
        Dhtilde(1,1)=MM1;
        for ii=2:dim-1
            MM1=2*sin(t+x(ii));
            MM2=2*sin(t+x(ii+1))+2*sin(t+x(ii-1));
            Chtilde(ii,1)=bb*MM1+cc*MM2;
            Dhtilde(ii,1)=MM1;
        end
        MM1=2*sin(t+1+y(1))+2*sin(t+x(dim));
        MM2=2*sin(t+1+y(2))+2*sin(t+1)+2*sin(t+x(dim-1));
        Chtilde(dim,1)=bb*MM1+cc*MM2;
        Dhtilde(dim,1)=MM1;           
        for m=2:dim-1
            MM1=2*sin(t+y(m));
            MM2=2*sin(t+y(m+1))+2*sin(t+y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;
            Dhtilde(1+(m-1)*dim,1)=MM1;            
            MM1=2*sin(t+1+y(m));
            MM2=2*sin(t+1+y(m+1))+2*sin(t+1+y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
            Dhtilde(dim+(m-1)*dim,1)=MM1;  
        end           
        MM1=2*sin(t+x(1)+1)+2*sin(t+y(dim));
        MM2=2*sin(t+x(2)+1)+2*sin(t+1)+2*sin(t+y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        Dhtilde(1+(dim-1)*dim,1)=MM1;
        for ii=2:dim-1
            MM1=2*sin(t+x(ii)+1);
            MM2=2*sin(t+x(ii+1)+1)+2*sin(t+x(ii-1)+1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
            Dhtilde(ii+(dim-1)*dim,1)=MM1;
        end
        MM1=2*sin(t+1+y(dim))+2*sin(t+x(dim)+1);
        MM2=2*sin(t+2)+2*sin(t+1+y(dim-1))+2*sin(t+x(dim-1)+1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        Dhtilde(dim2,1)=MM1;        
        ChAut=mult*multMm1(dim,Chtilde);         
        DhAut=multMm1(dim,Dhtilde);  
                 
        
        % Boundary Chf(tn+k2,u(tn)+k2ut(tn)) 
        % We use function fronfu in order to simplify the program
        MM1=fronfu(t,k2,0,y(1))+fronfu(t,k2,x(1),0);
        MM2=fronfu(t,k2,x(2),0)+fronfu(t,k2,0,y(2))+fronfu(t,k2,0,0);
        Chtilde(1,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=fronfu(t,k2,x(ii),0);
            MM2=fronfu(t,k2,x(ii+1),0)+fronfu(t,k2,x(ii-1),0);
            Chtilde(ii,1)=bb*MM1+cc*MM2;
        end
        MM1=fronfu(t,k2,1,y(1))+fronfu(t,k2,x(dim),0);
        MM2=fronfu(t,k2,1,y(2))+fronfu(t,k2,1,0)+fronfu(t,k2,x(dim-1),0);
        Chtilde(dim,1)=bb*MM1+cc*MM2;
           
        for m=2:dim-1
            MM1=fronfu(t,k2,0,y(m));
            MM2=fronfu(t,k2,0,y(m+1))+fronfu(t,k2,0,y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;
            
            MM1=fronfu(t,k2,1,y(m));
            MM2=fronfu(t,k2,1,y(m+1))+fronfu(t,k2,1,y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
        end
           
        MM1=fronfu(t,k2,x(1),1)+fronfu(t,k2,0,y(dim));
        MM2=fronfu(t,k2,x(2),1)+fronfu(t,k2,0,1)+fronfu(t,k2,0,y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=fronfu(t,k2,x(ii),1);
            MM2=fronfu(t,k2,x(ii+1),1)+fronfu(t,k2,x(ii-1),1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        end
        MM1=fronfu(t,k2,1,y(dim))+fronfu(t,k2,x(dim),1);
        MM2=fronfu(t,k2,1,1)+fronfu(t,k2,1,y(dim-1))+fronfu(t,k2,x(dim-1),1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        
        Chfk2=mult*multMm1(dim,Chtilde);
               
        
        % Boundaries ChAf and DhAf        
        % We use function funAf for the part that doesn't need numerical
        % approximation        
        a1=2*(ux0y(1)^2+uy0y(1)^2);
        a2=2*(uxx0(1)^2+uyx0(1)^2);
        MM1=a1+fronAf(t,0,y(1))+a2+fronAf(t,x(1),0);       
        a1=2*(uxx0(2)^2+uyx0(2)^2);
        a2=2*(ux0y(2)^2+uy0y(2)^2);
        a3=4*sin(t)^2;
        MM2=a1+fronAf(t,x(2),0)+a2+fronAf(t,0,y(2))+a3+fronAf(t,0,0);          
        Chtilde(1,1)=bb*MM1+cc*MM2;
        Dhtilde(1,1)=MM1;
        for ii=2:dim-1
            a1=2*(uxx0(ii)^2+uyx0(ii)^2);
            a2=2*(uxx0(ii+1)^2+uyx0(ii+1)^2);
            a3=2*(uxx0(ii-1)^2+uyx0(ii-1)^2);
            MM1=a1+fronAf(t,x(ii),0);
            MM2=a2+fronAf(t,x(ii+1),0)+a3+fronAf(t,x(ii-1),0);          
            Chtilde(ii,1)=bb*MM1+cc*MM2;
            Dhtilde(ii,1)=MM1;
        end
        a1=2*(ux1y(1)^2+uy1y(1)^2);
        a2=2*(uxx0(dim)^2+uyx0(dim)^2);
        MM1=a1+fronAf(t,1,y(1))+a2+fronAf(t,x(dim),0);               
        a1=2*(ux1y(2)^2+uy1y(2)^2);
        a2=4*sin(t+1)^2;
        a3=2*(uxx0(dim-1)^2+uyx0(dim-1)^2);
        MM2=a1+fronAf(t,1,y(2))+a2+fronAf(t,1,0)+a3+fronAf(t,x(dim-1),0);       
        Chtilde(dim,1)=bb*MM1+cc*MM2;
        Dhtilde(dim,1)=MM1;           
        for m=2:dim-1
            a1=2*(ux0y(m)^2+uy0y(m)^2);
            a2=2*(ux0y(m+1)^2+uy0y(m+1)^2);
            a3=2*(ux0y(m-1)^2+uy0y(m-1)^2);
            MM1=a1+fronAf(t,0,y(m));
            MM2=a2+fronAf(t,0,y(m+1))+a3+fronAf(t,0,y(m-1));                               
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;
            Dhtilde(1+(m-1)*dim,1)=MM1;            
            a1=2*(ux1y(m)^2+uy1y(m)^2);
            a2=2*(ux1y(m+1)^2+uy1y(m+1)^2);
            a3=2*(ux1y(m-1)^2+uy1y(m-1)^2);            
            MM1=a1+fronAf(t,1,y(m));
            MM2=a2+fronAf(t,1,y(m+1))+a3+fronAf(t,1,y(m-1));            
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
            Dhtilde(dim+(m-1)*dim,1)=MM1;
        end           
        a1=2*(uxx1(1)^2+uyx1(1)^2);
        a2=2*(ux0y(dim)^2+uy0y(dim)^2);        
        MM1=a1+fronAf(t,x(1),1)+a2+fronAf(t,0,y(dim));     
        a1=2*(uxx1(2)^2+uyx1(2)^2);
        a2=4*sin(t+1)^2;
        a3=2*(ux0y(dim-1)^2+uy0y(dim-1)^2);        
        MM2=a1+fronAf(t,x(2),1)+a2+fronAf(t,0,1)+a3+fronAf(t,0,y(dim-1));       
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        Dhtilde(1+(dim-1)*dim,1)=MM1;
        for ii=2:dim-1
            a1=2*(uxx1(ii)^2+uyx1(ii)^2);
            a2=2*(uxx1(ii+1)^2+uyx1(ii+1)^2);
            a3=2*(uxx1(ii-1)^2+uyx1(ii-1)^2);
            MM1=a1+fronAf(t,x(ii),1);
            MM2=a2+fronAf(t,x(ii+1),1)+a3+fronAf(t,x(ii-1),1);                         
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
            Dhtilde(ii+(dim-1)*dim,1)=MM1;
        end
        a1=2*(ux1y(dim)^2+uy1y(dim)^2);
        a2=2*(uxx1(dim)^2+uyx1(dim)^2);
        MM1=a1+fronAf(t,1,y(dim))+a2+fronAf(t,x(dim),1);        
        a1=4*sin(t+2)^2;
        a2=2*(ux1y(dim-1)^2+uy1y(dim-1)^2);
        a3=2*(uxx1(dim-1)^2+uyx1(dim-1)^2);
        MM2=a1+fronAf(t,1,1)+a2+fronAf(t,1,y(dim-1))+a3+fronAf(t,x(dim-1),1);         
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        Dhtilde(dim2,1)=MM1;
        
        ChAf=mult*multMm1(dim,Chtilde);        
        DhAf=multMm1(dim,Dhtilde);  
        
        
        % Boundaries DhA2u and ChA2u        
        a1=-2*(ux0y(1)^2+uy0y(1)^2);
        a2=-2*(uxx0(1)^2+uyx0(1)^2);
        MM1=a1+fronA2u(t,0,y(1))+a2+fronA2u(t,x(1),0);        
        a1=-2*(uxx0(2)^2+uyx0(2)^2);
        a2=-2*(ux0y(2)^2+uy0y(2)^2);
        a3=-4*sin(t)^2;
        MM2=a1+fronA2u(t,x(2),0)+a2+fronA2u(t,0,y(2))+a3+fronA2u(t,0,0);
        Chtilde(1,1)=bb*MM1+cc*MM2;
        Dhtilde(1,1)=MM1;
        for ii=2:dim-1
            a1=2*(uxx0(ii)^2+uyx0(ii)^2);
            a2=2*(uxx0(ii+1)^2+uyx0(ii+1)^2);
            a3=2*(uxx0(ii-1)^2+uyx0(ii-1)^2);           
            MM1=a1+fronA2u(t,x(ii),0);
            MM2=a2+fronA2u(t,x(ii+1),0)+a3+fronA2u(t,x(ii-1),0);
            Chtilde(ii,1)=bb*MM1+cc*MM2;
            Dhtilde(ii,1)=MM1;
        end
        a1=2*(ux1y(1)^2+uy1y(1)^2);
        a2=2*(uxx0(dim)^2+uyx0(dim)^2);
        MM1=a1+fronA2u(t,1,y(1))+a2+fronA2u(t,x(dim),0);               
        a1=2*(ux1y(2)^2+uy1y(2)^2);
        a2=4*sin(t+1)^2;
        a3=2*(uxx0(dim-1)^2+uyx0(dim-1)^2);
        MM2=a1+fronA2u(t,1,y(2))+a2+fronA2u(t,1,0)+a3+fronA2u(t,x(dim-1),0);
        Chtilde(dim,1)=bb*MM1+cc*MM2;
        Dhtilde(dim,1)=MM1;           
        for m=2:dim-1
            a1=2*(ux0y(m)^2+uy0y(m)^2);
            a2=2*(ux0y(m+1)^2+uy0y(m+1)^2);
            a3=2*(ux0y(m-1)^2+uy0y(m-1)^2);             
            MM1=a1+fronA2u(t,0,y(m));
            MM2=a2+fronA2u(t,0,y(m+1))+a3+fronA2u(t,0,y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;
            Dhtilde(1+(m-1)*dim,1)=MM1;            
            a1=2*(ux1y(m)^2+uy1y(m)^2);
            a2=2*(ux1y(m+1)^2+uy1y(m+1)^2);
            a3=2*(ux1y(m-1)^2+uy1y(m-1)^2); 
            MM1=a1+fronA2u(t,1,y(m));
            MM2=a2+fronA2u(t,1,y(m+1))+a3+fronA2u(t,1,y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
            Dhtilde(dim+(m-1)*dim,1)=MM1;
        end        
        a1=2*(uxx1(1)^2+uyx1(1)^2);
        a2=2*(ux0y(dim)^2+uy0y(dim)^2);        
        MM1=a1+fronA2u(t,x(1),1)+a2+fronA2u(t,0,y(dim));     
        a1=2*(uxx1(2)^2+uyx1(2)^2);
        a2=4*sin(t+1)^2;
        a3=2*(ux0y(dim-1)^2+uy0y(dim-1)^2); 
        MM2=a1+fronA2u(t,x(2),1)+a2+fronA2u(t,0,1)+a3+fronA2u(t,0,y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        Dhtilde(1+(dim-1)*dim,1)=MM1;
        for ii=2:dim-1
            a1=2*(uxx1(ii)^2+uyx1(ii)^2);
            a2=2*(uxx1(ii+1)^2+uyx1(ii+1)^2);
            a3=2*(uxx1(ii-1)^2+uyx1(ii-1)^2);  
            MM1=a1+fronA2u(t,x(ii),1);
            MM2=a2+fronA2u(t,x(ii+1),1)+a3+fronA2u(t,x(ii-1),1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
            Dhtilde(ii+(dim-1)*dim,1)=MM1;
        end
        a1=2*(ux1y(dim)^2+uy1y(dim)^2);
        a2=2*(uxx1(dim)^2+uyx1(dim)^2);
        MM1=a1+fronA2u(t,1,y(dim))+a2+fronA2u(t,x(dim),1);          
        a1=4*sin(t+2)^2;
        a2=2*(ux1y(dim-1)^2+uy1y(dim-1)^2);
        a3=2*(uxx1(dim-1)^2+uyx1(dim-1)^2);
        MM2=a1+fronA2u(t,1,1)+a2+fronA2u(t,1,y(dim-1))+a3+fronA2u(t,x(dim-1),1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        Dhtilde(dim2,1)=MM1;
        
        ChA2u=mult*multMm1(dim,Chtilde);        
        DhA2u=multMm1(dim,Dhtilde);  
                
        
        % Boundaries DhA2ut and ChA2ut
        a1=-4*(ux0y(1)*uxt0y(1)+uy0y(1)*uyt0y(1));
        a2=-4*(uxx0(1)*uxtx0(1)+uyx0(1)*uytx0(1));       
        MM1=a1+fronA2ut(t,0,y(1))+a2+fronA2ut(t,x(1),0);        
        a1=-4*(uxx0(2)*uxtx0(2)+uyx0(2)*uytx0(2));
        a2=-4*(ux0y(2)*uxt0y(2)+uy0y(2)*uyt0y(2));
        a3=-8*sin(t)*cos(t);        
        MM2=a1+fronA2ut(t,x(2),0)+a2+fronA2ut(t,0,y(2))+a3+fronA2ut(t,0,0);
        Chtilde(1,1)=bb*MM1+cc*MM2;
        Dhtilde(1,1)=MM1;
        for ii=2:dim-1
            a1=-4*(uxx0(ii)*uxtx0(ii)+uyx0(ii)*uytx0(ii));
            a2=-4*(uxx0(ii+1)*uxtx0(ii+1)+uyx0(ii+1)*uytx0(ii+1));
            a3=-4*(uxx0(ii-1)*uxtx0(ii-1)+uyx0(ii-1)*uytx0(ii-1));            
            MM1=a1+fronA2ut(t,x(ii),0);
            MM2=a2+fronA2ut(t,x(ii+1),0)+a3+fronA2ut(t,x(ii-1),0);
            Chtilde(ii,1)=bb*MM1+cc*MM2;
            Dhtilde(ii,1)=MM1;
        end
        a1=-4*(ux1y(1)*uxt1y(1)+uy1y(1)*uyt1y(1));
        a2=-4*(uxx0(dim)*uxtx0(dim)+uyx0(dim)*uytx0(dim));         
        MM1=a1+fronA2ut(t,1,y(1))+a2+fronA2ut(t,x(dim),0);        
        a1=-4*(ux1y(2)*uxt1y(2)+uy1y(2)*uyt1y(2));
        a2=-8*sin(1+t)*cos(1+t);
        a3=-4*(uxx0(dim-1)*uxtx0(dim-1)+uyx0(dim-1)*uytx0(dim-1));     
        MM2=a1+fronA2ut(t,1,y(2))+a2+fronA2ut(t,1,0)+a3+fronA2ut(t,x(dim-1),0);
        Chtilde(dim,1)=bb*MM1+cc*MM2;
        Dhtilde(dim,1)=MM1;           
        for m=2:dim-1
            a1=-4*(ux0y(m)*uxt0y(m)+uy0y(m)*uyt0y(m));
            a2=-4*(ux0y(m+1)*uxt0y(m+1)+uy0y(m+1)*uyt0y(m+1));
            a3=-4*(ux0y(m-1)*uxt0y(m-1)+uy0y(m-1)*uyt0y(m-1)); 
            MM1=a1+fronA2ut(t,0,y(m));
            MM2=a2+fronA2ut(t,0,y(m+1))+a3+fronA2ut(t,0,y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;
            Dhtilde(1+(m-1)*dim,1)=MM1;            
            a1=-4*(ux1y(m)*uxt1y(m)+uy1y(m)*uyt1y(m));
            a2=-4*(ux1y(m+1)*uxt1y(m+1)+uy1y(m+1)*uyt1y(m+1));
            a3=-4*(ux1y(m-1)*uxt1y(m-1)+uy1y(m-1)*uyt1y(m-1));  
            MM1=a1+fronA2ut(t,1,y(m));
            MM2=a2+fronA2ut(t,1,y(m+1))+a3+fronA2ut(t,1,y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
            Dhtilde(dim+(m-1)*dim,1)=MM1;
        end           
        a1=-4*(uxx1(1)*uxtx1(1)+uyx1(1)*uytx1(1));
        a2=-4*(ux0y(dim)*uxt0y(dim)+uy0y(dim)*uyt0y(dim));        
        MM1=a1+fronA2ut(t,x(1),1)+a2+fronA2ut(t,0,y(dim));        
        a1=-4*(uxx1(2)*uxtx1(2)+uyx1(2)*uytx1(2));
        a2=-8*sin(t+1)*cos(t+1);
        a3=-4*(ux0y(dim-1)*uxt0y(dim-1)+uy0y(dim-1)*uyt0y(dim-1));            
        MM2=a1+fronA2ut(t,x(2),1)+a2+fronA2ut(t,0,1)+a3+fronA2ut(t,0,y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        Dhtilde(1+(dim-1)*dim,1)=MM1;
        for ii=2:dim-1
            a1=-4*(uxx1(ii)*uxtx1(ii)+uyx1(ii)*uytx1(ii));
            a2=-4*(uxx1(ii+1)*uxtx1(ii+1)+uyx1(ii+1)*uytx1(ii+1));            
            a3=-4*(uxx1(ii-1)*uxtx1(ii-1)+uyx1(ii-1)*uytx1(ii-1));
            MM1=a1+fronA2ut(t,x(ii),1);
            MM2=a2+fronA2ut(t,x(ii+1),1)+a3+fronA2ut(t,x(ii-1),1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
            Dhtilde(ii+(dim-1)*dim,1)=MM1;
        end        
        a1=-4*(ux1y(dim)*uxt1y(dim)+uy1y(dim)*uyt1y(dim));
        a2=-4*(uxx1(dim)*uxtx1(dim)+uyx1(dim)*uytx1(dim));
        MM1=a1+fronA2ut(t,1,y(dim))+a2+fronA2ut(t,x(dim),1);               
        a1=-8*sin(t+2)*cos(t+2);
        a2=-4*(ux1y(dim-1)*uxt1y(dim-1)+uy1y(dim-1)*uyt1y(dim-1));
        a3=-4*(uxx1(dim-1)*uxtx1(dim-1)+uyx1(dim-1)*uytx1(dim-1));
        MM2=a1+fronA2ut(t,1,1)+a2+fronA2ut(t,1,y(dim-1))+a3+fronA2ut(t,x(dim-1),1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        Dhtilde(dim2,1)=MM1;
        
        ChA2ut=mult*multMm1(dim,Chtilde);        
        DhA2ut=multMm1(dim,Dhtilde);        
        
        
        % Boundaries DhAf(tn+k2,u(tn)+k2 ut(tn)) and ChA(tn+k2,u(tn)+k2 ut(tn))
        a1=2*(ux0y(1)+k2*uxt0y(1))^2+2*(uy0y(1)+k2*uyt0y(1))^2;
        a2=2*(uxx0(1)+k2*uxtx0(1))^2+2*(uyx0(1)+k2*uytx0(1))^2;        
        MM1=a1+fronAfu(t,k2,0,y(1))+a2+fronAfu(t,k2,x(1),0);        
        a1=2*(uxx0(2)+k2*uxtx0(2))^2+2*(uyx0(2)+k2*uytx0(2))^2; 
        a2=2*(ux0y(2)+k2*uxt0y(2))^2+2*(uy0y(2)+k2*uyt0y(2))^2;
        a3=4*(-sin(t)-k2*cos(t))^2;
        MM2=a1+fronAfu(t,k2,x(2),0)+a2+fronAfu(t,k2,0,y(2))+a3+fronAfu(t,k2,0,0);
        Chtilde(1,1)=bb*MM1+cc*MM2;
        Dhtilde(1,1)=MM1;
        for ii=2:dim-1
            a1=2*(uxx0(ii)+k2*uxtx0(ii))^2+2*(uyx0(ii)+k2*uytx0(ii))^2;
            a2=2*(uxx0(ii+1)+k2*uxtx0(ii+1))^2+2*(uyx0(ii+1)+k2*uytx0(ii+1))^2;
            a3=2*(uxx0(ii-1)+k2*uxtx0(ii-1))^2+2*(uyx0(ii-1)+k2*uytx0(ii-1))^2;
            MM1=a1+fronAfu(t,k2,x(ii),0);
            MM2=a2+fronAfu(t,k2,x(ii+1),0)+a3+fronAfu(t,k2,x(ii-1),0);
            Chtilde(ii,1)=bb*MM1+cc*MM2;
            Dhtilde(ii,1)=MM1;
        end
        a1=2*(ux1y(1)+k2*uxt1y(1))^2+2*(uy1y(1)+k2*uyt1y(1))^2;
        a2=2*(uxx0(dim)+k2*uxtx0(dim))^2+2*(uyx0(dim)+k2*uytx0(dim))^2;
        MM1=a1+fronAfu(t,k2,1,y(1))+a2+fronAfu(t,k2,x(dim),0);        
        a1=2*(ux1y(2)+k2*uxt1y(2))^2+2*(uy1y(2)+k2*uyt1y(2))^2;
        a2=4*(-sin(1+t)-k2*cos(1+t))^2;
        a3=2*(uxx0(dim-1)+k2*uxtx0(dim-1))^2+2*(uyx0(dim-1)+k2*uytx0(dim-1))^2; 
        MM2=a1+fronAfu(t,k2,1,y(2))+a2+fronAfu(t,k2,1,0)+a3+fronAfu(t,k2,x(dim-1),0);
        Chtilde(dim,1)=bb*MM1+cc*MM2;
        Dhtilde(dim,1)=MM1;           
        for m=2:dim-1
            a1=2*(ux0y(m)+k2*uxt0y(m))^2+2*(uy0y(m)+k2*uyt0y(m))^2;
            a2=2*(ux0y(m+1)+k2*uxt0y(m+1))^2+2*(uy0y(m+1)+k2*uyt0y(m+1))^2;
            a3=2*(ux0y(m-1)+k2*uxt0y(m-1))^2+2*(uy0y(m-1)+k2*uyt0y(m-1))^2;
            MM1=a1+fronAfu(t,k2,0,y(m));
            MM2=a2+fronAfu(t,k2,0,y(m+1))+a3+fronAfu(t,k2,0,y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;
            Dhtilde(1+(m-1)*dim,1)=MM1;            
            a1=2*(ux1y(m)+k2*uxt1y(m))^2+2*(uy1y(m)+k2*uyt1y(m))^2;
            a2=2*(ux1y(m+1)+k2*uxt1y(m+1))^2+2*(uy1y(m+1)+k2*uyt1y(m+1))^2;
            a3=2*(ux1y(m-1)+k2*uxt1y(m-1))^2+2*(uy1y(m-1)+k2*uyt1y(m-1))^2;
            MM1=a1+fronAfu(t,k2,1,y(m));
            MM2=a2+fronAfu(t,k2,1,y(m+1))+a3+fronAfu(t,k2,1,y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
            Dhtilde(dim+(m-1)*dim,1)=MM1;
        end           
        a1=2*(uxx1(1)+k2*uxtx1(1))^2+2*(uyx1(1)+k2*uytx1(1))^2; 
        a2=2*(ux0y(dim)+k2*uxt0y(dim))^2+2*(uy0y(dim)+k2*uyt0y(dim))^2;
        MM1=a1+fronAfu(t,k2,x(1),1)+a2+fronAfu(t,k2,0,y(dim));        
        a1=2*(uxx1(2)+k2*uxtx1(2))^2+2*(uyx1(2)+k2*uytx1(2))^2; 
        a2=4*(-sin(t+1)-k2*cos(t+1))^2;        
        a3=2*(ux0y(dim-1)+k2*uxt0y(dim-1))^2+2*(uy0y(dim-1)+k2*uyt0y(dim-1))^2;        
        MM2=a1+fronAfu(t,k2,x(2),1)+a2+fronAfu(t,k2,0,1)+a3+fronAfu(t,k2,0,y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        Dhtilde(1+(dim-1)*dim,1)=MM1;
        for ii=2:dim-1
            a1=2*(uxx1(ii)+k2*uxtx1(ii))^2+2*(uyx1(ii)+k2*uytx1(ii))^2;
            a2=2*(uxx1(ii+1)+k2*uxtx1(ii+1))^2+2*(uyx1(ii+1)+k2*uytx1(ii+1))^2;
            a3=2*(uxx1(ii-1)+k2*uxtx1(ii-1))^2+2*(uyx1(ii-1)+k2*uytx1(ii-1))^2;
            MM1=a1+fronAfu(t,k2,x(ii),1);
            MM2=a2+fronAfu(t,k2,x(ii+1),1)+a3+fronAfu(t,k2,x(ii-1),1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
            Dhtilde(ii+(dim-1)*dim,1)=MM1;
        end
        a1=2*(ux1y(dim)+k2*uxt1y(dim))^2+2*(uy1y(dim)+k2*uyt1y(dim))^2;
        a2=2*(uxx1(dim)+k2*uxtx1(dim))^2+2*(uyx1(dim)+k2*uytx1(dim))^2;
        MM1=a1+fronAfu(t,k2,1,y(dim))+a2+fronAfu(t,k2,x(dim),1);        
        a1=4*(-sin(2+t)-k2*cos(2+t))^2;
        a2=2*(ux1y(dim-1)+k2*uxt1y(dim-1))^2+2*(uy1y(dim-1)+k2*uyt1y(dim-1))^2;
        a3=2*(uxx1(dim-1)+k2*uxtx1(dim-1))^2+2*(uyx1(dim-1)+k2*uytx1(dim-1))^2; 
        MM2=a1+fronAfu(t,k2,1,1)+a2+fronAfu(t,k2,1,y(dim-1))+a3+fronAfu(t,k2,x(dim-1),1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        Dhtilde(dim2,1)=MM1;
        
        ChAfk2=mult*multMm1(dim,Chtilde);        
        DhAfk2=multMm1(dim,Dhtilde);                            
        
        % Boundaries DhAf(tn+k,u(tn)+k ut(tn)) and ChA(tn+k,u(tn)+k ut(tn))
        a1=2*(ux0y(1)+k*uxt0y(1))^2+2*(uy0y(1)+k*uyt0y(1))^2;
        a2=2*(uxx0(1)+k*uxtx0(1))^2+2*(uyx0(1)+k*uytx0(1))^2;        
        MM1=a1+fronAfu(t,k,0,y(1))+a2+fronAfu(t,k,x(1),0);        
        a1=2*(uxx0(2)+k*uxtx0(2))^2+2*(uyx0(2)+k*uytx0(2))^2; 
        a2=2*(ux0y(2)+k*uxt0y(2))^2+2*(uy0y(2)+k*uyt0y(2))^2;
        a3=4*(-sin(t)-k*cos(t))^2;
        MM2=a1+fronAfu(t,k,x(2),0)+a2+fronAfu(t,k,0,y(2))+a3+fronAfu(t,k,0,0);
        Chtilde(1,1)=bb*MM1+cc*MM2;
        Dhtilde(1,1)=MM1;
        for ii=2:dim-1
            a1=2*(uxx0(ii)+k*uxtx0(ii))^2+2*(uyx0(ii)+k*uytx0(ii))^2;
            a2=2*(uxx0(ii+1)+k*uxtx0(ii+1))^2+2*(uyx0(ii+1)+k*uytx0(ii+1))^2;
            a3=2*(uxx0(ii-1)+k*uxtx0(ii-1))^2+2*(uyx0(ii-1)+k*uytx0(ii-1))^2;
            MM1=a1+fronAfu(t,k,x(ii),0);
            MM2=a2+fronAfu(t,k,x(ii+1),0)+a3+fronAfu(t,k,x(ii-1),0);
            Chtilde(ii,1)=bb*MM1+cc*MM2;
            Dhtilde(ii,1)=MM1;
        end
        a1=2*(ux1y(1)+k*uxt1y(1))^2+2*(uy1y(1)+k*uyt1y(1))^2;
        a2=2*(uxx0(dim)+k*uxtx0(dim))^2+2*(uyx0(dim)+k*uytx0(dim))^2;
        MM1=a1+fronAfu(t,k,1,y(1))+a2+fronAfu(t,k,x(dim),0);        
        a1=2*(ux1y(2)+k*uxt1y(2))^2+2*(uy1y(2)+k*uyt1y(2))^2;
        a2=4*(-sin(1+t)-k*cos(1+t))^2;
        a3=2*(uxx0(dim-1)+k*uxtx0(dim-1))^2+2*(uyx0(dim-1)+k*uytx0(dim-1))^2; 
        MM2=a1+fronAfu(t,k,1,y(2))+a2+fronAfu(t,k,1,0)+a3+fronAfu(t,k,x(dim-1),0);
        Chtilde(dim,1)=bb*MM1+cc*MM2;
        Dhtilde(dim,1)=MM1;           
        for m=2:dim-1
            a1=2*(ux0y(m)+k*uxt0y(m))^2+2*(uy0y(m)+k*uyt0y(m))^2;
            a2=2*(ux0y(m+1)+k*uxt0y(m+1))^2+2*(uy0y(m+1)+k*uyt0y(m+1))^2;
            a3=2*(ux0y(m-1)+k*uxt0y(m-1))^2+2*(uy0y(m-1)+k*uyt0y(m-1))^2;
            MM1=a1+fronAfu(t,k,0,y(m));
            MM2=a2+fronAfu(t,k,0,y(m+1))+a3+fronAfu(t,k,0,y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;
            Dhtilde(1+(m-1)*dim,1)=MM1;            
            a1=2*(ux1y(m)+k*uxt1y(m))^2+2*(uy1y(m)+k*uyt1y(m))^2;
            a2=2*(ux1y(m+1)+k*uxt1y(m+1))^2+2*(uy1y(m+1)+k*uyt1y(m+1))^2;
            a3=2*(ux1y(m-1)+k*uxt1y(m-1))^2+2*(uy1y(m-1)+k*uyt1y(m-1))^2;
            MM1=a1+fronAfu(t,k,1,y(m));
            MM2=a2+fronAfu(t,k,1,y(m+1))+a3+fronAfu(t,k,1,y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
            Dhtilde(dim+(m-1)*dim,1)=MM1;
        end           
        a1=2*(uxx1(1)+k*uxtx1(1))^2+2*(uyx1(1)+k*uytx1(1))^2; 
        a2=2*(ux0y(dim)+k*uxt0y(dim))^2+2*(uy0y(dim)+k*uyt0y(dim))^2;
        MM1=a1+fronAfu(t,k,x(1),1)+a2+fronAfu(t,k,0,y(dim));        
        a1=2*(uxx1(2)+k*uxtx1(2))^2+2*(uyx1(2)+k*uytx1(2))^2; 
        a2=4*(-sin(t+1)-k*cos(t+1))^2;        
        a3=2*(ux0y(dim-1)+k*uxt0y(dim-1))^2+2*(uy0y(dim-1)+k*uyt0y(dim-1))^2;        
        MM2=a1+fronAfu(t,k,x(2),1)+a2+fronAfu(t,k,0,1)+a3+fronAfu(t,k,0,y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        Dhtilde(1+(dim-1)*dim,1)=MM1;
        for ii=2:dim-1
            a1=2*(uxx1(ii)+k*uxtx1(ii))^2+2*(uyx1(ii)+k*uytx1(ii))^2;
            a2=2*(uxx1(ii+1)+k*uxtx1(ii+1))^2+2*(uyx1(ii+1)+k*uytx1(ii+1))^2;
            a3=2*(uxx1(ii-1)+k*uxtx1(ii-1))^2+2*(uyx1(ii-1)+k*uytx1(ii-1))^2;
            MM1=a1+fronAfu(t,k,x(ii),1);
            MM2=a2+fronAfu(t,k,x(ii+1),1)+a3+fronAfu(t,k,x(ii-1),1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
            Dhtilde(ii+(dim-1)*dim,1)=MM1;
        end
        a1=2*(ux1y(dim)+k*uxt1y(dim))^2+2*(uy1y(dim)+k*uyt1y(dim))^2;
        a2=2*(uxx1(dim)+k*uxtx1(dim))^2+2*(uyx1(dim)+k*uytx1(dim))^2;
        MM1=a1+fronAfu(t,k,1,y(dim))+a2+fronAfu(t,k,x(dim),1);        
        a1=4*(-sin(2+t)-k*cos(2+t))^2;
        a2=2*(ux1y(dim-1)+k*uxt1y(dim-1))^2+2*(uy1y(dim-1)+k*uyt1y(dim-1))^2;
        a3=2*(uxx1(dim-1)+k*uxtx1(dim-1))^2+2*(uyx1(dim-1)+k*uytx1(dim-1))^2; 
        MM2=a1+fronAfu(t,k,1,1)+a2+fronAfu(t,k,1,y(dim-1))+a3+fronAfu(t,k,x(dim-1),1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        Dhtilde(dim2,1)=MM1;
        
        ChAfk=mult*multMm1(dim,Chtilde);        
        DhAfk=multMm1(dim,Dhtilde);                         
         
        % Boundary Chf(tn+k2,u(tn)+k2ut(tn)+k^2/8 Aut(tn))
        MM1=fronfu2(t,k2,0,y(1))+fronfu2(t,k2,x(1),0);
        MM2=fronfu2(t,k2,x(2),0)+fronfu2(t,k2,0,y(2))+fronfu2(t,k2,0,0);
        Chtilde(1,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=fronfu2(t,k2,x(ii),0);
            MM2=fronfu2(t,k2,x(ii+1),0)+fronfu2(t,k2,x(ii-1),0);
            Chtilde(ii,1)=bb*MM1+cc*MM2;
        end
        MM1=fronfu2(t,k2,1,y(1))+fronfu2(t,k2,x(dim),0);
        MM2=fronfu2(t,k2,1,y(2))+fronfu2(t,k2,1,0)+fronfu2(t,k2,x(dim-1),0);
        Chtilde(dim,1)=bb*MM1+cc*MM2;           
        for m=2:dim-1
            MM1=fronfu2(t,k2,0,y(m));
            MM2=fronfu2(t,k2,0,y(m+1))+fronfu2(t,k2,0,y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;            
            MM1=fronfu2(t,k2,1,y(m));
            MM2=fronfu2(t,k2,1,y(m+1))+fronfu2(t,k2,1,y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
        end           
        MM1=fronfu2(t,k2,x(1),1)+fronfu2(t,k2,0,y(dim));
        MM2=fronfu2(t,k2,x(2),1)+fronfu2(t,k2,0,1)+fronfu2(t,k2,0,y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=fronfu2(t,k2,x(ii),1);
            MM2=fronfu2(t,k2,x(ii+1),1)+fronfu2(t,k2,x(ii-1),1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        end
        MM1=fronfu2(t,k2,1,y(dim))+fronfu2(t,k2,x(dim),1);
        MM2=fronfu2(t,k2,1,1)+fronfu2(t,k2,1,y(dim-1))+fronfu2(t,k2,x(dim-1),1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        
        Chfk2_11=mult*multMm1(dim,Chtilde);        
               
        % Boundary Chf(tn+k,u(tn)+k ut(tn)+k^2/2 utt(tn))
        MM1=fronf(t,k,0,y(1))+fronf(t,k,x(1),0);
        MM2=fronf(t,k,x(2),0)+fronf(t,k,0,y(2))+fronf(t,k,0,0);
        Chtilde(1,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=fronf(t,k,x(ii),0);
            MM2=fronf(t,k,x(ii+1),0)+fronf(t,k,x(ii-1),0);
            Chtilde(ii,1)=bb*MM1+cc*MM2;
        end
        MM1=fronf(t,k,1,y(1))+fronf(t,k,x(dim),0);
        MM2=fronf(t,k,1,y(2))+fronf(t,k,1,0)+fronf(t,k,x(dim-1),0);
        Chtilde(dim,1)=bb*MM1+cc*MM2;           
        for m=2:dim-1
            MM1=fronf(t,k,0,y(m));
            MM2=fronf(t,k,0,y(m+1))+fronf(t,k,0,y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;            
            MM1=fronf(t,k,1,y(m));
            MM2=fronf(t,k,1,y(m+1))+fronf(t,k,1,y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
        end           
        MM1=fronf(t,k,x(1),1)+fronf(t,k,0,y(dim));
        MM2=fronf(t,k,x(2),1)+fronf(t,k,0,1)+fronf(t,k,0,y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=fronf(t,k,x(ii),1);
            MM2=fronf(t,k,x(ii+1),1)+fronf(t,k,x(ii-1),1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        end
        MM1=fronf(t,k,1,y(dim))+fronf(t,k,x(dim),1);
        MM2=fronf(t,k,1,1)+fronf(t,k,1,y(dim-1))+fronf(t,k,x(dim-1),1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        
        Chfk2_12=mult*multMm1(dim,Chtilde);           
        
        % Boundary Chf(tn+k2,u(tn)+k2 ut(tn)+k^2/8 (utt(tn)+ht(tn)...)
        MM1=fronf2(t,k2,0,y(1))+fronf2(t,k2,x(1),0);
        MM2=fronf2(t,k2,x(2),0)+fronf2(t,k2,0,y(2))+fronf2(t,k2,0,0);
        Chtilde(1,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=fronf2(t,k2,x(ii),0);
            MM2=fronf2(t,k2,x(ii+1),0)+fronf2(t,k2,x(ii-1),0);
            Chtilde(ii,1)=bb*MM1+cc*MM2;
        end
        MM1=fronf2(t,k2,1,y(1))+fronf2(t,k2,x(dim),0);
        MM2=fronf2(t,k2,1,y(2))+fronf2(t,k2,1,0)+fronf2(t,k2,x(dim-1),0);
        Chtilde(dim,1)=bb*MM1+cc*MM2;
           
        for m=2:dim-1
            MM1=fronf2(t,k2,0,y(m));
            MM2=fronf2(t,k2,0,y(m+1))+fronf2(t,k2,0,y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;
            
            MM1=fronf2(t,k2,1,y(m));
            MM2=fronf2(t,k2,1,y(m+1))+fronf2(t,k2,1,y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
        end
           
        MM1=fronf2(t,k2,x(1),1)+fronf2(t,k2,0,y(dim));
        MM2=fronf2(t,k2,x(2),1)+fronf2(t,k2,0,1)+fronf2(t,k2,0,y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=fronf2(t,k2,x(ii),1);
            MM2=fronf2(t,k2,x(ii+1),1)+fronf2(t,k2,x(ii-1),1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        end
        MM1=fronf2(t,k2,1,y(dim))+fronf2(t,k2,x(dim),1);
        MM2=fronf2(t,k2,1,1)+fronf2(t,k2,1,y(dim-1))+fronf2(t,k2,x(dim-1),1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        
        Chfk2_13=mult*multMm1(dim,Chtilde);
        
        ChAu=Chut-Chf;
        
        % Evaluations of function h
        for ii=1:dim
            for jj=1:dim
                pos=ii+(jj-1)*dim;
                point=t+x(ii)+y(jj);              
                vech(pos,1)=-sin(point)+2*cos(point)-cos(point)^2;
                vechk2(pos,1)=-sin(point+k2)+2*cos(point+k2)-cos(point+k2)^2;
                vechk(pos,1)=-sin(point+k)+2*cos(point+k)-cos(point+k)^2; 
                
            end
        end     
          
        % Auxiliary vectors that are used
        F1=U.^2+vech;              
        vecb4=zeros(dim2,4);           
        vecb4(:,1)=U;
        vecb4(:,2)=F1+Chu-DhAu;
        vecb4(:,3)=Chut-DhAut;
        vecb4(:,4)=ChAut;
               
        % Stage K2       
        K2=phipmM9points(k2,JJ,vecb4,10^(-13),1,1);                 
                     
        % Stage K3        
        F2=K2.^2+vechk2;         
        vecb4(:,3)=(4*F2-4*F1)/k+Chut-DhAut;
        vecb4(:,4)=(-4*Chf+4*Chfk2)/k+ChAut;
        K3=phipmM9points(k2,JJ,vecb4,10^(-13),1,1);     
                        
        % Stage K4           
        F3=K3.^2+vechk2;                
        vecb4(:,3)=(-2*F1+2*F3)/k+Chut-DhAut;
        vecb4(:,4)=(-2*Chf+2*Chfk2)/k+ChAut;                 
        K4=phipmM9points(k,JJ,vecb4,10^(-13),1,1);
        
        
        % Approximation to the exact solution at time t_{n+1}        
        F4=K4.^2+vechk;             
        vecbdef=zeros(dim2,6);        
        vecbdef(:,1)=U;
        vecbdef(:,2)=F1+Chu-DhAu;        
        vecbdef(:,3)=(-3*F1+2*F2+2*F3-F4)/k+Chut-DhAut;        
        aux1=4*F1-4*F2-4*F3+4*F4;
        aux2=-3*Chf+2*Chfk2_11+2*Chfk2_13-Chfk2_12+3*DhAf-4*DhAfk2+DhAfk;
        aux3=ChAut-DhA2ut;
        vecbdef(:,4)=aux1/k^2+aux2/k+aux3;        
        aux1=4*Chf-4*Chfk2_11-4*Chfk2_13+4*Chfk2_12-4*DhAf+8*DhAfk2-4*DhAfk;
        aux2=-3*ChAf+4*ChAfk2-ChAfk;
        aux3=ChA2ut;
        vecbdef(:,5)=aux1/k^2+aux2/k+aux3;        
        vecbdef(:,6)=(4*ChAf-8*ChAfk2+4*ChAfk)/k^2;      
        
        U=phipmM9points(k,JJ,vecbdef,10^(-13),1,1);    
         
        % New value of t                                            
        t=t+k;  
        
    end
   
    % ERROR1 contains the exact solution at time T  
    for mm=1:dim
        for ii=1:dim
            ERROR1(ii,mm)=cos(t+x(ii)+y(mm))-U(ii+(mm-1)*dim,1);
        end
    end    
  
    % Error in the infinite norm
    errdef1=norm(max(abs(ERROR1)),inf);
    % The order in the infinite norm is calculated as log2(e0/errdef1), 
    % with e0 the error obtained with k and errdef1 the error obtained with k/2. 
    % When ll=1, as there is not a previous error, it can't  be calculated. 
    % Two consecutive errors are compared
    if ll==1
        errdef1
        e0=errdef1;
    else
        [errdef1 log2(e0/errdef1)]
        e0=errdef1;   
    end

    % The new value of n is 2*n, and the new value of k_n=k/2.
    k=k/2;
    n=2*n;      
end
