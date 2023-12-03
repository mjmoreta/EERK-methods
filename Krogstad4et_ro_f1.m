% Program that implements Krogstad method, given by tableau
%
% 0    |   
% 1/2  |  phi_{1,2}/2
% 1/2  |  phi_{1,3}/2-phi_{2,3}   phi_{2,3}
%  1   |  phi_{1,4}-2phi_{2,4}       0           2phi_{2,4}   
% -------------------------------------------------------------------------
%      |  phi_1-3phi_2+4phi_3   2phi_2-4phi_3   2phi_2-4phi_3   4phi_3-phi2  
%
% without avoiding the order reduction for problem 
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
    
    k=1/n;
    k2=k/2;
    t=0;
            
    % U is the initial value of u(0,x,y)
    U=zeros(dim2,1);
    for ii=1:dim
        for jj=1:dim
            U(ii+(jj-1)*dim,1)=cos(x(ii)+y(jj));
        end
    end
          
    % For the local error r=1. For the global one, r=n      
    for kk=1:1
 	          
        % Frontera Chu. Se calcula como (12/h^2)M^{-1} Chtilde \parcial
        % u.
        
        bb=2/3;
        cc=1/6; 
        
        % Ch and Dh are calculated by solving a system M x=u. Here, u is
        % calculated u and then, system M x=u is solved. 
        % There are two types of boundary values, Ch g and D_h g. 
        % Boundaries of the form C_h g are of the form (12/h^2) M^{-1} bound Ag 
        % and the ones of the form D_h g are M^{-1} fron M g. 
        
        % All the vector are (JJ-1)*(JJ-1). THe result corresponding to
        % (x(i),y(JJ)) is at position ii+JJ(jj-1)
          
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
        
        % Boundary Chu(tn+k2)
        MM1=cos(t+k2+y(1))+cos(t+k2+x(1));
        MM2=cos(t+k2+x(2))+cos(t+k2+y(2))+cos(t+k2);
        Chtilde(1,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=cos(t+k2+x(ii));
            MM2=cos(t+k2+x(ii+1))+cos(t+k2+x(ii-1));
            Chtilde(ii,1)=bb*MM1+cc*MM2;
        end
        MM1=cos(t+k2+1+y(1))+cos(t+k2+x(dim));
        MM2=cos(t+k2+1+y(2))+cos(t+k2+1)+cos(t+k2+x(dim-1));
        Chtilde(dim,1)=bb*MM1+cc*MM2;           
        for m=2:dim-1
            MM1=cos(t+k2+y(m));
            MM2=cos(t+k2+y(m+1))+cos(t+k2+y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;
            
            MM1=cos(t+k2+1+y(m));
            MM2=cos(t+k2+1+y(m+1))+cos(t+k2+1+y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
        end           
        MM1=cos(t+k2+x(1)+1)+cos(t+k2+y(dim));
        MM2=cos(t+k2+x(2)+1)+cos(t+k2+1)+cos(t+k2+y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=cos(t+k2+x(ii)+1);
            MM2=cos(t+k2+x(ii+1)+1)+cos(t+k2+x(ii-1)+1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        end
        MM1=cos(t+k2+1+y(dim))+cos(t+k2+x(dim)+1);
        MM2=cos(t+k2+2)+cos(t+k2+1+y(dim-1))+cos(t+k2+x(dim-1)+1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        
        Chuk2=mult*multMm1(dim,Chtilde); 
                
        % Boundary Chu(tn+k)
        MM1=cos(t+k+y(1))+cos(t+k+x(1));
        MM2=cos(t+k+x(2))+cos(t+k+y(2))+cos(t+k);
        Chtilde(1,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=cos(t+k+x(ii));
            MM2=cos(t+k+x(ii+1))+cos(t+k+x(ii-1));
            Chtilde(ii,1)=bb*MM1+cc*MM2;
        end
        MM1=cos(t+k+1+y(1))+cos(t+k+x(dim));
        MM2=cos(t+k+1+y(2))+cos(t+k+1)+cos(t+k+x(dim-1));
        Chtilde(dim,1)=bb*MM1+cc*MM2;           
        for m=2:dim-1
            MM1=cos(t+k+y(m));
            MM2=cos(t+k+y(m+1))+cos(t+k+y(m-1));
            Chtilde(1+(m-1)*dim,1)=bb*MM1+cc*MM2;
            
            MM1=cos(t+k+1+y(m));
            MM2=cos(t+k+1+y(m+1))+cos(t+k+1+y(m-1));
            Chtilde(dim+(m-1)*dim,1)=bb*MM1+cc*MM2;
        end           
        MM1=cos(t+k+x(1)+1)+cos(t+k+y(dim));
        MM2=cos(t+k+x(2)+1)+cos(t+k+1)+cos(t+k+y(dim-1));
        Chtilde(1+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        for ii=2:dim-1
            MM1=cos(t+k+x(ii)+1);
            MM2=cos(t+k+x(ii+1)+1)+cos(t+k+x(ii-1)+1);
            Chtilde(ii+(dim-1)*dim,1)=bb*MM1+cc*MM2;
        end
        MM1=cos(t+k+1+y(dim))+cos(t+k+x(dim)+1);
        MM2=cos(t+k+2)+cos(t+k+1+y(dim-1))+cos(t+k+x(dim-1)+1);
        Chtilde(dim2,1)=bb*MM1+cc*MM2;
        
        Chuk=mult*multMm1(dim,Chtilde); 
                   
        % Boundary DhAu(tn+k2)=Dh ut-f     
        Dhtilde(1,1)=-2*cos(t+k2+y(1))-2*cos(t+k2+x(1));
        for ii=2:dim-1
            Dhtilde(ii,1)=-2*cos(t+k2+x(ii));
        end
        Dhtilde(dim,1)=-2*cos(t+k2+1+x(1))-2*cos(t+k2+x(dim));        
        for m=2:dim-1
            Dhtilde(1+(m-1)*dim,1)=-2*cos(t+k2+y(m));
            Dhtilde(dim+(m-1)*dim,1)=-2*cos(t+k2+1+y(m));            
        end        
        Dhtilde(1+(dim-1)*dim,1)=-2*cos(t+k2+y(dim))-2*cos(t+k2+x(1)+1);
        for ii=2:dim-1
            Dhtilde(ii+(dim-1)*dim,1)=-2*cos(t+k2+x(ii)+1);            
        end
        Dhtilde(dim2,1)=-2*cos(t+k2+1+y(dim))-2*cos(t+k2+x(dim)+1);

        DhAuk2=multMm1(dim,Dhtilde);
             
        % Boundary Dh Au(tn+k))               
        Dhtilde(1,1)=-2*cos(t+k+y(1))-2*cos(t+k+x(1));
        for ii=2:dim-1
            Dhtilde(ii,1)=-2*cos(t+k+x(ii));
        end
        Dhtilde(dim,1)=-2*cos(t+k+1+x(1))-2*cos(t+k+x(dim));        
        for m=2:dim-1
            Dhtilde(1+(m-1)*dim,1)=-2*cos(t+k+y(m));
            Dhtilde(dim+(m-1)*dim,1)=-2*cos(t+k+1+y(m));            
        end
        Dhtilde(1+(dim-1)*dim,1)=-2*cos(t+k+y(dim))-2*cos(t+k+x(1)+1);
        for ii=2:dim-1
            Dhtilde(ii+(dim-1)*dim,1)=-2*cos(t+k+x(ii)+1);            
        end
        Dhtilde(dim2,1)=-2*cos(t+k+1+y(dim))-2*cos(t+k+x(dim)+1);

        DhAuk=multMm1(dim,Dhtilde);
                                      
        % Boundary DhAu(tn)
        Dhtilde(1,1)=-2*cos(t+y(1))-2*cos(t+x(1));
        for ii=2:dim-1
            Dhtilde(ii,1)=-2*cos(t+x(ii));
        end
        Dhtilde(dim,1)=-2*cos(t+1+x(1))-2*cos(t+x(dim));        
        for m=2:dim-1
            Dhtilde(1+(m-1)*dim,1)=-2*cos(t+y(m));
            Dhtilde(dim+(m-1)*dim,1)=-2*cos(t+1+y(m));            
        end
        Dhtilde(1+(dim-1)*dim,1)=-2*cos(t+y(dim))-2*cos(t+x(1)+1);
        for ii=2:dim-1
            Dhtilde(ii+(dim-1)*dim,1)=-2*cos(t+x(ii)+1);            
        end
        Dhtilde(dim2,1)=-2*cos(t+1+y(dim))-2*cos(t+x(dim)+1);
 
        DhAu=multMm1(dim,Dhtilde); 
        
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
        
    
        % Stage K2        
        G1=U.^2+vech+Chu-DhAu;              
        vecb2=zeros(dim2,2);           
        vecb2(:,1)=U;
        vecb2(:,2)=G1;      
        K2=phipmM9points(k2,JJ,vecb2,10^(-13),1,1);                                  
        
        % Stage K3         
        G2=K2.^2+vechk2+Chuk2-DhAuk2;
        vecb3=zeros(dim2,3);           
        vecb3(:,1)=U;
        vecb3(:,2)=G1;
        vecb3(:,3)=(4*G2-4*G1)/k;
        K3=phipmM9points(k2,JJ,vecb3,10^(-13),1,1);     
                                
        % Stage K4           
        G3=K3.^2+vechk2+Chuk2-DhAuk2;                
        vecb3(:,3)=(-2*G1+2*G3)/k;          
        K4=phipmM9points(k,JJ,vecb3,10^(-13),1,1);
                
        % Approximation to the exact solution at time t_{n+1}        
        G4=K4.^2+vechk+Chuk-DhAuk;           
        vecbdef=zeros(dim2,4);        
        vecbdef(:,1)=U;
        vecbdef(:,2)=G1;
        vecbdef(:,3)=(-3*G1+2*G2+2*G3-G4)/k;
        vecbdef(:,4)=(4*G1-4*G2-4*G3+4*G4)/k^2;
         
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
         
    % Error in thee infinite norm
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
