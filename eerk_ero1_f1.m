% Program that implements method given by tableau
%
% 0    |   
% 1/2  |  phi_{1,2}/2
% ------------------------------
%      |  0              phi_1 
%
% avoiding the order reduction with p=1 for problem 
% u_t=u_xx+f(t,u), x\in [0,1]
% u(t,0)=g0(t), u(t,1)=g1(t) 
% with f(t,u)=u^2+h(t,u) and with exact solution u(t,x)=cos(t+x)
% For the spatial discretization we have used the standard second-order
% difference squeme. 
% This methods appears in 
% M. Hochbruck and A. Ostermann, Explicit exponential Runge-Kutta methods 
% for semilinear parabolic problems, SIAM J. Num. Anal. 43 (2005), 1069–1090.

% Some initial values, such as h, matrix A_h0 and some more that are used
% several times along the program
% N is such that h/1/N is the grid diameter in [0,1]
N=1000;
A=-2*diag([ones(N-1,1)])+diag (ones (N-2,1),1)+diag (ones (N-2,1),-1);
h=1/N;
h2=h^2;
A=A./(h^2);
x=[h:h:1-h]';
hdiv=1/h2;

% n is such that the time step size is k=1/n.
n=5;

% The program runs for 8 different values of k, from k=1/5, in order to calculate the
% error and the order of the method
for ll=1:8

    k=1/n;
    k2=k/2;

    % U is the initial value of u(0,x)
    U=zeros(N-1,1);
    for ii=1:N-1
        U(ii)=cos(x(ii));
    end

    t=0;
    A=sparse(A);

     % For the local error r=1. For the global one, r=n
    for r=1:1
    
        % vecb and vecb2 contain the vectors that are multiplied by the
        % exponential funcions. The first column is multiplied by exp(kA),
        % the second one by phi_1(kA)...
        vecb=zeros(N-1,3);
        vecb2=zeros(N-1,2);         
    
        for ii=1:N-1     
            vecb2(ii,1)=U(ii);           
            funh1=-sin(x(ii)+t)+cos(x(ii)+t)-(cos(x(ii)+t))^2;
            vecb2(ii,2)=U(ii)^2+funh1;        
        end
            
        % Boundary conditions. They are 0 except for ii=1 and ii=N.
        % They are Dirichlet.
        vecb2(1,2)=vecb2(1,2)+hdiv*cos(t);
        vecb2(N-1,2)=vecb2(N-1,2)+hdiv*cos(1+t);
    
        % Stage
        U1=phipm(k2,A,vecb2,10^(-10),1,1);     
    
        for ii=1:N-1           
            vecb(ii,1)=U(ii);
            funh2=-sin(x(ii)+t+k2)+cos(x(ii)+t+k2)-(cos(x(ii)+t+k2))^2;
            vecb(ii,2)=U1(ii)^2+funh2;                    
        end    
    
        vecb(1,2)=vecb(1,2)+hdiv*cos(t);
        vecb(N-1,2)=vecb(N-1,2)+hdiv*cos(1+t);    

        vecb(1,3)=-hdiv*sin(t);
        vecb(N-1,3)=-hdiv*sin(1+t);
    
        % Aproximation to the exact solution at time t_{n+1}
        U=phipm(k,A,vecb,10^(-10),1,1);
      
        % New value of t
        t=t+k;        
             
    end

    % Sol contains the exact solution at time T
    sol=zeros(N-1,1);
    for ii=1:N-1
        sol(ii)=cos(x(ii)+t);
    end

    % Error in the infinite norm
    err=norm(sol-U,inf);

    % Order. When ll=1, as there is not a previous error, it can't
    % be calculated. Two consecutive errors are compared
    if ll==1
        err
        err0=err;
    else
        [err log2(err0/err)]
        err0=err;
    end
        
    % The new value of n is 2*n, and the new value of k_n=k/2. It is
    % calculated at the beginning os the next iteration
    n=2*n;

end
