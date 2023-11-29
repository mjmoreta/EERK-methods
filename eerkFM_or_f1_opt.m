% Program that implements method given by tableau
% 0    |   
% 1/2  |  phi_{1,2}/2
% ------------------------------
%      |  0              phi_1 
%
% without avoing the order reduction for problem 
% u_t=u_xx+f(t,u), x\in [0,1]
% u(t,0)=g0(t), u(1)=g1(t) 
% with f(t,u)=u^2+h(t,u) and with exact solution u(t,x)=cos(t+x)

% N is such that h/1/N is the grid diameter in [0,1]

N=400;
h=1/N;
h2=h^2;
dim=2*N-1;

% Matrices A and M are (2N-1) x (2N-1). dim=2*N-1 is their dimension.
% Matrix A_{h,0} is given by M^{-1}A. We do not compute M^{-1}, when
% necessary, we solve a system or we use phipmM

A=zeros(dim);
M=zeros(dim);

% We have N odd rows and N-1 even rows
% The first 2 rows and the last 2 rows are different. We calculate firstly odd
% rows and then even rows 

for j=1:N
    if j==1
        A(1,1)=-16;
        A(1,2)=8;        
        M(1,1)=8/15;
        M(1,2)=1/15;
    elseif j==N
        A(dim,dim-1)=8;
        A(dim,dim)=-16;         
        M(dim,dim-1)=1/15;
        M(dim,dim)=8/15;
    else
        A(2*j-1,2*j-2)=8;
        A(2*j-1,2*j-1)=-16;
        A(2*j-1,2*j)=8;       
        M(2*j-1,2*j-2)=1/15;
        M(2*j-1,2*j-1)=8/15;
        M(2*j-1,2*j)=1/15;     
    end
end

for j=1:N-1
    if j==1
        A(2,1)=8;
        A(2,2)=-14;
        A(2,3)=8;
        A(2,4)=-1;
        
        M(2,1)=1/15;
        M(2,2)=4/15;
        M(2,3)=1/15;
        M(2,4)=-1/30;
    elseif j==N-1
        A(dim-1,dim-3)=-1;
        A(dim-1,dim-2)=8;
        A(dim-1,dim-1)=-14;
        A(dim-1,dim)=8;
        
        M(dim-1,dim-3)=-1/30;
        M(dim-1,dim-2)=1/15;
        M(dim-1,dim-1)=4/15;
        M(dim-1,dim)=1/15;
    else
        A(2*j,2*j-2)=-1;
        A(2*j,2*j-1)=8;
        A(2*j,2*j)=-14;
        A(2*j,2*j+1)=8;
        A(2*j,2*j+2)=-1;       

        M(2*j,2*j-2)=-1/30;
        M(2*j,2*j-1)=1/15;
        M(2*j,2*j)=4/15;
        M(2*j,2*j+1)=1/15;
        M(2*j,2*j+2)=-1/30;        
    end
end

A=A./(3*h^2);
A=sparse(A);
M=sparse(M);

x=[h/2:h/2:1-h/2]';
hdiv=1/h2;

% The initial value of k is k=1/5. Later, k is divided by 2.
n=5;

% There are 8 iterations for 8 different values of k, starting in k=1/5.
% Later, k is divided by 2

for ll=1:8
    
    % Initial values that are needed at each iteration. It is used that 
    % u(0,x)=cos(x). 
    
    k=1/n;
    k2=k/2;    
    U=zeros(dim,1);
    for ii=1:dim       
        U(ii)=cos(x(ii));
    end
    
    t=0;
    
    % CPU time starts afer all the initial values have been calculated    
    tstart=tic;
    
    % For the local error r=1. For the global one, r=n    
    for r=1:1
        
        % vecb and vecb2 contain the vectors that are multiplied by the
        % exponential funcions. The first column is multiplied by exp(kA),
        % the second one by phi_1(kA)...
        vecb=zeros(dim,2);
        vecb2=zeros(dim,2);
            
        for ii=1:dim
            vecb2(ii,1)=U(ii);    
            funh1=-sin(x(ii)+t)+cos(x(ii)+t)-(cos(x(ii)+t))^2;
            vecb2(ii,2)=U(ii)^2+funh1;                
        end
        
        % Ch and Dh are calculated by solving a system M x=u. Here, u is
        % calculated u and then, system M x=u is solved. 
        % U=C_h bound u(t_n,x)+ D_h bound (f(t_n,u(t_n))-u_t(t_n,x)), where 
        % bound means boundary.
        % The values of u(t,0)=cos(t), u(t,1)=cos(t+1),
        % f(t,0)-u_t(t,0)=cos(t) and f(t,1)-u_t(t,0)=-cos(1+t) are used,
        % where u_t(t,0)=-sin(t), u_t(t,1)=-sin(t+1), 
        % f(t,0)=u(t,0)^2+h(t,0)=-sin(t)+cos(t) and
        % f(t,1)=u(t,1)^2+h(t,1)=-sin(t+1)+cos(t+1).

        vecind=zeros(dim,1);
        
        u0=cos(t);
        u1=cos(t+1)
        fut0=cos(t);
        fut1=cos(1+t);
        vecind(1)=8*hdiv*u0/3+fut0/15;
        vecind(2)=-hdiv*u0/3-fut0/30;
        vecind(dim-1)=-hdiv*u1/3-fut1/30;
        vecind(dim)=8*hdiv*u1/3+fut1/15;        
        vecind=sparse(vecind);     

        vecb2(:,2)=vecb2(:,2)+M\vecind;   
        
        % Stage
        U1=phipmM(k2, A, M, vecb2, 10^(-7), 1, 1);
        
        for ii=1:dim
            vecb(ii,1)=U(ii);   
            funh2=-sin(x(ii)+t+k2)+cos(x(ii)+t+k2)-(cos(x(ii)+t+k2))^2;
            vecb(ii,2)=U1(ii)^2+funh2;
        end
       
        % The boundary values are calculated at t+k/2
        u0=cos(t+k2);
        u1=cos(t+k2+1)
        fut0=cos(t+k2);
        fut1=cos(1+t+k2);
        vecind=zeros(dim,1);
        vecind(1)=8*hdiv*u0/3+fut0/15;
        vecind(2)=-hdiv*u0/3-fut0/30;
        vecind(dim-1)=-hdiv*u1/3-fut1/30;
        vecind(dim)=8*hdiv*u1/3+fut1/15;
        vecind=sparse(vecind);        
        
        vecb(:,2)=vecb(:,2)+M\vecind;

        % Approximation to the exact solution at time t_{n+1}
        U=phipmM(k, A, M, vecb, 10^(-7), 1, 1);
        
        % New value of t
        t=t+k;    
    end
    
    % CPU time comparison finishes
    telapsed=toc(tstart)

    % Sol contains the exact solution at time T
    sol=zeros(dim,1);

    for ii=1:dim
        sol(ii)=cos(x(ii)+t);
    end
    
    % Error in the discrete L2-norm and the infinite norm

    err2=sqrt(h/2)*norm(sol-U);
    err=norm(sol-U,inf);

    % The order in the discrete L2-norm is calculated as log2(err02/err2), 
    % with err02 the error obtained with k and err2 the error obtained with k/2. 
    % The order in the infinite norm is calculated in a similar way.
    % When ll=1, as there is not a previous error, it can't  be calculated. 
    % Two consecutive errors are compared
    if ll==1
        [err2, err log2(err0/err)]
        err02=err2;
        err0=err;
    else
        [err2 log2(err02/err2) err log2(erro0,err)]
        err02=err2;
        err0=err;
    end

    % The new value of n is 2*n, and the new value of k_n=k/2. It is
    % calculated at the beginning os the next iteration
    n=2*n;

end
