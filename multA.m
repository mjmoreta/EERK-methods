% Program that calculated A*U, donde with A the matrix that appears in the 
% 9 point formulas. Matrix A is very large, so it is not calculated

function Au=multA(dim,u)

Au=zeros(dim*dim,1);

A11=-(10/3)*diag([ones(dim,1)])+(2/3)*(diag(ones(dim-1,1),1)+diag(ones(dim-1,1),-1));
A12=(2/3)*diag([ones(dim,1)])+(1/6)*(diag(ones(dim-1,1),1)+diag(ones(dim-1,1),-1));

A11=sparse(A11);
A12=sparse(A12);


% Rows 1 to dim
Au(1:dim)=A11*u(1:dim)+A12*u(dim+1:2*dim);
for ll=2:dim-1
    % Values of the diagonals. THey are used for the rest of the terms.
    % Rows (ll-1)(J-1)+1 to ll(J-1)       
    Dm1=(ll-2)*dim;
    D=(ll-1)*dim;
    Dp1=ll*dim;        
    Au(D+1:D+dim)=A11*u(D+1:D+dim)+A12*(u(Dm1+1:Dm1+dim)+u(Dp1+1:Dp1+dim));   
end
% Last J-1 rows, from (J-1)*(J-2)+1 to (J-1)*(J-1)
Dm1=(dim-2)*dim;
D=(dim-1)*dim;
Au(D+1:D+dim)=A11*u(D+1:D+dim)+A12*u(Dm1+1:Dm1+dim);