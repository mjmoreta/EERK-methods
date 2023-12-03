% Program that computes M^-1}*U, with M the mass matrix that appears in the
% 9 point formula. Instead of calculating M^{-1}U, we solve system Mx=U 
% by using Gauss-Seidel

function Mu=multMm1(dim,u)

Tol=10^(-15);

Mu=zeros(dim*dim,1);

Mul=1/8;
Ep=10;
n=1;

while Ep >= Tol && n<=10000
    Ep=0;
    Soln=Mul*(u(1)-Mu(2)-Mu(dim+1));
    Ep=Ep+(Soln-Mu(1))^2;
    Mu(1)=Soln;
    
    for ll=2:dim+1-2
        Soln=Mul*(u(ll)-Mu(ll-1)-Mu(ll+1)-Mu(dim+ll));
	    Ep=Ep+(Soln-Mu(ll))^2;
	    Mu(ll)=Soln;
    end
    
    Soln=Mul*(u(dim)-Mu(dim-1)-Mu(2*dim));
	Ep=Ep+(Soln-Mu(dim))^2;
    Mu(dim)=Soln;
    
    for ll=2:dim-1
        Dm1=(ll-2)*dim;
        D=(ll-1)*dim;
        Dp1=ll*dim;
        
        Soln=Mul*(u(D+1)-Mu(Dm1+1)-Mu(D+2)-Mu(Dp1+1));
	    Ep=Ep+(Soln-Mu(D+1))^2;
	    Mu(D+1)=Soln;
	    for jj=2:dim-1
	      Soln=Mul*(u(D+jj)-Mu(Dm1+jj)-Mu(D+jj-1)-Mu(D+jj+1)-Mu(Dp1+jj));
	      Ep=Ep+(Soln-Mu(D+jj))^2;
	      Mu(D+jj)=Soln;
        end
        Soln=Mul*(u(D+dim)-Mu(Dm1+dim)-Mu(D+dim-1)-Mu(Dp1+dim));
	    Ep=Ep+(Soln-Mu(D+dim))^2;
	    Mu(D+dim)=Soln;
    end
           
    Dm1=(dim-2)*dim;
    D=(dim-1)*dim;
    Soln=Mul*(u(D+1)-Mu(Dm1+1)-Mu(D+2));
    Ep=Ep+(Soln-Mu(D+1))^2;
    Mu(D+1)=Soln;
    for ll=2:dim-1
        Soln=Mul*(u(D+ll)-Mu(Dm1+ll)-Mu(D+ll-1)-Mu(D+ll+1));
	    Ep=Ep+(Soln-Mu(D+ll))^2;
	    Mu(D+ll)=Soln;
    end
	Soln=Mul*(u(D+dim)-Mu(Dm1+dim)-Mu(D+dim-1));
	Ep=Ep+(Soln-Mu(D+dim))^2;
    Mu(D+dim)=Soln;
    Ep=sqrt(Ep);
    n=n+1;
end



