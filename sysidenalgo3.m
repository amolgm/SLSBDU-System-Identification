function [xk fak alphak] = sysidenalgo3(rho, A, y, Ai, yi)

p = size(A,1);
global k
k = 1;
alphak = zeros(p,k);
fak = zeros(p,k);
W = eye(p);

ya = @(ak) y + sum(diag(ak).*yi,2);
Aa = @(ak) A + diag(ak)*Ai{1} + diag(ak)*Ai{2} + diag(ak)*Ai{3} + diag(ak)*Ai{4} + ...
    diag(ak)*Ai{5} + diag(ak)*Ai{6} + diag(ak)*Ai{7} + diag(ak)*Ai{8} + diag(ak)*Ai{9} + ...
    diag(ak)*Ai{10};
Pa = @(ak) eye(p) - A*pinv(Aa(ak));
falphak = @(ak,i) ya(ak)'*Pa(ak)*(Ai{i}*pinv(Aa(ak))*ya(ak)-yi(i));

converged = false;
epsilon = 1e-6;

while ~converged
    
    for i = 1:p
        fak(i,k) = falphak(alphak(:,k),i);
    end
    for i = 1:p 
        alphak(i,k+1) = rho * fak(i,k) / norm(W*fak(:,k),2);
    end
    
    k = k + 1;
    
    if norm(alphak(:,k) - alphak(:,k-1)) < epsilon
        converged = true;
        k = k - 1;
    end
    
end

xk = pinv(Aa(alphak(:,k)))*ya(alphak(:,k));

end