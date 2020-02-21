% calculates the lagrange basis functions psi and the first derivative
% dpsi on a set of nodes xgl
% psi(k,l) contains the k-th lagrange polynomial evaluated at x_l
function [psi, dpsi] = lagrange_basis(xgl)
    ngl = size(xgl, 2);
    
    psi = eye(ngl,ngl);
    dpsi = ones(ngl,ngl);
    
    % calculate the entries in dpsi(j,i)
    for i=1:ngl
        for j=1:ngl
            y = 0;
            for l=1:ngl
                if not(l==j)
                    k = 1/(xgl(j)-xgl(l));
                    for m=1:ngl
                        if not(m==j) && not(m==l)
                            k = k*(xgl(i)-xgl(m))/(xgl(j)-xgl(m));
                        end
                    end
                    y = y + k;
                end
            end
            dpsi(j,i) = y;
        end
    end
end