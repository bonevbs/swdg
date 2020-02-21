% calculates the vandermonde matrix of the lagrange polynomials on xfrom
% evaluated at xto
function vdm = vandermonde(xto,xin)
    nto = size(xto,2);
    nin = size(xin,2);
    
    vdm = zeros(nto,nin);
    for j=1:nin
        for i=1:nto
            y = 1;
            for k=1:nin
                if (k~=j)
                    y = y*(xto(i) - xin(k))/(xin(j) - xin(k));
                end
            end
            vdm(i,j) = y;
        end
    end
end