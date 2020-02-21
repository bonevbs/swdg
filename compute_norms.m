% computes the error using the ki
function [E1, E2] = compute_norms(q,qe,grid_obj)
    ngl = grid_obj.ngl;
    nel = grid_obj.nel;
    
    E1=0.0; E2=0.0; N1=0.0; N2=0.0;
    for ie=1:nel
        for i=1:ngl
            wq = grid_obj.wq(i,ie);
            
            E1 = E1 + wq*(q(1,i,ie) - qe(1,i,ie))^2;
            E2 = E2 + wq*(q(2,i,ie) - qe(2,i,ie))^2;
            N1 = N1 + wq*qe(1,i,ie)^2;
            N2 = N2 + wq*qe(2,i,ie)^2;
        end
    end
    E1 = sqrt(E1/N1);
    E2 = sqrt(E2/N1);
end