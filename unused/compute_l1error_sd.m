function err = compute_l1error_sd(q,qe,grid_obj)
    ngl = grid_obj.ngl;
    nel = grid_obj.nel;
    
    err=0.0; nrm=0.0;
    for ie=1:nel
        q1 = grid_obj.normvdm*q(1,:,ie)';
        q2 = grid_obj.normvdm*q(2,:,ie)';
        for i=1:2*ngl+1
            x = grid_obj.normgrd(i,ie);
            wq = grid_obj.normwq(i,ie);
            
            if (x >= 0.0 && x <=1.0) 
                err = err + wq*sqrt( (q1(i) - qe(1,i,ie))^2 + (q2(i) - qe(2,i,ie))^2 );
                nrm = nrm + wq*sqrt( qe(1,i,ie)^2 + qe(2,i,ie)^2 );
            end
        end
    end
    err = err/nrm;
end