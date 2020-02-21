function lshore = compute_lshore(q,grid_obj,tol)
    ngl = grid_obj.ngl;
    nel = grid_obj.nel;
    
    lshore=0;
    ibreak = 0;
    for ie=1:nel
        q1 = grid_obj.normvdm*q(1,:,ie)';
        for i=1:2*ngl+1
            x = grid_obj.normgrd(i,ie);
            if ( q1(i) > tol ); ibreak=1; break; end
            lshore = x;
        end
        if (ibreak); break; end
    end
end