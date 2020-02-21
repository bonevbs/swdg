function err = compute_l2error(q,qe,grid_obj)
    ngl = grid_obj.ngl;
    nel = grid_obj.nel;
    
    err=0.0; nrm=0.0;
    for ie=1:nel
        %q1 = q(1,:,ie);
        %q2 = q(2,:,ie);
        qe1 = qe(1,:,ie);
        qe2 = qe(2,:,ie);
        q1 = grid_obj.normvdm*q(1,:,ie)';
        q2 = grid_obj.normvdm*q(2,:,ie)';
        for i=1:2*ngl+1
            x = grid_obj.normgrd(i,ie);
            wq = grid_obj.normwq(i,ie);
            
            if ( x >=-0.5 && x<=0.5 )
            err = err + wq*( (q1(i) - qe1(i))^2 + (q2(i) - qe2(i))^2 );
            nrm = nrm + wq*( qe1(i)^2 + qe2(i)^2 );
            end
        end
    end
    err = sqrt(err/nrm);
end