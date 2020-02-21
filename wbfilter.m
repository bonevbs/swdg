function q = wbfilter(q, q_ref, grid_obj, tol)
    ngl = grid_obj.ngl;
    nel = grid_obj.nel;
    
    iselective = 0;
    
    r = zeros(ngl,1);
    u = zeros(ngl,1);
    
    for ie=1:nel
        if (iselective && (min(q(1,:,ie)) > tol) )
            continue
        elseif (min(q(1,:,ie)) > tol)
            r(:,1) = q(1,:,ie) + q_ref(1,:,ie);
            idry = 0;
        else
            dind = find(q(1,:,ie) < tol);
            r(:,1) = q(1,:,ie) + q_ref(1,:,ie);
            r(dind, 1) = 0.1005*9.81;
            idry = 1;
        end
        u(:,1)  = q(2,:,ie);
        b(:,1)  = q_ref(1,:,ie);
        
        r = grid_obj.filter(:,:,ie)*r;
        u = grid_obj.filter(:,:,ie)*u;
        
        if (idry)
            q(1,:,ie) = r(:,1) - b(:,1);
            q(2,:,ie) = u(:,1);
            q(1, q(1,:,ie)<0, ie) = 0.0;
        else
            q(1,:,ie) = r(:,1) - b(:,1);
            q(2,:,ie) = u(:,1);
        end

    end
end