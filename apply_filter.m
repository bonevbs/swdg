function q = apply_filter(q, q_ref, grid_obj, tol)
    ngl = grid_obj.ngl;
    nel = grid_obj.nel;
    
    r = zeros(ngl,1);
    u = zeros(ngl,1);
    
    for ie=1:nel
        r(:,1)  = q(1,:,ie);
        u(:,1)  = q(2,:,ie);
        
        r = grid_obj.filter(:,:,ie)*r;
        u = grid_obj.filter(:,:,ie)*u;

        q(1,:,ie) = r(:,1);
        q(2,:,ie) = u(:,1);
        
    end
end