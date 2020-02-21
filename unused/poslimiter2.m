function q = poslimiter2(q, grid_obj, tol)
    ngl = grid_obj.ngl;
    nel = grid_obj.nel;
    
    for ie=1:nel

        for i=1:ngl
            if (q(1,i,ie) < tol); q(2,i,ie) = 0.0; q(1,i,ie) = 0.0; end
        end

    end
end