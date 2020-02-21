function xl = find_lshore(q, grid_obj, tol)
    ngl = grid_obj.ngl;
    nel = grid_obj.nel;
    xl = grid_obj.x0;
    ibreak = 0;
    
    for ie=1:nel
        for i=1:ngl
            if (q(1,i,ie) > tol)
                ibreak = 1;
                break
            else
                xl = grid_obj.coords(i,ie);
            end
        end
        if (ibreak); break; end
    end
end