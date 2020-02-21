function q = poslimiter(q, grid_obj, tol)
    ngl = grid_obj.ngl;
    nel = grid_obj.nel;
    
    for ie=1:nel
        
        % calculate the average and find the minimum
        ra = 0.0;
        ua = 0.0;
        vol = 0.0;
        rm = q(1,1,ie);
        for i=1:ngl
            % we are directly dividing by the volume of the cell
            wq = grid_obj.wq(i,ie);
            ra = ra + wq*q(1,i,ie);
            ua = ua + wq*q(2,i,ie);
            vol = vol + wq;
            if (q(1,i,ie) < rm); rm = q(1,i,ie); end
        end
        ra = ra/vol;
        ua = ua/vol;
        
        if (ra < 0.0)
            q(1,:,ie) = 0.0;
            q(2,:,ie) = 0.0;
        elseif (rm < 0.0)
            % determine scaling parameter
            theta = min(1, ra/(ra-rm));
        
            % apply the scaling to the solution
            q(1,:,ie) = theta*(q(1,:,ie) - ra) + ra;
            q(2,:,ie) = theta*(q(2,:,ie) - ua) + ua;
        end
        
        % make sure u = 0 if water height is too low
        for i=1:ngl
            if (q(1,i,ie) < tol); q(2,i,ie) = 0.0; end
        end
    end
end