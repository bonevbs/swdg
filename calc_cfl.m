% just a function to get the local cfl
function cfl = calc_cfl(q, dt, grid_obj)
    ngl = grid_obj.ngl;
    nel = grid_obj.nel;
    
    cfl = zeros(ngl,nel);
    
    for ie=1:nel
        for k=1:ngl
            r = q(1,k,ie);
            ru = q(2,k,ie);
            
            u = safe_div(ru, r, 0.0);
            
            %cfl(k,ie) = (abs(u) + sqrt(r))*dt*grid_obj.wgl(k)/grid_obj.wq(k,ie);
            cfl(k,ie) = (abs(u) + sqrt(max([r,0])))*dt/grid_obj.wq(k,ie);
        end
    end
end