% initialize solution and bottom topography
function rhs = create_rhs_forcing(rhs, grid_obj, t, gravity)
    ngl = grid_obj.ngl;
    nel = grid_obj.nel;
    
    h0 = 0.3; ha = 0.1; om = pi; kp = pi;

    for iel=1:nel
        for igl=1:ngl
            x = grid_obj.coords(igl,iel);
            wq = grid_obj.wq(igl,iel);

            h = gravity*(h0 + ha*cos(om*t)*sin(kp*x));
            hu = gravity*(-1*ha*om/kp*sin(om*t)*cos(kp*x));
            hut = gravity*(-1*ha*om^2/kp*cos(om*t)*cos(kp*x));
            hx = gravity*(ha*kp*cos(om*t)*cos(kp*x));
            hux = gravity*(ha*om*sin(om*t)*sin(kp*x));
            
            rhs(2,igl,iel) = rhs(2,igl,iel) + wq*(hut + h*hx + 2*hu*hux/h - hx*(hu/h)^2);
                    
        end
    end
end
