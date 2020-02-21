% create RHS for the DG weak and strong form of the shallow water equations
function rhs = create_rhs_fd(q, q_ref, q_bc, grid_obj, time, form_method, visc, tol)
    nel = grid_obj.nel;
    ngl = grid_obj.ngl;
    nvar = size(q,1);
    rhs = zeros(size(q));
    
    % calculate the volume integrals
    if strcmp(form_method, 'weak')
        rhs = create_rhs_volume_weak(rhs, q, q_ref, grid_obj, visc, tol);
    elseif strcmp(form_method, 'strong')
        rhs = create_rhs_volume_strong(rhs, q, q_ref, grid_obj, visc, tol);
    else
        error('Unknown form method.')
    end
    
    % calculate the surface integrals
    rhs = create_rhs_flux(rhs, q, q_bc, grid_obj, form_method, tol);
    
    % calculate the forcing term
    % this is a hack for creating results
    % rhs = create_rhs_forcing(rhs, grid_obj, time, 9.81);
    
    % rescale with the mass matrix
    for i = 1:nvar
        rhs(i,:,:) = reshape(rhs(i,:,:),[ngl,nel])./grid_obj.wq;
    end
end

% calculate the volume integral of the weak form
function rhs = create_rhs_volume_weak(rhs, q, q_ref, grid_obj, visc, tol)
    ngl = grid_obj.ngl;
    nel = grid_obj.nel;
    
    % iterate over the elements
    for ie=1:nel
        for i=1:ngl
            % copy local variables
            wq = grid_obj.wq(i,ie);
            e_x = grid_obj.ksi_x(i, ie);
            r_i = q(1,i,ie);
            ru_i = q(2,i,ie);
            u_i = ru_i/r_i;
            if (r_i <= tol); r_i = 0; u_i = 0; ru_i = 0; end
            
            % calculate the flux vector
            f_r = ru_i;
            f_u = ru_i*u_i + 0.5*r_i*r_i;
            f_re = e_x*f_r;
            f_ue = e_x*f_u;
            
            % add the contribution of the i-th flux F(q_i) through dpsi(k,i)
            % also add up all influences of the bottom topography from all
            % nodes within the cell
            b_e = 0;
            for j=1:ngl
                % the flux at node i influences contributes to all other
                % nodes k
                dpsi = grid_obj.dpsi(j,i);
                rhs(1,j,ie) = rhs(1,j,ie) + f_re*dpsi*wq;
                rhs(2,j,ie) = rhs(2,j,ie) + f_ue*dpsi*wq;
                
                % it is exactly the other way around for the source terms
                % all other nodes k influence the source term through the
                % derivatives of the k-th Lagrange polynomials ath the i-th
                % node
                b_e = b_e + q_ref(1,j,ie)*grid_obj.dpsi(j,i);
            end
            bup = r_i*b_e*e_x;
            %bup = r_k*0.5;
            rhs(2,i,ie) = rhs(2,i,ie) + wq*(-bup);
        end
    end
end

% calculate the volume integral of the strong form
function rhs = create_rhs_volume_strong(rhs, q, q_ref, grid_obj, visc, tol)
    ngl = grid_obj.ngl;
    nel = grid_obj.nel;
    xgl = grid_obj.xgl;
    
    % iterate over the elements
    for ie=1:nel
        % determine whether element is dry
        idry = 0;
        for i=1:ngl
            if (q(1,i,ie) <= tol); idry=1; end
        end
        
        for i=1:ngl
            % copy local variables
            wq      = grid_obj.wq(i,ie);
            e_x     = grid_obj.ksi_x(i, ie);
            r_i     = q(1,i,ie);
            u_i     = q(2,i,ie);
            b_i     = q_ref(1,i,ie);
            
            % calculate derivatives at x_i
            b_e     = 0.0;
            r_e     = 0.0;
            u_e     = 0.0;
            t_e     = 0.0;
            % switch to finite differences in partly dry cells
            if (idry == 1)
                il=max(i-1,1); ir=min(i+1,ngl);
                if (q(1,il,ie) <= tol); il=i; end
                if (q(1,ir,ie) <= tol); ir=i; end
                if (il ~= ir)
                    b_e = (q_ref(1,ir,ie)-q_ref(1,il,ie))/(xgl(ir)-xgl(il));
                    r_e = (q(1,ir,ie)-q(1,il,ie))/(xgl(ir)-xgl(il));
                end
                if (q(1,i,ie) <= tol)
                    b_e = 0.0;
                    r_e = 0.0;
                end
                for j=1:ngl
                    if (q(1,j,ie) > tol) 
                        h_e = grid_obj.dpsi(j,i);
                        u_e = u_e + h_e*q(2,j,ie);
                        t_e = t_e + h_e*q(2,j,ie)*q(2,j,ie)/q(1,j,ie);
                    end
                end
            else
                for j=1:ngl
                    h_e = grid_obj.dpsi(j,i);
                    b_e = b_e + h_e*q_ref(1,j,ie);
                    r_e = r_e + h_e*q(1,j,ie);
                    u_e = u_e + h_e*q(2,j,ie);
                    t_e = t_e + h_e*q(2,j,ie)*q(2,j,ie)/q(1,j,ie);
                end
            end
            
            % calc perturbation variables
            b_x     = b_e*e_x;
            r_x     = r_e*e_x;
            u_x     = u_e*e_x;
            t_x     = t_e*e_x;
            
            % calculate the flux vector
            f_r     = u_x;
            f_u     = t_x + r_i*r_x;
            
            % calculate bottom pressure gradient
            bup     = r_i*b_x;
            
            % calculate Manning friction
            %vu      = visc*safe_div(u_i*abs(u_i)*9.81^(4/3),r_i^(7/3)*10E+5,0.0)*10E+5;
            vu      = 0.0;
            
            % assemble RHS
            br      = f_r;
            bu      = f_u + bup + vu;
            
            % add contribution to RHS
            rhs(1,i,ie) = rhs(1,i,ie) - wq*br;
            rhs(2,i,ie) = rhs(2,i,ie) - wq*bu;
            
            % add viscosity to all nodes
            if (visc > 0)
                uvisc_e = 0.0;
                for l=1:ngl
                    h_e = grid_obj.dpsi(l,i);
                    uvisc_e = uvisc_e+h_e*safe_div(q(2,j,ie),q(1,j,ie)*10E+5,0.0)*10E+5;
                end
                uvisc_x = uvisc_e*e_x;
                
                % distribute onto all nodes in the cell
                for l=1:ngl
                    h_e = grid_obj.dpsi(l,i);
                    dhudx = h_e*e_x*uvisc_x;
                    rhs(2,l,ie) = rhs(2,l,ie) - wq*visc*r_i*dhudx;
                end
            end
        end
    end
end

% calculate the fluxes at cell interfaces both for the strong and weak form
function rhs = create_rhs_flux(rhs, q, q_bc, grid_obj, form_method, tol)
    ngl = grid_obj.ngl;
    nel = grid_obj.nel;
    nfa = grid_obj.nfa;
    
    if strcmp(form_method, 'weak')
        iflux=0;
    elseif strcmp(form_method, 'strong')
        iflux=1;
    else
        error('Unknown form_method.')
    end
    
    % iterate over the faces
    for ifa=1:nfa
        % pointers to left and right elements
        iel = grid_obj.faces(1,ifa);
        ier = grid_obj.faces(2,ifa);
        
        % for now we only use analytic bcs (for correct error norms)
        if (iel < 1)
            rl = q_bc(1,1);
            ul = safe_div(q_bc(2,1), rl, 0.0);
            rr = q(1,1,ier);
            ur = safe_div(q(2,1,ier), rr, 0.0);
        elseif (ier < 1)
            rl = q(1,ngl,iel);
            ul = safe_div(q(2,ngl,iel), rl, 0.0);
            rr = q_bc(1,2);
            ur = safe_div(q_bc(2,2), rr, 0.0);
        else
            rl = q(1,ngl,iel);
            ul = safe_div(q(2,ngl,iel), rl, 0.0);
            rr = q(1,1,ier);
            ur = safe_div(q(2,1,ier), rr, 0.0);
        end
        
        if (rl <= tol); ul = 0.0; end
        if (rr <= tol); ur = 0.0; end
        
        %[flux_r, flux_u]    = roe_flux(rr, rl, ur, ul);
        [flux_r, flux_u]    = rusanov_flux(rr, rl, ur, ul);
        [flux_rr, flux_ur]  = elemental_flux(rr, ur);
        [flux_rl, flux_ul]  = elemental_flux(rl, ul);
        
        % update left and ride side of the face
        if (iel ~= -1)
            rhs(1,ngl,iel) = rhs(1,ngl,iel) + 1*(iflux*flux_rl - flux_r);
            rhs(2,ngl,iel) = rhs(2,ngl,iel) + 1*(iflux*flux_ul - flux_u);
        end
        if (ier ~= -1)
            rhs(1,1,ier) = rhs(1,1,ier) + 1*(-iflux*flux_rr + flux_r);
            rhs(2,1,ier) = rhs(2,1,ier) + 1*(-iflux*flux_ur + flux_u);
        end
        
    end
end

% Rusanov numerical flux
function [flux_r, flux_u] = rusanov_flux(rr, rl, ur, ul)
    rur = rr*ur;
    rul = rl*ul;
    
    f_rr = rur;
    f_ur = rur*ur + 0.5*rr*rr;
    f_rl = rul;
    f_ul = rul*ul + 0.5*rl*rl;
    
    % maximum signal velocity
    claml = abs(ul) + sqrt(max([rl, 0.0]));
    clamr = abs(ur) + sqrt(max([rr, 0.0]));
    clam = max(claml, clamr);
    
    % average terms
    flux_r = f_rr + f_rl;
    flux_u = f_ur + f_ul;
    
    % dissipative terms
    diss_r = clam*(rr - rl);
    diss_u = clam*(rur - rul);
    
    % assemble the Rusanov flux
    flux_r = 0.5*(flux_r - diss_r);
    flux_u = 0.5*(flux_u - diss_u);
end

% Roe numerical flux
function [flux_r, flux_u] = roe_flux(rr, rl, ur, ul)
    rur = rr*ur;
    rul = rl*ul;
    
    f_rr = rur;
    f_ur = rur*ur + 0.5*rr*rr;
    f_rl = rul;
    f_ul = rul*ul + 0.5*rl*rl;
    
    % construct intermediate state and linearized jacobian
    c_bar = sqrt( max([0.5*(rr + rl),0.0]) );
    u_bar = safe_div(sqrt(max([rl, 0.0]))*ul + sqrt(max([rr, 0.0]))*ur, sqrt(max([rl, 0.0])) + sqrt(max([rr, 0.0])), 0.0 );
    l1 = u_bar - c_bar;
    l2 = u_bar + c_bar;
    a1 = safe_div( -(rur - rul) + (u_bar + c_bar)*(rr - rl), 2*c_bar, 0.0);
    a2 = safe_div(  (rur - rul) - (u_bar - c_bar)*(rr - rl), 2*c_bar, 0.0);
    
    % average terms
    flux_r = f_rr + f_rl;
    flux_u = f_ur + f_ul;

    % dissipative terms
    diss_r = (abs(l1)*a1 + abs(l2)*a2);
    diss_u = (abs(l1)*a1*l1 + abs(l2)*a2*l2);
    
    % assemble the Rusanov flux
    flux_r = 0.5*(flux_r - diss_r);
    flux_u = 0.5*(flux_u - diss_u);
end

% analytical flux
function [flux_r, flux_u] = elemental_flux(r, u)
    ru = r*u;
    
    flux_r = ru;
    flux_u = ru*u + 0.5*r*r;
end