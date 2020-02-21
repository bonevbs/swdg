% create RHS for the DG weak and strong form of the shallow water equations
function rhs = create_rhs_min(q, q_ref, grid_obj, form_method, visc, tol)
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
    rhs = create_rhs_flux(rhs, q, grid_obj, form_method, tol);
    
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
    
    % iterate over the elements
    for ie=1:nel
        % collect indices of dry nodes
        qd = zeros(3,ngl);
        wnodes = find(q(1,:,ie) > tol);
        dnodes = find(q(1,:,ie) <= tol);
        if (size(dnodes, 2) > 0 && size(wnodes, 2) > 0)
            D = grid_obj.dpsi';
            DD = D'*D;
            % assemble RHS
            sw = q(1,wnodes,ie)' + q_ref(1,wnodes,ie)';
            %uw = q(2,wnodes,ie)';
            %DD = DD(dnodes, dnodes);
            qd(1,dnodes) = (DD(dnodes, dnodes)\(-1*DD(dnodes, wnodes) * sw))';
            qd(1,wnodes) = q(1,wnodes,ie)' + q_ref(1,wnodes,ie)';
            %qd(2,dnodes) = (DD(dnodes, dnodes)\(-1*D(:,dnodes)'*D(:,wnodes) * uw))';
            %qd(2,wnodes) = q(2,wnodes,ie)';
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
            s_e     = 0.0;
            r_e     = 0.0;
            u_e     = 0.0;
            t_e     = 0.0;
            for j=1:ngl
                h_e = grid_obj.dpsi(j,i);
                if (q(1,j,ie) <= tol)
                    b_e = b_e + h_e*q_ref(1,j,ie);
                    s_e = s_e + h_e*(qd(1,j) - q_ref(1,j,ie));
                    %u_e = u_e + h_e*qd(2,j);
                else
                    b_e = b_e + h_e*q_ref(1,j,ie);
                    s_e = s_e + h_e*q(1,j,ie);
                    u_e = u_e + h_e*q(2,j,ie);
                    t_e = t_e + h_e*q(2,j,ie)*q(2,j,ie)/q(1,j,ie);
                end
            end
            
            % calc perturbation variables
            b_x     = b_e*e_x;
            s_x     = s_e*e_x;
            u_x     = u_e*e_x;
            t_x     = t_e*e_x;
            
            % calculate the flux vector
            f_r     = u_x;
            f_u     = t_x + r_i*s_x;
            
            % calculate bottom pressure gradient
            bup     = r_i*b_x;
            br      = f_r;
            bu      = f_u + bup;
            
            % add contribution to RHS
            rhs(1,i,ie) = rhs(1,i,ie) - wq*br;
            rhs(2,i,ie) = rhs(2,i,ie) - wq*bu;
        end
    end
end

% calculate the fluxes at cell interfaces both for the strong and weak form
function rhs = create_rhs_flux(rhs, q, grid_obj, form_method, tol)
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
        
        % for now we have reflecting bcs
        if (iel < 1)
            rl = q(1,1,1);
            ul = -1*q(2,1,1)/rl;
            rr = q(1,1,ier);
            ur = q(2,1,ier)/rr;
        elseif (ier < 1)
            rl = q(1,ngl,iel);
            ul = q(2,ngl,iel)/rl;
            rr = q(1,ngl,nel);
            ur = -1*q(2,ngl,nel)/rr;
        else
            rl = q(1,ngl,iel);
            ul = q(2,ngl,iel)/rl;
            rr = q(1,1,ier);
            ur = q(2,1,ier)/rr;
        end
        
        if (rl <= tol); ul = 0.0; end
        if (rr <= tol); ur = 0.0; end
        
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

% rusanov numerical flux
function [flux_r, flux_u] = rusanov_flux(rr, rl, ur, ul)
    rur = rr*ur;
    rul = rr*ul;
    
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

% analytical flux
function [flux_r, flux_u] = elemental_flux(r, u)
    ru = r*u;
    
    flux_r = ru;
    flux_u = ru*u + 0.5*r*r;
end