% initialize solution and bottom topography
function [q, q_ref] = initial_conditions(icase, coords, t, gravity)
    [ngl, nel] = size(coords);
    
    % initialize solution arrays
    q       = zeros(2, ngl, nel);
    q_ref   = zeros(2, ngl, nel);

    for iel=1:nel
        for igl=1:ngl
            
            % select initial condition
            switch icase
                
                % still water and no bottom topography
                case 1
                    x = coords(igl,iel);
                    
                    q(1,igl,iel)      = gravity*1.0;
                    q(2,igl,iel)      = gravity*0.0;
                    q_ref(1,igl,iel)  = 0.0;
                    q_ref(2,igl,iel)  = 0.0;
                    
                % SOMETHING
                case 2
                    
                    
                % simple dam break on a wet domain at 0
                case 3
                    % if hl or hr are changed, cm has to be recalculated!
                    hl = 2.0; hr = 1.0; x0 = 0.0; hm = 1.45384; cm = sqrt(gravity*hm);
                    xA = x0-t*sqrt(gravity*hl);
                    xB = x0+t*(2*sqrt(gravity*hl)-3*cm);
                    xC = x0+t*2*cm^2*(sqrt(gravity*hl)-cm)/(cm^2-gravity*hr);
                    x = coords(igl,iel);
                    
                    if (x <= xA)
                        q(1,igl,iel)  = gravity*hl;
                        q(2,igl,iel)  = gravity*0.0;
                    elseif (x <= xB)
                        q(1,igl,iel)  = 4/9*( sqrt(gravity*hl) - (x-x0)/2/t )^2;
                        q(2,igl,iel)  = q(1,igl,iel)*2/3*( (x-x0)/t + sqrt(gravity*hl));
                    elseif (x <= xC)
                        q(1,igl,iel)  = cm^2;
                        q(2,igl,iel)  = q(1,igl,iel)*2*( sqrt(gravity*hl) - cm );
                    else
                        q(1,igl,iel)  = gravity*hr;
                        q(2,igl,iel)  = gravity*0.0;
                    end
                
                % dam break on dry domain at 0
                case 4
                    % if hl or hr are changed, cm has to be recalculated!
                    hl = 0.1; x0 = 0.0;
                    xA = x0-t*sqrt(gravity*hl);
                    xB = x0+2*t*sqrt(gravity*hl);
                    x = coords(igl,iel);
                    
                    if (x <= xA)
                        q(1,igl,iel)  = gravity*hl;
                        q(2,igl,iel)  = gravity*0.0;
                    elseif (x <= xB)
                        q(1,igl,iel)  = 4/9*( sqrt(gravity*hl) - (x-x0)/2/t )^2;
                        q(2,igl,iel)  = q(1,igl,iel)*2/3*( (x-x0)/t + sqrt(gravity*hl));
                    else
                        q(1,igl,iel)  = gravity*0.0;
                        q(2,igl,iel)  = gravity*0.0;
                    end
                    
                % sloping beach
                case 5
                    x = coords(igl,iel);
                    x0 = -0.0;
                    h0 = 0.1005;
                    
                    q_ref(1,igl,iel)  = gravity*0.5*(x-x0);
                    q_ref(2,igl,iel)  = 0.0;
                    
                    q(1,igl,iel)  = max([gravity*h0 - q_ref(1,igl,iel), 0.0]);
                    %q(1,igl,iel)  = gravity*h0 - q_ref(1,igl,iel);
                    q(2,igl,iel)  = gravity*0.0;
                    
                % oscillating lake
                case 6
                    x = coords(igl,iel);
                    h0 = 0.1005;
                    l = 1.0;
                    amp = 0.1;
                    freq = sqrt(2*gravity*h0)/l;
                    b = h0*(x/l)^2;
                    hb = h0 + 2*amp*h0*l^-2*cos(freq*t)*(x-amp/2*cos(freq*t));
                    
                    u = -amp*freq*sin(freq*t);
                    
                    q(1,igl,iel)        = max([gravity*(hb - b), 0.0]);
                    q(2,igl,iel)        = q(1,igl,iel)*u;
                    
                    q_ref(1,igl,iel)    = gravity*b;
                    q_ref(2,igl,iel)    = 0.0;
                    
                % test case for testing the lagrange basis functions
                case 100
                    q(1,igl,iel)      = grid_obj.psi(1,igl);
                    
                case 200
                    h0 = 0.3; ha = 0.1; om = pi; kp = pi;
                    x = coords(igl,iel);
                    
                    q(1,igl,iel)        = gravity*(h0 + ha*cos(om*t)*sin(kp*x));
                    q(2,igl,iel)        = gravity*(-1*ha*om/kp*sin(om*t)*cos(kp*x));
                    
                    q_ref(1,igl,iel)    = 0.0;
                    q_ref(2,igl,iel)    = 0.0;
                    
                otherwise
                    error('Unknown icase.')
            end
            
        end
    end
end
