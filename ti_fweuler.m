% computes the fourth order, three-stage ssp runge kutta method
% according to Gottlieb, Ketcheson, Shu - SSP Runge Kutta and multistep
% discretizations
% the cfl coefficient for this method is 2
function q2 = ti_fweuler(q0, q_ref, grid_obj, time, tol, form_method, visc, dt, iposlimit, ifilter, RHS, BC)
    q1 = q0;
    tt = time;
    
    q_bc = BC(q1, tt);
    q2 = q1 + dt * RHS(q1, q_ref, q_bc, grid_obj, tt, form_method, visc, tol);
    if (ifilter); q2 = wbfilter(q2, q_ref, grid_obj, tol); end
    if (iposlimit); q2 = poslimiter(q2, grid_obj, tol); end

end