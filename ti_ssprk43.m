% computes the third order, three-stage ssp runge kutta method
% according to Gottlieb, Ketcheson, Shu - SSP Runge Kutta and multistep
% discretizations
% the cfl coefficient for this method is 2
% should be ssprk33 really
function q2 = ti_ssprk43(q0, q_ref, grid_obj, time, tol, form_method, visc, dt, iposlimit, ifilter, RHS, BC)
    q1 = q0;
    tt = time;
    
    q_bc = BC(q1, tt);
    q2 = q1 + dt/2 * RHS(q1, q_ref, q_bc, grid_obj, tt, form_method, visc, tol);
    if (ifilter); q2 = wbfilter(q2, q_ref, grid_obj, tol); end
    if (iposlimit); q2 = poslimiter(q2, grid_obj, tol); end
    tt = tt+dt/2;
    
    q_bc = BC(q2, tt);
    q2 = q2 + dt/2 * RHS(q2, q_ref, q_bc, grid_obj, tt, form_method, visc, tol);
    if (ifilter); q2 = wbfilter(q2, q_ref, grid_obj, tol); end
    if (iposlimit); q2 = poslimiter(q2, grid_obj, tol); end
    tt = tt+dt/2;
    
    q_bc = BC(q2, tt);
    q2 = 2/3 * q1 + 1/3 * (q2 + dt/2 * RHS(q2, q_ref, q_bc, grid_obj, tt, form_method, visc, tol));
    if (ifilter); q2 = wbfilter(q2, q_ref, grid_obj, tol); end
    if (iposlimit); q2 = poslimiter(q2, grid_obj, tol); end
    tt = (tt+dt/2)/3 + 2*time/3;
    
    q_bc = BC(q2, tt);
    q2 = q2 + dt/2 * RHS(q2, q_ref, q_bc, grid_obj, tt, form_method, visc, tol);
    if (ifilter); q2 = wbfilter(q2, q_ref, grid_obj, tol); end
    if (iposlimit); q2 = poslimiter(q2, grid_obj, tol); end
end