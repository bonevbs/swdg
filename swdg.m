function [err1, err2, errm, erre] = swdg(icase,time_initial,time_final,dt,x0,x1,nop,nel,form_method,bc_type,ifilter,iposlimit,tol,iverbose)
    % constants
    % wetting/drying threshold
    %tol         = 1E-6;

    % scenario parameters
    gravity     = 9.81;
    visc        = 0.0;
    % plotting parameters
    iplot       = 0;
    plot_method = 'interp';
    istep       = 0;

    %%
    % initialize grid
    grid_obj = grid1d(nel, nop, x0, x1);

    % create and plot initial conditions
    time = time_initial;
    [q0, q_ref] = initial_conditions(icase, grid_obj.coords, time, 9.81);
    if (iplot); plot_solution(q0, q0, q_ref, grid_obj, gravity, plot_method, tol, 0, 0); end
    if (istep); pause(); end
    
    % create function handle for boundary conditions
    if strcmp(bc_type, 'exact')
        calc_bc = @(q, t) initial_conditions(icase, [x0,x1], t, gravity);
    elseif strcmp(bc_type, 'reflecting')
        calc_bc     = @(q, t) [q(:,1,1), q(:,end,end)];
    end

    %%
    % set errors for the case of termination
    err1 = NaN; err2 = NaN;
    
    %%
    % time loop
    q = q0;
    while (time < time_final)
        % calculate the cfl number
        cfl = calc_cfl(q0, dt, grid_obj);
        maxcfl = max(max(cfl));
        %while (maxcfl > 0.5*grid_obj.wgl(1)); dt=dt/2; maxcfl=maxcfl/2; fprintf('Reducing timestep: %f -> %f \n', 2*dt, dt); end
        if (iverbose); fprintf('time: %f maxcfl: %f \n', time, maxcfl); end
        % perform the actual timestepping
        %q = ti_ssprk43(q, q_ref, grid_obj, time, tol, form_method, visc, dt, iposlimit, ifilter, @create_rhs, calc_bc);
        q = ti_ssprk43(q, q_ref, grid_obj, time, tol, form_method, visc, dt, iposlimit, ifilter, @create_rhs_fd, calc_bc);
        %q = ti_fweuler(q, q_ref, grid_obj, time, tol, form_method, visc, dt, iposlimit, ifilter, @create_rhs_fd, calc_bc);
        if any(q(1,:,:) < 0.0); warning('negative value encountered.'); end
        % warn if NaN is created
        if any(isnan(reshape(q,[],1))); return; end
        time = time + dt;
        % plot during simulation
        if (iplot); q_ex = initial_conditions(icase, grid_obj.coords, time, gravity); plot_solution(q, q_ex, q_ref, grid_obj, gravity, plot_method, 0, ''); end
        % calculate exact solution and compute errors
        if (iverbose)
            q_ex = initial_conditions(icase, grid_obj.normgrd, time, gravity);
            err1 = compute_l1error(q, q_ex, grid_obj);
            err2 = compute_l2error(q, q_ex, grid_obj);
            fprintf('err1: %f \t err2: %f \n', err1, err2)
        end
        % stepwise
        if (istep); pause(); end
    end

    %%
    % compute errors
    q_ex = initial_conditions(icase, grid_obj.normgrd, time, gravity);
    err1 = compute_l1error(q, q_ex, grid_obj);
    err2 = compute_l2error(q, q_ex, grid_obj);
end