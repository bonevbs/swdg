% Timestepping parameters
time_initial= 0.0;
time_final  = 1.0;
dt          = 0.00005;

% DG parameters
form_method = 'strong';
% h/p parameters
nop         = 3;
nel         = 100;
% wetting/drying threshold
tol         = 1E-4;
% post-processing parameters
ifilter     = 0;
iposlimit   = 1;

% scenario parameters
icase       = 6;
x0          = -1.5;
x1          = 1.5;
gravity     = 9.81;
visc        = 0;

% plotting parameters
iplot       = 1;
plot_method = 'interp';
plot_times  = time_initial:0.01:time_final;
istep       = 0;
iverbose    = 1;

% write out plots
iwrite      = 1;
write_times = time_initial:0.5:time_final;

% boundary conditions
bc_type     = 'exact';

%%
% initialize grid
grid_obj = grid1d(nel, nop, x0, x1);

% initialize VideoWriter object
fname = sprintf('case%d_nop%d_nel%d.avi',icase,nop,nel);
if (iplot); writerObj = VideoWriter(fname); writerObj.FrameRate = 25; open(writerObj); end

% create and plot initial conditions
time = time_initial;
[q0, q_ref] = initial_conditions(icase, grid_obj.coords, time, 9.81);
if (iplot); plot_solution(q0, q_ref, grid_obj, icase, time, gravity, plot_method, 0, 0); frame = getframe(gcf); writeVideo(writerObj, frame); end
if (istep); pause(); end

% create function handle for boundary conditions
if strcmp(bc_type, 'exact')
    calc_bc = @(q, t) initial_conditions(icase, [x0,x1], t, gravity);
elseif strcmp(bc_type, 'reflecting')
    calc_bc = @(q, t) [q(:,1,1), q(:,end,end)];
end

ae1 = []; ae2 = []; timings=[];

%%
% time loop
q = q0;
while (time <= time_final+dt)
    % calculate the cfl number
    cfl = calc_cfl(q0, dt, grid_obj);
    maxcfl = max(max(cfl));
    %while (maxcfl > 0.5*grid_obj.wgl(1)); dt=dt/2; maxcfl=maxcfl/2; fprintf('Reducing timestep: %f -> %f \n', 2*dt, dt); end
    fprintf('time: %f maxcfl: %f \n', time, maxcfl)
    % calculate norms
    if (iverbose)
        q_ex = initial_conditions(icase, grid_obj.normgrd, time, gravity);
        err1 = compute_l1error(q, q_ex, grid_obj);
        err2 = compute_l2error(q, q_ex, grid_obj);
        fprintf('err1: %e \t err2: %e \n', err1, err2)
    end
    % calculate energy errors
     if (ismembertol(time, plot_times, 0.5*dt))
        timings = [timings, time];
        q_ex = initial_conditions(icase, grid_obj.coords , 0, gravity);
        [merr, eerr] = compute_energy_errors(q,q_ex,q_ref,grid_obj);
        ae1 = [ae1,merr]; ae2 = [ae2, eerr];
        fprintf('errm: %e \t erre: %e \n', merr, eerr)
    end
    % write out to file
    if (iwrite && ismembertol(time, write_times, 0.5*dt))
        q_ex = initial_conditions(icase, grid_obj.coords, time, 9.81);
        %write output file for pgfplots
        fname = sprintf('case%d_t%.2f_nop%d_nel%d.dat',icase,time,nop,nel);
        plot_solution(q, q_ref, grid_obj, icase, time, gravity, plot_method, 1, fname)
    end
    % perform the actual timestepping
    %q = ti_fweuler(q, q_ref, grid_obj, time, tol, form_method, visc, dt, iposlimit, ifilter, @create_rhs, calc_bc);
    q = ti_ssprk43(q, q_ref, grid_obj, time, tol, form_method, visc, dt, iposlimit, ifilter, @create_rhs_fd, calc_bc);
    % warn if NaN is created
    if any(isnan(reshape(q,[],1))); warning('NaN encountered.'); end
    time = time + dt;
    % plot during simulation
    if (iplot && ismembertol(time, plot_times, 0.5*dt)); q_ex = initial_conditions(icase, grid_obj.coords, time, 9.81); plot_solution(q, q_ref, grid_obj, icase, time, gravity, plot_method, 0, 0); frame = getframe(gcf); writeVideo(writerObj, frame); end
    % move step by step
    if (istep); pause(); end
end

%%one-time for writing shores
fileID = fopen('eng_errors_case'+string(icase)+'_pd_nop'+string(nop)+'.dat', 'w+');
fprintf(fileID,'%6d\t%e\t%e\n',[timings;ae1;ae2]);
fclose(fileID);

%%
% close writer
if (iplot);close(writerObj);end