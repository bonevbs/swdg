% function for plotting.
% currently two display modes are possible.
% direct renders only the values at the nodes
% interp renders the Lagrange polynomials on an equidistant grid
function plot_solution(q,q_ref,grid_obj,icase,time,gravity,plot_method,iwrite,fname)
    nel = grid_obj.nel;
    ngl = grid_obj.ngl;
    npl = grid_obj.npl;
    x0 = grid_obj.x0;
    x1 = grid_obj.x1;
    
    if strcmp(plot_method, 'direct')
        x        = reshape(grid_obj.coords, [1,nel*ngl]);
        [qe,qb]  = initial_conditions(icase, x, time, gravity);
        qe       = (qe +qb)./gravity;
        qh(1,:)  = reshape(q(1,:,:) + q_ref(1,:,:), [1,nel*ngl])./gravity;
        qh(2,:)  = reshape(q(2,:,:), [1,nel*ngl])./gravity;
        qb       = reshape(q_ref, [2,nel*ngl])./gravity;

        xticks = unique(x);
    elseif strcmp(plot_method, 'interp')
        x        = reshape(grid_obj.plotgrd, [1,nel*npl]);
        [qe,qb]  = initial_conditions(icase, x, time, gravity);
        qe       = (qe+qb)./gravity;
        qh       = zeros(2, npl, nel);
        qb       = zeros(2, npl, nel);
        for ie=1:nel
            qh(1,:,ie) = grid_obj.plotvdm*(q(1,:,ie) + q_ref(1,:,ie))'./gravity;
            qh(2,:,ie) = grid_obj.plotvdm*q(2,:,ie)'./gravity;
            qb(1,:,ie) = grid_obj.plotvdm*q_ref(1,:,ie)'./gravity;
        end
        qh  = reshape(qh, [2,nel*npl]);
        qb  = reshape(qb, [2,nel*npl]);
        
        xticks = grid_obj.coords(1,grid_obj.faces(2,1:end-1));
    end

    % first plot waterheight
    subplot(1,2,1)
    plot(x,qh(1,:),'-b',x,qe(1,:),'-g',x,qb(1,:),'-r')
    grid on
    xlim([x0, x1])
    %axis([x0, x1, min(qb), 3])
    set(gca,'XTick', xticks)
    set(gca,'XTickLabel','')
    title('Geopotential height \phi and topography \tau')
    
    % second plot velocity
    subplot(1,2,2)
    plot(x,qh(2,:),'-b',x,qe(2,:),'-g')
    grid on
    xlim([x0, x1])
    set(gca,'XTick', xticks)
    set(gca,'XTickLabel','')
    title('Discharge (\phi u)')
    
    %pause(0.05)
    
    if (iwrite)
        % write to file
        fileID = fopen(fname, 'w+');
        fprintf(fileID,'%e\t%e\t%e\t%e\t%e\t%e\n',[x;qh(1,:);qe(1,:);qb(1,:);qh(2,:);qe(2,:)]);
        fclose(fileID);
    end
end