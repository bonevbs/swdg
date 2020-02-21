N = 10;

nel = [1,2,5,10,20,50,100,200,500,1000];
nel = nel(1:N);
err1 = zeros(1,N);
err2 = zeros(1,N);

tol = 1E-4;


for nop=1:4
    for i=1:N
        fprintf('Running\tnop=%d\tnel=%d \n', nop, nel(i))
        % manufactured solution, no dry areas, no bottom topography analytic forcing term
        [err1(i),err2(i)] = swdg(200,0,0.01,0.00001,-0.5,0.5,nop,nel(i),'strong','periodic',0,0,tol,0);
        % dambreak on wet domain, filtering on, poslimit off
        %[err1(i),err2(i)] = swdg(3,0,0.1,0.0001,-2,2,nop,nel(i),'strong','exact',1,0,tol,0);
        % dambreak on dry domain, filtering on, poslimit on
        %tol = 1E-4/(nel(i));
        %[err1(i),err2(i)] = swdg(4,0.1,1.0,0.00005,-2,2,nop,nel(i),'strong','exact',1,1,tol,0);
        %[err1(i),err2(i)] = swdg(4,0.1,1.0,0.0001,0,0.5,nop,nel(i),'strong','exact',1,0,tol,0);
        % oscillating lake, filtering on, poslimit on 
        %tol = 1E-2/(nel(i));
        %[err1(i),err2(i)] = swdg(6,0,1,0.0001,-1.5,1.5,nop,nel(i),'strong','exact',1,1,tol,0);
        % oscillating lake, no dry areas, filtering off, poslimit off 
        %[err1(i),err2(i)] = swdg(6,0,0.5,0.0001,-0.5,0.5,nop,nel(i),'strong','exact',0,0,tol,0);
        fprintf('err1: %e err2: %e \n', err1(i), err2(i))
    end

    % write to file
    fileID = fopen('conv_l2wet_case6_nop'+string(nop)+'.dat', 'w+');
    fprintf(fileID,'%6d\t%e\t%e\n',[nel;err1;err2]);
    fclose(fileID);
    
    % plot
%     slope1 = nel.^(-1);
%     slope2 = nel.^(-2);
%     loglog(nel,err2)
%     hold on
%     loglog(nel,slope1)
%     hold on
%     loglog(nel,slope2)
end