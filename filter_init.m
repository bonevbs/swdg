% assemble filter matrix for each element
function filter = filter_init(ngl, xgl, wq, nel)
    filter = zeros(ngl,ngl,nel);
    fbasis = zeros(ngl,ngl);

    nop = ngl-1;
    order = 10;
    alpha = 30;

    % build the weight matrix
    weights = diag(exp(-alpha.*((0:nop)./nop).^order));
    %weights = zeros(ngl,ngl);
    %weights(1,1) = 1.0;

    % build filter basis per element
    for ie=1:nel
        % construct modal filter basis functions on nodes
        for i=1:ngl
            % construct monome
            fbasis(:,i) = xgl.^(i-1);

            % perform numerical Gram-Schmidt
            for j=i-1:-1:1
                fbasis(:,i) = fbasis(:,i) - sum(fbasis(:,j).*fbasis(:,i).*wq(:,ie))*fbasis(:,j);
            end
            fbasis(:,i) = fbasis(:,i)/sqrt(sum(fbasis(:,i).*wq(:,ie).*fbasis(:,i)));
            
            %plot(xgl, fbasis(i,:));
            %pause();
        end

        filter(:,:,ie) = fbasis*weights/fbasis;
        %filter(:,:,ie) = fbasis\fbasis;
    end
end