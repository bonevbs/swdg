% grid1d class carries all the properties of the grid including geometry
% and integration weights. Since this object is being passed all the time
% we force call-by-reference by deriving it from the handle class
classdef grid1d < handle
	properties
        % grid parameters
        nel
        nfa
        ngl
        nop
        npl
        % grid geometry
        x0
        x1       
        % global arrays
        coords
        wq
        faces
        % local arrays
        xgl
        wgl
        % local basis functions
        psi
        dpsi
        % metric terms
        ksi_x
        % stuff for plotting
        plotgrd
        plotvdm
        % stuff for error norms
        normgrd
        normwq
        normvdm
        % filter matrix
        filter
    end
    
    methods
        
        % initializes grid
        function obj = grid1d(nelem, nop, x0, x1)
            
            % copy properties
            obj.x0 = x0;
            obj.x1 = x1;
            
            % generate Legendre-Gauss-Lobatto points
            obj.npl = 20;
            obj.nop = nop;
            obj.ngl = nop+1;
            obj.nel = nelem;
            obj.nfa = nelem+1;
            [obj.xgl, obj.wgl] = lgl_points(obj.ngl);
            % generate the basis functions
            [obj.psi, obj.dpsi] = lagrange_basis(obj.xgl);

            % initialize arrays
            obj.coords  = zeros(obj.ngl, obj.nel);
            obj.wq      = zeros(obj.ngl, obj.nel);
            obj.faces   = zeros(2, obj.nfa);
            obj.plotgrd = zeros(obj.npl, obj.nel);
            obj.normgrd = zeros(2*obj.ngl+1, obj.nel);
            obj.normwq  = zeros(2*obj.ngl+1, obj.nel);

            % create regular grid
            assert(x0 < x1, 'x0 has to be smaller than x1');
            jac     = 0.5*(x1-x0)/obj.nel;
            center  = linspace(x0+jac, x1-jac, obj.nel);
            
            % create integration points for error norms
            [nxgl, nwgl] = lgl_points(2*obj.ngl+1);
            
            % metric terms. Actually quite boring but to keep it consistent
            % with numa
            obj.ksi_x = ones(obj.ngl, obj.nel)./jac;

            % create the coordinate arrays and faces
            for ie=1:nelem
                obj.coords(:,ie)    = center(ie) + jac*obj.xgl;
                obj.wq(:,ie)        = obj.wgl*jac;
                
                obj.plotgrd(:,ie)   = center(ie) + jac*linspace(-1,1,obj.npl);
                
                obj.normgrd(:,ie)   = center(ie) + jac*nxgl;
                obj.normwq(:,ie)    = nwgl*jac;
            end

            % create faces with pointers to left and right elements
            for ifa=1:obj.nfa
                % pointers to adjacent elements
                obj.faces(1,ifa)	= ifa-1;
                obj.faces(2,ifa)    = ifa;
            end
            % create -1 pointers for boundary conditions
            obj.faces(1,1) = -1;
            obj.faces(2,obj.nfa) = -1;
            
            % finally calculate the vandermonde matrix for plotting and
            % error norms
            obj.plotvdm = vandermonde(linspace(-1,1,obj.npl), obj.xgl);
            obj.normvdm = vandermonde(nxgl, obj.xgl);
            
            % calculate filter matrices
            obj.filter = filter_init(obj.ngl, obj.xgl, obj.wq, obj.nel);
        end
    end
end