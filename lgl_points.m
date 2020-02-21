% Calculates the Legendre-Gauss-Lobatto points with the corresponding 
% weights. The algorithm which we utilize is described in "Calculation
% of Gauss Quadrature Rules" by Gene H. Golub and John H. Welsh
function [xgl, wgl] = lgl_points(ngl)
    if (ngl < 2)
        error('Cannot generate less than 2 Legendre Gauss Lobatto Points!')
    end
    % For Legendre polynomials we have the recurrence formula
    % k P_k(x) = (2k-1) x P_{k-1}(x) - (k-1) P_{k-2}(x)
    % with P_0(x) = 1 and P_1(x) = x on [-1,1]
    % For Legendre-Gauss points we have the recurrence coefficients
    % a_k = (2k-1)/k
    % b_k = 0
    % c_k = (k-1)/k
    % For (k-1) P'_k(x) = (2k-1) x P'_{k-1}(x) - k P'_{k-2}(x)
    % a_k = (2k-1)/(k-1)
    % b_k = 0
    % c_k = k/(k-1)
    
    % N is the order of the legendre polynomial who's roots we are looking
    % for
    N = ngl-1;
    xgl = zeros(1,ngl);
    wgl = zeros(1,ngl);
    xgl(1) = -1.0; xgl(ngl) = 1.0;
    
    % function handles for upper and lower diagonal of the T matrix
    % see Golub et al. for details
    a = @(k) (2*k-1)/(k-1);
    b = @(k) 0;
    c = @(k) k/(k-1);
    ld = @(k) c(k)/a(k);
    ud = @(k) 1/a(k);
    
    % solve eigenvalue problem x*p(x) = T*p(x) + (1/a_N)p_N(x)e_N
    % if p_N(t_j) = 0, t_j has to be an eigenvalue of t_j*p(t_j) = T*p(t_j)
    T = diag(arrayfun(ld,3:N),-1) + diag(arrayfun(ud,2:N-1),1);
    [~,D] = eig(T);
    [D,~] = sort(diag(D));  
    xgl(2:ngl-1) = D';
    
    % calculate the weights using 2/(n*(n-1)*P_{n-1}^2(x))
    P = legendre(ngl-1, xgl);
    P = P(1, :).*P(1, :);
    wgl = 2/(ngl*(ngl-1))./P;

end