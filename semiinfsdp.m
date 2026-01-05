function [info] = semiinfsdp(X, Ds, tol, maxit, verbose, A1, A2)
% Apply semi-infinite programming to solve
% maximize              tau
% subject to
%                                  X >= diag(Ds * alpha)
%         diag(Ds * alpha) - X * tau >= 0
%                         Ds * alpha >= 0
if nargin < 8
    A1 = [];
    A2 = [];
end % End if

if nargin < 4
    tol = 1e-18;
    maxit = 100;
    verbose = false;
end % End if

[n, p] = size(Ds);

tlp = 0.0;
teig = 0.0;

if n >= 1e+05
    neigiter = 20;
else
    neigiter = 30;
end % End if

e1s = ones(maxit, 1);
e2s = ones(maxit, 1);
condvals = ones(maxit, 1);

if ~isempty(A1) && ~isempty(A2)
    if verbose
        fprintf("Warm starting from existing cuts\n");
    end % End if
    k1 = size(A1, 1);
    k2 = size(A2, 1);
else
    % Build initial constraints
    if verbose
        fprintf("Starting from scratch\n");
    end % End if
    k1 = 10;
    k2 = 10;
    A1 = randn(k1, n);
    A2 = randn(k2, n);
end % End if

B1 = zeros(k1, p + 1);
B2 = zeros(k2, p + 1);
rhs1 = zeros(k1, 1);
rhs2 = zeros(k2, 1);

for i = 1:k1
    B1(i, end) = -A1(i, :) * X * A1(i, :)';
    for j = 1:p
        B1(i, j) = A1(i, :) * (Ds(:, j) .* A1(i, :)');
    end % End for
end % End for

for i = 1:k2
    rhs2(i) = A2(i, :) * X * A2(i, :)';
    for j = 1:p
        B2(i, j) = A2(i, :) * (Ds(:, j) .* A2(i, :)');
    end % End for
end % End for

e1 = -1.0;
e2 = -1.0;

for iter = 1:maxit
    
    rng(iter * n);
    [alpha, tau, subt, pi] = getgrbsol(B1, B2, Ds, rhs1, rhs2, 2);
    
    if e1 >= -tol && e2 >= -tol && min(y) > 0.0
        break;
    end % End if
    
    if iter >= maxit - 1
        break;
    end % End if
    
    tlp = tlp + subt;
    
    y = Ds * alpha;
    S1 = sparse(1:n, 1:n, y) - X * tau;
    S2 = X - sparse(1:n, 1:n, y);
    
    tic;
    [v1, e1] = eigs(S1, 1, 'smallestreal', 'FailureTreatment', 'keep',...
        'MaxIterations', neigiter, 'Tolerance', 1e-15, 'StartVector', randn(n, 1), 'IsFunctionSymmetric', 1);
    [v2, e2] = eigs(S2, 1, 'smallestreal', 'FailureTreatment', 'keep',...
        'MaxIterations', neigiter, 'Tolerance', 1e-15, 'StartVector', randn(n, 1), 'IsFunctionSymmetric', 1);
    teig = teig + toc;
    
    if e1 > -tol
        e1 = tol / 10;
    end % End if
    
    if e2 > -tol
        e2 = tol / 10;
    end % End if
    
    e1s(iter) = e1;
    e2s(iter) = e2;
    
    if e1 < -tol
        A1 = [A1; v1'];
        B1 = [B1;
            zeros(1, p), -v1' * X * v1];
        rhs1 = [rhs1; 0.0];
        for j = 1:p
            B1(end, j) = A1(end, :) * (Ds(:, j) .* A1(end, :)');
        end % End for
    end % End if
    
    if e2 < -tol
        A2 = [A2; v2'];
        B2 = [B2;
            zeros(1, p), 0.0];
        rhs2 = [rhs2; v2' * X * v2];
        for j = 1:p
            B2(end, j) = A2(end, :) * (Ds(:, j) .* A2(end, :)');
        end % End for
    end % End if
    
    y = y / norm(y);
    
    try
        if n <= 1000
            if verbose
                condvals(iter) = evalcond(X, y);
                fprintf("%3d E1: %+5.1e E2: %+5.1e *Kappa* %10.8e tLp: %f tEig: %f \n", iter, e1, e2,...
                    condvals(iter), tlp, teig);
            end % End if
        else
            if verbose
                fprintf("%3d E1: %+5.1e E2: %+5.1e Est. Kappa %10.8e tLp: %f tEig: %f \n", iter, e1, e2,...
                    1 / tau, tlp, teig);
            end % End if
        end % End if
    catch
        
    end % End try
    
end % End for

[~, tau, ~, ~] = getgrbsol(B1, B2, Ds, rhs1, rhs2, 2);

info.niter = iter;
info.alpha = alpha;
info.y = y;
info.e1s = e1s;
info.e2s = e2s;
info.condvals = condvals;
info.kappaest = 1 / tau;
info.pi = pi;

% Warm start
% if scratch
info.A1 = A1;
info.A2 = A2;
% end % End if

end % End function