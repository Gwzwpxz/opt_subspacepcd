function [prob] = optprecond(prob)
% Compute the optimal diagonal preconditioner for X

ptype = upper(prob.ptype);

if prob.cond > 1e+10
    warning("The matrix is likely to be rank deficient.");
end % End if

if ptype == 'R'
    prob.E = getcvxdiag(prob.M, ptype);
    prob.E = prob.E.^(-0.5);
    prob.pX = prob.X * diag(prob.E);
elseif ptype == 'L'
    prob.D = getcvxdiag(prob.X, ptype);
    prob.D = sqrt(prob.D);
    prob.pX = diag(prob.D) * prob.X;
elseif ptype == 'T'
    [D1, D2] = gettwosidedprecond(prob.X);
    prob.D = sqrt(diag(D2));
    prob.E = sqrt(diag(D1).^-1);
    prob.pX = diag(prob.D) * prob.X * diag(prob.E);
else
    if prob.params.method == "SINF"
        warning off;
        fprintf("Using the semi-infinite approach for SDP\n");
        info = semiinfsdp(prob.M, prob.params.subspace);
        warning on;
        prob.E = info.y;
    else
        prob.E = getcvxsubspace(prob.M, prob.params.subspace);
    end % End if
    prob.E = prob.E.^(-0.5);
    prob.pX = prob.X * diag(prob.E);
end % End if

% Upperbound may not be improved.
% Use cond instead of condest to get true condition number
prob.optcond = condest(prob.pX' * prob.pX);

end % End function