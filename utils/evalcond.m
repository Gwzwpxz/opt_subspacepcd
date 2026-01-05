function [condval] = evalcond(A, D)

if size(A, 1) > 5000 || min(D) <= 0
    condval = -1.0;
    return
end % End if

A = full(A);
Dhalf = diag(sqrt(1./D));
condval = cond(Dhalf * A * Dhalf);

end % End function