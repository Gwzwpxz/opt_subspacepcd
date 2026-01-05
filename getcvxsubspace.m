function [d] = getcvxsubspace(X, Ds)
% Get a pre-conditioner that is the linear/affine combination of base pre-conditioners

X = X / trace(X);
k = size(Ds, 2); %#ok

cvx_begin sdp quiet
cvx_precision high
cvx_solver mosek
variable b(k, 1)
variable tau nonnegative
maximize tau
subject to
X - diag(Ds * b) >= 0; %#ok
diag(Ds * b) - X * tau >= 0; %#ok
cvx_end

d = Ds * b;

end % End function