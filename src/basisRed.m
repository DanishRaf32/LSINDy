function [ Vk, Sigmakplusone ] = basisRed( X,k )
% Performs a basis reduction of arbitrary matrix X, if k given according to
% k dominant singular values
%
% Vk: reduced basis
% Sigmakplusone: first neglected singular value

% SVD and plot
[V,Sigma,~] = svd(X,'econ');

% Selection procedure for k
if ~exist('k','var') || isempty(k)
    % Relative Sigma values
    Sigma_rel = Sigma./Sigma(1,1);

    % Plot
    sigmaplot = figure;
    semilogy(diag(Sigma_rel));
    xlim([1 size(Sigma,1)]);
    grid on;
    title('Relative Singular Values of X');

    % User Input
    k = input('Last significant Singular Value: ');
    close(sigmaplot);
end

% Creation of reduced basis
Vk = V(:,1:k);

% First neglected singular value for error bound
if k < size(Sigma,2)
    Sigmakplusone = Sigma(k+1,k+1);
else
    Sigmakplusone = [];
end
end