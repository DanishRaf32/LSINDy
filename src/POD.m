function [ xdotr,x0r,Vk ] = POD( X,xdot,x0,Opts )
% Performs POD Reduction of ODE system xdot with initial state x0 using
% training matrix X for creation of reduced basis
%
% xdotr: reduced system
% x0r: projected initial state
% Vk: reduced basis

% Parse options for number of singular values
% Def.NoSV = [];
% if ~exist('Opts','var') || isempty(Opts)
%     Opts = Def;
% else
%     Opts = parseOpts(Opts,Def);
% end

% Create reduced basis and apply to system
[Vk,~] = basisRed(X,Opts.NoSV);
x0r = Vk'*x0;
xdotr = @(t,xr) Vk'*xdot(t,Vk*xr);

end