clc
clear all
close all

%% Desc

% 29.05.22
% TODO - MPC wrapper for solver
% TODO - 

% 01.05.22
% This file implements optimization for changes in input value
% at a discrete raster - with zoh discrete changes
% input changest with the same sampling freq as sim
% cost fun is a combination of energy minimalization, trajectory deviation
% and final position deviation

%% setup

% Define sampling & final simulation time
params.sfreq = 10;
params.T = 60;
params = getDefaultParams(params);

% Weight manipulation
W = 1*[1; 1; 1; 1];
Wstat = 1e2*[1; 1; 1; 1];
params.R = 1e-2;

% update params with appropriate horizon
% horizon = 10;
% step = 1;

% uncomment to reduce to last task
horizon = params.T;
step = params.T;

%% get model parameters

params.W = eye(length(params.x0)).*W;
params.Wstat = eye(length(params.x0)).*Wstat;

% exponential accel
% pos(:,1) = ones(length(t),1);
% pos(:,2) = -1+0.1*exp(3*t/T);

pos =getMPCTrajectory('wobble', params.sfreq, params.T, [1 0], 0.3);
params.posref = pos;
[params.t,params.xref] = inverseKin(pos, params.T, "below", params.J1, params.J2);

params.x0 = [params.xref(1,1:2) 0 0 0];
params.xf = [params.xref(end,1:2) 0 0];
params.xref = [params.xref; params.xf];
%params.x0 = [params.xref(1,:) 0];
%params.xf = params.xref(end,:);

%% optimize

[uopt, ~] = MPC(params, horizon, step);


%% get stats

% disp(costFun(params,uopt))
 
% solve with optimized control 
[t, x] = qsolve45(uopt, params.T, params);

% antisolve with optimized control
pf =  - (params.Wstat)*(x(end,1:end-1) - params.xf)';
pf = [pf;0;0];
p = qantisolve45(pf, t, x, uopt, params);

pos = simpleKin(x, params.J1, params.J2);
Hu = getHuvec(t, x, p, uopt, params);

%% checc gradient

%checkGradient(uopt ,params, 0.001);
% is ok

%% plots, plots, plots

plotAdjoint(1,t,p,params);
plotAngles(2,t,x,uopt,params);
plotPos(3, pos, params.posref, params);
plotSwitch(5,t,uopt,Hu,params);

animate(6,pos, params);
