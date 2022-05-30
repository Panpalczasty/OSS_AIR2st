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
params.T = 30;
params = getDefaultParams(params);

% Weight manipulation
W = 1*[1; 1; 1; 1]*1;
Wstat = 1e2*[1; 1; 1; 1];
params.R = 1;

% update params with appropriate horizon
horizon = 1;
step = 0.2;

% Define trajectory here - as x = f(t), y = g(t)
% Parametric curve definition
t = linspace(0, params.T, params.T*params.sfreq);
T = params.T;
w = 2*pi/T;

sf = params.sfreq
pos = zeros(length(t),2);

% circle
pos(:,1) = 1 + 0.5*cos(w*t);
pos(:,2) = 0.5 * sin(w*t);

% S shaped curve
%pos(1:sf*T/2,1) = 0.3 + 0.3*cos(w*t(1:sf*T/2));
%pos(1:sf*T/2,2) = 1+0.3*sin(w*t(1:sf*T/2));
%pos(sf*T/2+1:end,1) = -0.3 + 0.3*cos(w*t(1:sf*T/2));
%pos(sf*T/2+1:end,2) = 1-0.3*sin(w*t(1:sf*T/2));

% sharp angle
% pos(1:sf*T/2,1) = 0.7 + t(1:sf*T/2)/T;
% pos(1:sf*T/2,2) = t(1:sf*T/2)/T;
% pos(sf*T/2+1:end,1) = 1.2;
% pos(sf*T/2+1:end,2) = -t(1:sf*T/2)/T + 0.5;

%% get model parameters

params.W = eye(length(params.x0)).*W;
params.Wstat = eye(length(params.x0)).*Wstat;

% exponential accel
% pos(:,1) = ones(length(t),1);
% pos(:,2) = -1+0.1*exp(3*t/T);

params.posref = pos;

[params.t,params.xref] = inverseKin(pos, params.T, "below", params.J1, params.J2);

params.x0 = [params.xref(1,1:2) 0 0 0];
params.xf = [params.xref(end,1:2) 0 0];
params.xref(end,:) = params.xf;
%params.x0 = [params.xref(1,:) 0];
%params.xf = params.xref(end,:);

%params.u = uopt;
params.u = zeros((length(t)),2);


%% optimize

[uopt, ~] = MPC(params, horizon, step);


%% get stats

disp(costFun(params,uopt))
 
% solve with optimized control 
[t, x] = qsolve45(uopt, params.T, params);

% antisolve with optimized control
pf =  - (params.Wstat)*(x(end,1:end-1) - params.xf)';
pf = [pf;0;0];
p = qantisolve45(pf, t, x, uopt, params);

pos = simpleKin(x, params.J1, params.J2);
Hu = getHuvec(t, x, p, uopt, params);

%% plots, plots, plots

plotAdjoint(1,t,p,params);
plotAngles(2,t,x,uopt,params);
plotPos(3, pos, params.posref, params);
plotSwitch(5,t,uopt,Hu,params);

animate(6,pos, params);
