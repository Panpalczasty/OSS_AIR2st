clc
clear all
close all

%% Description

% 29.05.22
% MPC wrapper for solver
% Cleanup and preparation for 

% 01.05.22
% This file implements optimization for changes in input value
% at a discrete raster - with zoh discrete changes
% input changest with the same sampling freq as sim
% cost fun is a combination of energy minimalization, trajectory deviation
% and final position deviation

%% setup
mode = "MPC";
%mode = "std";

curve = "circle";
horizon = 1;
step = .1;
time = 8;
sampling_freq = 10;
center = [0 -1];
radius = 0.2;
R = 1e-2;
Wx = 1e1;
Wv = 1e-1;
Wstat = 1e2;

tname = sprintf("%s_T%d_h%d_s%.1f_%s_R%.2f", mode, time, horizon, step, curve, radius);

%% get model parameters
params = getDefaultParams(time,sampling_freq);
[~, params] = getTrajectory(curve, params, center, radius);
params = getWeights(params,R,Wx,Wv,Wstat);

%% optimize

if mode=="std"
    horizon = params.T;
    step = params.T;
end

[uopt, ~] = MPC(params, horizon, step);
[x,t,p,hu,pos] = getStats(params, uopt);

[u_qa, e_qa, pos_qa] = getQA(uopt, x, params.xref, params)

%% plot

plotError(1, t, x, params.xref, params, tname);
%plotAdjoint(1,t,p,params, tname);
plotAngles(2,t,x,uopt,params, tname);
plotPos(3, pos, params.posref, params, tname);
plotSwitch(5,t,uopt,hu,params, tname);

animate(6,pos, params, tname);
