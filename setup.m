% SETUP.m
% Author: Edgar Hernandez
% Date: June 13, 2025
% ------------
% Purpose:
%   Script to setup project for developmental use, may need to run script
%   at each MATLAB instance. 

clear
clc
setupPaths

function setupPaths
trunkPath = pwd;
addpath(genpath(pwd))
cd (trunkPath)
end