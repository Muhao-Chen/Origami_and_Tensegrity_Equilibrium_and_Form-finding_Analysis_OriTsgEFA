%% Welcome to OriTsgEFA: Origami and Tensegrity Equilibrium and Form-finding Analysis software
% This setup file initializes the environment for OriTsgEFA and should be run only once.

% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.

%% Initialize the MATLAB Workspace
% Clear variables, close open figures, and clear the command window
clear;
close all;
clc;

%% Add Folders to MATLAB Path
% Add function libraries
addpath(genpath('Function_library'));

% Add software verification and examples
addpath(genpath('Software_Verification_and_Examples'));

% Add user guide
addpath(genpath('User_Guide'));

% Add videos folder
addpath(genpath('Videos'));

% Add JOSS paper
addpath(genpath('JOSS_paper'));

%% Display Welcome Message and User Guide Instructions
% Uncomment the next lines if you'd like to automatically open the User Guide PDF.
% cd User_Guide;
% open('User_Guide_Origami and Tensegrity Equilibrium and Form-finding Analysis(OriTsgEFA).pdf');
% cd ..;

% Display a welcome message
disp('Welcome! Please follow the step-by-step instructions from the User Guide.');
