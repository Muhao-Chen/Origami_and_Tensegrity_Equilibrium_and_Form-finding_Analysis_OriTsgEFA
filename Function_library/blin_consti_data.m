function [data_b1,data_b2]=blin_consti_data(Eb,sigma_b)
% This function is written by Ma, S., Chen, M. and Skelton, R., 2022. TsgFEM: Tensegrity finite element method. Journal of Open Source Software, 7(75), p.3390.
% The code license is Mozilla Public License, v. 2.0.
% The source code is here:
% https://github.com/Muhao-Chen/Tensegrity_Finite_Element_Method_TsgFEM/tree/main/Function_library/bin_consti_data.m
%
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/.
%
% This function output the bilinear constitutive law given Young's modulus and
% yielding stress,
%
% Inputs:
%   Eb: tangent Young's modulus in zero strain
%   sigma_b: yielding stress(the first turing point of stress-strain curve)
% Outputs:
%	data_b1: strain information( tension and compression)
%	data_b2: stress information( tension and compression)
%% material
strain_b1=[sigma_b/Eb,2];
stress_b1=[sigma_b,sigma_b];
[data_b1,data_b2,Eb,sigma_b]=point2consti_data(strain_b1,stress_b1);
end

