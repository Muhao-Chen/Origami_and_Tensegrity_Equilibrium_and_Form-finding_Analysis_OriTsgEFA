function Gp=tenseg_str_gp(gr,C)
% This function is written by Ma, S., Chen, M. and Skelton, R., 2022. TsgFEM: Tensegrity finite element method. Journal of Open Source Software, 7(75), p.3390. 
% The code license is Mozilla Public License, v. 2.0.
% The source code is here:
% https://github.com/Muhao-Chen/Tensegrity_Finite_Element_Method_TsgFEM/tree/main/Function_library/tenseg_str_gp.m
%
% /* This Source Code Form is subject to the terms of the Mozilla Public
% * License, v. 2.0. If a copy of the MPL was not distributed with this
% * file, You can obtain one at http://mozilla.org/MPL/2.0/. 
%
% This function transfer the group infomation, gr, into group matrix Gp.
%
% Inputs:
%	gr: a cell containing vectors, each vector is the number of elements in
%	one group
%	C: the connectivity matrix
% Outputs:
%	Gp: group matrix, ne(number of elements) rows, ng(number of groups) columns

%% 
Gp=eye(size(C,1));
E=eye(size(C,1));
Gp1=[];
%% method 1  keep the origional order of members

if ~isempty(gr)
   for i=1:numel(gr)       % this is to combine members in one group
Gp(:,gr{i}(1))=sum(E(:,gr{i}),2);
   end
      for i=1:numel(gr)         % set duplicate columns into 0
Gp(:,gr{i}(2:end))=0*Gp(:,gr{i}(2:end));
   end
end
% delete zero column
Gp(:,~max(Gp))=[];
%% method 2
% num=[];      %give index for group string number
% for i=1:size(gr,1)
%     num=[num,gr{i}];
% end
% 
% if ~isempty(gr)
% for i=1:size(gr,1)
%     s=sum(E(:,gr{i}),2);
%     Gp1=[Gp1,s];
% end
% Gp(:,num)=[];
% Gp=[Gp1,Gp];
% end
end