%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% explicitMatrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function y = explicitMatrix(A,x,mode)
%
% Function that can be used to create anonymous function A(x,mode) for 
% explicit matrices.
% Copyright (C) 2007-2008 Elaine Hale, Wotao Yin and Yin Zhang
% FPC is distributed under the GNU GPL, see README.txt and COPYING.txt.

function y = explicitMatrix(A,x,mode)


switch mode
    case 1
        y = A*x;
    case 2
        y = A'*x;
    otherwise
        error('Unknown mode passed to explicitMatrix in fpc.m');
end

end % explicitMatrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%