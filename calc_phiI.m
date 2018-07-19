% Author: Susanne Schindler
% Date: 9.6.2017
% Copyright: Susanne Schindler (susanne.schindler2@web.de)

% script returns involvement of intruder
function phiI = calc_phiI(sizeDiff_breederIntruder)

epsilon = get_params; 

phiI = 1/(1+sizeDiff_breederIntruder/epsilon);