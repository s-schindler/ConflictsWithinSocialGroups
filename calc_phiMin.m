% Author: Susanne Schindler
% Date: 21.3.2017
% Copyright: Susanne Schindler (susanne.schindler2@web.de)

% script calculates the amount of group involvement needed to
% defeat intruder
function phiMin = calc_phiMin(k,sizeDiff_breederIntruder)

% parameters
epsilon = get_params;
	       
phiMin = exp(log(epsilon)*exp(-k*sizeDiff_breederIntruder));
