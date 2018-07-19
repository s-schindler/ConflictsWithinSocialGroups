% Author: Susanne Schindler
% Date: 9.6.2017
% Copyright: Susanne Schindler (susanne.schindler2@web.de)

% script calculates the optimal involvement of helper (full-sib,
% half-sib staying, half-sib defeating, unrelated), breeder
% (same-sex, other sex than intruder) for a given infanticide risk,
% size difference between intruder and same-sex breeder, and
% efficiancy of the helper to chase away the intuder
% 
% code can be run both with and  without input parameters: 
% input parameters: k1 ... how fast involvement needed for defense
%                   increases with size difference between same-sex 
%                    
%                   k2 ... inefficency of helpers
%                   
%                   k3 ... remaining brood dependency
%                   
%                   sizeDiff_breederIntruder ... factor that
%                   intruder is larger than breeder
%
%                   probOfInfanticide...risk that intruder kills
%                   infants upon take-over
%
%                   adultSurvival ... extrinsic (non-contest
%                   related survival probability) 
function [opt_phi_H_defeat,opt_phi_BS_wHdefeat,opt_phi_BO_wHdefeat,opt_phi_H_full,opt_phi_BS_wHfull,opt_phi_BO_wHfull,opt_phi_H_stay,opt_phi_BS_wHstay,opt_phi_BO_wHstay,opt_phi_BO_wHunrel,phi_I,rel_payoff_BS_wHstay,rel_payoff_BO_wHstay,rel_payoff_Hstay,rel_payoff_BS_wHdefeat,rel_payoff_BO_wHdefeat,rel_payoff_Hdefeat,rel_payoff_BS_wHfull,rel_payoff_BO_wHfull,rel_payoff_Hfull,rel_payoff_BS_wHunrel,rel_payoff_BO_wHunrel,rel_nf_payoff_BS_wHstay,rel_nf_payoff_BO_wHstay,rel_nf_payoff_Hstay,rel_nf_payoff_BS_wHdefeat,rel_nf_payoff_BO_wHdefeat,rel_nf_payoff_Hdefeat,rel_nf_payoff_BS_wHfull,rel_nf_payoff_BO_wHfull,rel_nf_payoff_Hfull,rel_nf_payoff_BS_wHunrel,rel_nf_payoff_BO_wHunrel,payoff_won_H_stay,payoff_lost_H_stay,payoff_won_BS,payoff_lost_BS,payoff_won_H_full,payoff_lost_H_full,payoff_won_H_defeat,payoff_lost_H_defeat] = calc_involvement(k1,k2,k3,sizeDiff_breederIntruder,probOfInfanticide,adultSurvival)

% parameters
if nargin ~= 6
  sizeDiff_breederIntruder = 0:0.005:0.25;
  probOfInfanticide = 0:0.02:1;
  adultSurvival = 0.5:0.1:1;
  k1 = [1 5 10 20];  % steepness of curve giving minimal involvement needed
  k2 = [1 1.5 2]; % efficiancy of helper
  k3 = [0 5 10];  % length of remaining dependence of current brood
end

noSizeDiffs = length(sizeDiff_breederIntruder);
noInfanticideVals = length(probOfInfanticide);
noSurvivalRates = length(adultSurvival);
noK1 = length(k1);
noK2 = length(k2);
noK3 = length(k3);

% initialization
phi_H = zeros(noK1,noK2,noSizeDiffs);
phi_BO = zeros(noK1,noK2,noSizeDiffs);
phi_BS = zeros(noK1,noK2,noSizeDiffs);
phi_I = zeros(1,noSizeDiffs);

payoff_won_BS = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
payoff_won_BO = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
payoff_won_H_stay =  zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
payoff_won_H_defeat = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
payoff_won_H_full =  zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);

payoff_lost_BS = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
payoff_lost_BO = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
payoff_lost_H_stay = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
payoff_lost_H_defeat = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
payoff_lost_H_full =zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);

for s=1:noSizeDiffs
  % involvement of intruder
  phi_I(s) = calc_phiI(sizeDiff_breederIntruder(s));
  
  for i1=1:noK1
    % minimal involvement neccessary to defeat intruder
    needed_I(s) = calc_phiMin(k1(i1),sizeDiff_breederIntruder(s));

    for i2=1:noK2
      % calc involvement of helper needed for winning
      phi_H(i1,i2,s) = min(1,k2(i2)*needed_I(s));
      % calc involvement of breeders needed for winning
      phi_BO(i1,i2,s) = (needed_I(s) - phi_H(i1,i2,s)/k2(i2))/2;
      phi_BS(i1,i2,s) = phi_BO(i1,i2,s);
    
      for i3=1:noK3
	for surv=1:noSurvivalRates
	  % calc payoff in case of winning (accounting for needed
	  % involvement)
	  [payoff_won_BS(i1,i2,i3,s,:,surv),payoff_won_BO(i1,i2,i3,s,:,surv),payoff_won_H_stay(i1,i2,i3,s,:,surv),payoff_won_H_defeat(i1,i2,i3,s,:,surv),payoff_won_H_full(i1,i2,i3,s,:,surv)] =  calc_payoff_won(phi_BS(i1,i2,s),phi_BO(i1,i2,s),phi_H(i1,i2,s),k3(i3),adultSurvival(surv));

	  for m=1:noInfanticideVals
	    % calc payoff in case of loosing (assuming no involvement of helper
	    % and other-sex breeder, but with involvement of intruder)
	    [payoff_lost_BS(i1,i2,i3,s,m,surv),payoff_lost_BO(i1,i2,i3,s,m,surv),payoff_lost_H_stay(i1,i2,i3,s,m,surv),payoff_lost_H_defeat(i1,i2,i3,s,m,surv),payoff_lost_H_full(i1,i2,i3,s,m,surv)] =  calc_payoff_lost(phi_I(s),probOfInfanticide(m),k3(i3),adultSurvival(surv));
	  end
	end
      end
    end
  end
end

% check whether outcome of winning (1) is bigger than loosing
% (0) and save preferred outcome
preferredOutcome_BS = (payoff_won_BS > payoff_lost_BS);
preferredOutcome_BO = (payoff_won_BO > payoff_lost_BO);
preferredOutcome_H_stay = (payoff_won_H_stay > payoff_lost_H_stay);
preferredOutcome_H_defeat = (payoff_won_H_defeat > payoff_lost_H_defeat);
preferredOutcome_H_full = (payoff_won_H_full > payoff_lost_H_full);

% use preferences of contest outcome to determine optimal
% involvement of helper
% initialisation
opt_phi_H_stay = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
opt_phi_BS_wHstay = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
opt_phi_BO_wHstay = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
opt_phi_H_defeat = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
opt_phi_BS_wHdefeat = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
opt_phi_BO_wHdefeat = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
opt_phi_H_full = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
opt_phi_BS_wHfull = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
opt_phi_BO_wHfull = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
opt_phi_BO_wHunrel = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);

rel_payoff_BS_wHstay = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
rel_payoff_BO_wHstay = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
rel_payoff_Hstay = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
rel_payoff_BS_wHdefeat = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
rel_payoff_BO_wHdefeat = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
rel_payoff_Hdefeat = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
rel_payoff_BS_wHfull = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
rel_payoff_BO_wHfull = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
rel_payoff_Hfull = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
rel_payoff_BS_wHunrel = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
rel_payoff_BO_wHunrel = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
rel_nf_payoff_BS_wHstay = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
rel_nf_payoff_BO_wHstay = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
rel_nf_payoff_Hstay = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
rel_nf_payoff_BS_wHdefeat = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
rel_nf_payoff_BO_wHdefeat = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
rel_nf_payoff_Hdefeat = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
rel_nf_payoff_BS_wHfull = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
rel_nf_payoff_BO_wHfull = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
rel_nf_payoff_Hfull = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
rel_nf_payoff_BS_wHunrel = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);
rel_nf_payoff_BO_wHunrel = zeros(noK1,noK2,noK3,noSizeDiffs,noInfanticideVals,noSurvivalRates);

% calculation
for i1=1:noK1
  for i2=1:noK2
    for i3=1:noK3
      for s=1:noSizeDiffs
	for m=1:noInfanticideVals
	  for l=1:noSurvivalRates
	    
	    %% calculate optimal involvements
	    % 1st case helper only related to staying breeder
	    if preferredOutcome_H_stay(i1,i2,i3,s,m,l)
	      opt_phi_H_stay(i1,i2,i3,s,m,l) = phi_H(i1,i2,s);
	      opt_phi_BS_wHstay(i1,i2,i3,s,m,l) = phi_BS(i1,i2,s);
	      opt_phi_BO_wHstay(i1,i2,i3,s,m,l) = phi_BO(i1,i2,s);
	  
	      % realised payoff (record two matrices for cases that
	      % BS steps in or not)
	      [rel_payoff_BS_wHstay(i1,i2,i3,s,m,l),rel_payoff_BO_wHstay(i1,i2,i3,s,m,l),rel_payoff_Hstay(i1,i2,i3,s,m,l),p1,p1] = calc_payoff_won(opt_phi_BS_wHstay(i1,i2,i3,s,m,l),opt_phi_BO_wHstay(i1,i2,i3,s,m,l),opt_phi_H_stay(i1,i2,i3,s,m,l),k3(i3),adultSurvival(l));
	      [rel_nf_payoff_BS_wHstay(i1,i2,i3,s,m,l),rel_nf_payoff_BO_wHstay(i1,i2,i3,s,m,l),rel_nf_payoff_Hstay(i1,i2,i3,s,m,l),p1,p1] = calc_payoff_won(opt_phi_BS_wHstay(i1,i2,i3,s,m,l),opt_phi_BO_wHstay(i1,i2,i3,s,m,l),opt_phi_H_stay(i1,i2,i3,s,m,l),k3(i3),adultSurvival(l));
	    
	    else
	      % if helper prefers to breed with intruder BS might take on the
	      % defense by itself, even if close to 1
	      opt_phi_BS_wHstay(i1,i2,i3,s,m,l) = calc_phiMin(k1(i1),sizeDiff_breederIntruder(s));
	  
	      %% calculate realised payoffs when BS steps in for H
	      %(i.e. conflict won)
	      [rel_payoff_BS_wHstay(i1,i2,i3,s,m,l),rel_payoff_BO_wHstay(i1,i2,i3,s,m,l),rel_payoff_Hstay(i1,i2,i3,s,m,l),p1,p1] = calc_payoff_won(opt_phi_BS_wHstay(i1,i2,i3,s,m,l),opt_phi_BO_wHstay(i1,i2,i3,s,m,l),opt_phi_H_stay(i1,i2,i3,s,m,l),k3(i3),adultSurvival(l));

	      %% calculate realised payoffs when BS doesn't step in for H
	      %(i.e. conflict lost, note: no-one gets involved)
	      [rel_nf_payoff_BS_wHstay(i1,i2,i3,s,m,l),rel_nf_payoff_BO_wHstay(i1,i2,i3,s,m,l),rel_nf_payoff_Hstay(i1,i2,i3,s,m,l),p1,p1] = calc_payoff_lost(phi_I(s),probOfInfanticide(m),k3(i3),adultSurvival(l));
	  
	    end

	    % 2nd case helper only related to challenged (would-be defeated) breeder
	    if preferredOutcome_H_defeat(i1,i2,i3,s,m,l)
	      opt_phi_H_defeat(i1,i2,i3,s,m,l) = phi_H(i1,i2,s);
	      opt_phi_BS_wHdefeat(i1,i2,i3,s,m,l) = phi_BS(i1,i2,s);
	      opt_phi_BO_wHdefeat(i1,i2,i3,s,m,l) = phi_BO(i1,i2,s);

	      % realised payoff (record two matrices for cases that
	      % BS steps in or not)
	      [rel_payoff_BS_wHdefeat(i1,i2,i3,s,m,l),rel_payoff_BO_wHdefeat(i1,i2,i3,s,m,l),p1,rel_payoff_Hdefeat(i1,i2,i3,s,m,l),p1] = calc_payoff_won(opt_phi_BS_wHdefeat(i1,i2,i3,s,m,l),opt_phi_BO_wHdefeat(i1,i2,i3,s,m,l),opt_phi_H_defeat(i1,i2,i3,s,m,l),k3(i3),adultSurvival(l));
	      [rel_nf_payoff_BS_wHdefeat(i1,i2,i3,s,m,l),rel_nf_payoff_BO_wHdefeat(i1,i2,i3,s,m,l),p1,rel_nf_payoff_Hdefeat(i1,i2,i3,s,m,l),p1] = calc_payoff_won(opt_phi_BS_wHdefeat(i1,i2,i3,s,m,l),opt_phi_BO_wHdefeat(i1,i2,i3,s,m,l),opt_phi_H_defeat(i1,i2,i3,s,m,l),k3(i3),adultSurvival(l));
	  
	    else
	      % if helper prefers to breed with intruder BS might take on the
	      % defense by itself, even if close to 1
	      opt_phi_BS_wHdefeat(i1,i2,i3,s,m,l) = calc_phiMin(k1(i1),sizeDiff_breederIntruder(s));

	      %% calculate realised payoffs when BS steps in for H
	      %(i.e. conflict won)
	      [rel_payoff_BS_wHdefeat(i1,i2,i3,s,m,l),rel_payoff_BO_wHdefeat(i1,i2,i3,s,m,l),p1,rel_payoff_Hdefeat(i1,i2,i3,s,m,l),p1] = calc_payoff_won(opt_phi_BS_wHdefeat(i1,i2,i3,s,m,l),opt_phi_BO_wHdefeat(i1,i2,i3,s,m,l),opt_phi_H_defeat(i1,i2,i3,s,m,l),k3(i3),adultSurvival(l));
	      %% calculate realised payoffs when BS doesn't step in for H
	      %(i.e. conflict lost, note: no-one gets involved)
	      [rel_nf_payoff_BS_wHdefeat(i1,i2,i3,s,m,l),rel_nf_payoff_BO_wHdefeat(i1,i2,i3,s,m,l),p1,rel_nf_payoff_Hdefeat(i1,i2,i3,s,m,l),p1] = calc_payoff_lost(phi_I(s),probOfInfanticide(m),k3(i3),adultSurvival(l));

	    end

	    % 3rd case helper related to both breeders
	    if preferredOutcome_H_full(i1,i2,i3,s,m,l)
	      opt_phi_H_full(i1,i2,i3,s,m,l) = phi_H(i1,i2,s);
	      opt_phi_BS_wHfull(i1,i2,i3,s,m,l) = phi_BS(i1,i2,s);
	      opt_phi_BO_wHfull(i1,i2,i3,s,m,l) = phi_BO(i1,i2,s);

	      % realised payoff (record two matrices for cases that
	      % BS steps in or not)
	      [rel_payoff_BS_wHfull(i1,i2,i3,s,m,l),rel_payoff_BO_wHfull(i1,i2,i3,s,m,l),p1,p1,rel_payoff_Hfull(i1,i2,i3,s,m,l)] = calc_payoff_won(opt_phi_BS_wHfull(i1,i2,i3,s,m,l),opt_phi_BO_wHfull(i1,i2,i3,s,m,l),opt_phi_H_full(i1,i2,i3,s,m,l),k3(i3),adultSurvival(l));	    
	      [rel_nf_payoff_BS_wHfull(i1,i2,i3,s,m,l),rel_nf_payoff_BO_wHfull(i1,i2,i3,s,m,l),p1,p1,rel_nf_payoff_Hfull(i1,i2,i3,s,m,l)] = calc_payoff_won(opt_phi_BS_wHfull(i1,i2,i3,s,m,l),opt_phi_BO_wHfull(i1,i2,i3,s,m,l),opt_phi_H_full(i1,i2,i3,s,m,l),k3(i3),adultSurvival(l));	 
	    else
	      % if helper prefers to breed with intruder BS might take on the
	      % defense by itself, even if close to 1
	      opt_phi_BS_wHfull(i1,i2,i3,s,m,l) = calc_phiMin(k1(i1),sizeDiff_breederIntruder(s));
	    
	      %% calculate realised payoffs when BS steps in for H
	      %(i.e. conflict won)
	      [rel_payoff_BS_wHfull(i1,i2,i3,s,m,l),rel_payoff_BO_wHfull(i1,i2,i3,s,m,l),p1,p1,rel_payoff_Hfull(i1,i2,i3,s,m,l)] = calc_payoff_won(opt_phi_BS_wHfull(i1,i2,i3,s,m,l),opt_phi_BO_wHfull(i1,i2,i3,s,m,l),opt_phi_H_full(i1,i2,i3,s,m,l),k3(i3),adultSurvival(l));	    
	      %% calculate realised payoffs when BS doesn't step in for H
	      %(i.e. conflict lost, note: no-one gets involved)
	      [rel_payoff_BS_wHfull(i1,i2,i3,s,m,l),rel_payoff_BO_wHfull(i1,i2,i3,s,m,l),p1,p1,rel_payoff_Hfull(i1,i2,i3,s,m,l)] = calc_payoff_lost(phi_I(s),probOfInfanticide(m),k3(i3),adultSurvival(l));	    
	    end

	    % 4th case helper unrelated to either of the breeders
	    opt_phi_BS_wHunrel(i1,i2,i3,s,m,l) = calc_phiMin(k1(i1),sizeDiff_breederIntruder(s));
	
	    % realised payoff (record two matrices for cases that
	    % BS steps in or not)
	    [rel_payoff_BS_wHunrel(i1,i2,i3,s,m,l),rel_payoff_BO_wHunrel(i1,i2,i3,s,m,l),p1,p1,p1] = calc_payoff_won(opt_phi_BS_wHunrel(i1,i2,i3,s,m,l),0,0,k3(i3),adultSurvival(l));
	    [rel_nf_payoff_BS_wHunrel(i1,i2,i3,s,m,l),rel_nf_payoff_BO_wHunrel(i1,i2,i3,s,m,l),p1,p1,p1] = calc_payoff_lost(phi_I(s),probOfInfanticide(m),k3(i3),adultSurvival(l));


	  end
	end
      end
    end
  end
end

end

function [p_BS,p_BO,p_HSstay,p_HSdefeat,p_HFS] = calc_payoff_won(inv_BS,inv_BO,inv_H,k3,adultSurvival)
  % % parameters
  % [p1,adultSurvival,offspringSurvival,numberEggs] = get_params;
  % relatedness coefficients
  r_BS = 0.5;
  r_BO = 0.5;
  r_HSstay = 0.25;
  r_HSdefeat = 0.25;
  r_HFS = 0.5;
  % % calculate surviving brood size
  % broodSize = offspringSurvival*numberEggs;
  % calculate payoff from current reproduction
  currentBroodFullyDep = 1-mean([inv_BS inv_BO inv_H]);
  % account for length of remaining dependence of young
  currentBrood = calc_factorDependence(k3,currentBroodFullyDep);
  % calculate payoff from future reproduction
  futureBrood = (1-inv_BS)*(1-inv_BO)*adultSurvival^2;
  % calculate payoff from current and future reproduction
  payoff =  currentBrood + futureBrood;
  % weigh payoff with relatedness and brood size
  p_BS = r_BS*payoff;
  p_BO = r_BO*payoff;
  p_HSstay =r_HSstay*payoff;
  p_HSdefeat =r_HSdefeat*payoff;
  p_HFS =r_HFS*payoff;
end

function [p_BS,p_BO,p_HSstay,p_HSdefeat,p_HFS] = calc_payoff_lost(inv_I,probOfInfanticide,k3,adultSurvival)
  % % parameters
  % [p1,adultSurvival,offspringSurvival,numberEggs] = get_params;
  % relatedness coefficients
  r_BS = 0.5;
  r_BO = 0.5;
  r_HSstay = 0.25;
  r_HSdefeat = 0.25;
  r_HFS = 0.5;
  % relatedness coefficients after takeover
  r_new_BO = 0.5;
  r_new_HSstay = 0.25;
  r_new_HSdefeat = 0;
  r_new_HFS = 0.25;
  % % calculate size of surviving brood
  % broodSize = offspringSurvival*numberEggs;
  % calculate payoff from current reproduction and account for length of remaining dependence of young
  effectCurrentBrood = calc_factorDependence(k3,2/3*(1-probOfInfanticide));
  % calculate effect on payoff from future reproduction
  effectFutureBrood = (1-inv_I)*adultSurvival^2;
  % weigh payoffs with relatedness and brood size
  p_BS = r_BS*effectCurrentBrood;
  p_BO = r_BO*effectCurrentBrood + r_new_BO*effectFutureBrood;
  p_HSstay = r_HSstay*effectCurrentBrood + r_new_HSstay*effectFutureBrood;
  p_HSdefeat = r_HSdefeat*effectCurrentBrood + r_new_HSdefeat*effectFutureBrood;
  p_HFS = r_HFS*effectCurrentBrood + r_new_HFS*effectFutureBrood;
end