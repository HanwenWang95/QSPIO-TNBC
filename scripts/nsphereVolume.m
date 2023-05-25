% Reference: Allen, R J et al. “Efficient Generation and Selection of Virtual Populations
%                               in Quantitative Systems Pharmacology Models.”
%             CPT: pharmacometrics & systems pharmacology vol. 5,3 (2016): 140-6.

function [volume] = nsphereVolume(n,radius)

volume=pi^(n/2).*(1/gamma(n/2+1)).*radius.^n;


end
