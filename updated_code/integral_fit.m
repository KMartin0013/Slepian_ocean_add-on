function [ESTtotal_CC]=integral_fit(N,functionintegrals,ESTsignaldelta)
% The last part of the function slept2resid to multiple the Slepian 
% coefficient by the integral (found using INTEGRATEBASIS) of the corresponding
% function over the region
%
% INPUT:
%
% N          Number of largest eigenfunctions in which to expand.  By default
%             rounds to the Shannon number.
% functionintegrals   This is a vector of the integrals of the Slepian
%                      functions, up to the Shannon number.  With this we
%                      can multiply by the data or ESTSignal to get the
%                      changes in total mass over time for each function.
%                      This also has the conversion from kg to Gt already
%                      in it, so don't do that again.
% ESTsignaldelta   Coefficients with reference to some mean
%
% OUTPUT:
%
% ESTtotal_CC   Total sum of the data for each Slepian
%                 coefficient evaluated at those months in the same format 
%
% Last modified by zhongtian.ma@connect.polyu.hk, 07/07/2024

% calculate for individual or all slepian coefficients
if length(N)==1
    loop_range=1:N; % The total number
else
    loop_range=N; % selective number
end

for i=loop_range

    % Since Int should have units of (fn * m^2), need to go from fractional
    % sphere area to real area.  If the fn is surface density, this output is
    % in kilograms.  Then change the units from kg to Gt in METRIC tons

    %         [eigfunINT_CC] = integratebasis(CC{i},TH,N);

    %     eigfunINT = eigfunINT*4*pi*6370000^2/10^3/10^9;
    %     functionintegrals = eigfunINT;

    % Now multiply by the appropriate slepcoffs to get the months
    % This becomes alpha by months
    %functimeseries=repmat(eigfunINT',1,nmonths).*sleptdelta(:,1:N)';
    %functimeseries = sleptdelta(:,1:N)';

    % Here do the total sum of the data
    %  Original code for integrating of all slepian functions
    % >total=eigfunINT*sleptdelta(:,1:N)'; (Line 526 in in slept2resid)
    %  Updated code for each slepian function
    % because >functionintegrals=eigfunINT; (Line 518 in in slept2resid)

    ESTtotal_CC(i,:)=functionintegrals(i)*ESTsignaldelta(:,i)';


end



end
