% Main purpose: use slepian functions to calculate the sea level anomaly in
%   a regional ocean.
% This code is based on the Free Software from the Simons Laboratories
%   (https://geoweb.princeton.edu/people/simons/software.html).
%   with some changes for oceanic application (e.g., GAD, GIA, IB correction).
% It is strongly recommended to familiarize yourself with the original code
%   before using this add-on code.
%
% Required software:
%   slepian_alpha (by csdms-contrib)
%   slepian_bravo (by csdms-contrib)
%   slepian_delta (by csdms-contrib)
%   slepian_zero (actually only guyotphysics.m) (by csdms-contrib)
%   m_map (by https://www.eoas.ubc.ca/~rich/map.html)

% Main Reference:
% Ref1: Harig, C., & Simons, F. J. (2012). Mapping GreenlandAs mass loss in space and time. Proceedings of the National Academy of Sciences, 109(49), 19934-19937.
% Ref2: Zhongtian, Ma et al., A Novel Slepian Approach for Determining Mass-term Sea Level from GRACE over the
%           South China Sea
