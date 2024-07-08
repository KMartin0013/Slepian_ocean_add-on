function lmcosig=plm_Gausssmooth(lmcosi,radius,a)
% plm_Gauss=plm_Gaussian_smoothing(lmcosi,order,radius)
%
% Converts the real spherical harmonic coefficients expressing a geoPOTENTIAL
% field into a smoothed version with Gaussian filter of the radius R
%
% INPUT:
%
% lmcosi     [l m Ccos Csin] degrees, order, coefficients, as in PLM2XYZ
%            Units are properly those of potential, J/kg or m^2/s^2
% radius     The averaging radius, in km [default: 500]
% a          Reference radius [default: a_EGM96]
%
% OUPUT:
%
% lmcosig    [l m Ccos Csin] degrees, order, smoothed coefficients of 
%            output in the standard units of what's requested [J/kg, m/s^2, m]
%
% See Also: 
% 
% Wahr et al. 1998, for Gaussian Kernel
%
% Last modified by zhongtian.ma@connect.polyu.hk, 06/05/2024


% compute Gaussian kernel W
defval('a',fralmanac('a_EGM96','Earth'))
defval('radius',500)

degree_max=lmcosi(end,1);

W=ones(1,degree_max+1); 
b=log(2)/(1-cos(radius*1000/a));
W(1)=1;
W(2)=((1+exp(-2*b))/(1-exp(-2*b))-1/b);
for I=2:degree_max
   W(I+1)=-(2*I-1)*W(I)/b+W(I-1);
end

lmcosW=[ones(size(lmcosi,1),2) repmat(W(lmcosi(:,1)+1)',1,2)];

lmcosig=lmcosi.*lmcosW;

end


