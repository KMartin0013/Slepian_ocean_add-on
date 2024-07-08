function varargout=SCSpTH(res,buf)  
% XY=MEKONG(res,buf)
% MEKONG(...) % Only makes a plot
%
% Finds the coordinates of Mekong, potentially buffered by some amount.
%
% INPUT:
%
% res      0 The standard, default values
%          N Splined values at N times the resolution
% buf      Distance in degrees that the region outline will be enlarged
%          by BUFFERM, not necessarily integer, possibly negative
%          [default: 0]
%
% OUTPUT:
%
% XY       Closed-curved coordinates of the continent
%
% Last modified by charig-at-princeton.edu, 11/23/2011
% Last modified by fjsimons-at-alum.mit.edu, 11/23/2011

defval('res',10)
defval('buf',0)

% Parameters that make this the region in question
regn='SCSpTH';
% c11=[280 5];
% cmn=[310 -25];
xunt=[];

% This admittedly is a special case
XY=load(fullfile(getenv('IFILES'),'COASTS','SCSpTH'));

XY=XY.XY;
% Do it! Make it, load it, save it

XY=regselect(regn,XY(:,1),XY(:,2),xunt,res,buf);

% The XY could contain not only the study area and the buffering area, 
% but some small regions (e.g., islands). 
% This procedure is to pick up the study area and the buffering area, which
% are properly arranged in clockwise.
if sum(isnan(XY(:,1)))>1
    % ignore small regions
    k=find(isnan(XY(:,1)));
    kk=[0;k;size(XY,1)];
    diff_kk=diff(kk);

    B = sort(diff_kk, 'descend');

    % The outer boundary and inner boundary
    idx1 = find(diff_kk == B(1), 1);
    idx2 = find(diff_kk == B(2), 1);

    XY_1=XY(kk(idx1)+1:kk(idx1+1)-1,:); % No NaN
    XY_2=XY(kk(idx2)+1:kk(idx2+1)-1,:); % No NaN

    [in,on] = inpolygon(XY_2(:,1),XY_2(:,2),XY_1(:,1),XY_1(:,2));

    if all(in) && ~any(on)
        % the first polygon contains the second polygon
        XY_boundary=XY_1; % No NaN
        XY_inner=XY_2; % No NaN
    else
        % the second polygon do not contain the first polygon
        XY_boundary=XY_2; % No NaN
        XY_inner=XY_1; % No NaN
    end

    % Boundary should be definitely clockwise
    [X2,Y2]=poly2cw(XY_boundary(:,1),XY_boundary(:,2));
    XY_boundary=[X2 Y2]; % be careful that polycw could remove NaN
    clear X2 Y2
    % Inner should be definitely counterclockwise
    [X2,Y2]=poly2ccw(XY_inner(:,1),XY_inner(:,2));
    XY_inner=[X2 Y2]; clear X2 Y2

    if buf<0
        XY=XY_inner(end:-1:1,:);
    else
        XY=XY_boundary;
    end

    if XY(end,1)~=XY(1,1) || XY(end,2)~=XY(1,2)
        XY(end+1,:)=XY(1,:);

    end
%     XY=[XY_boundary; NaN, NaN ;XY_inner ];
%     XY=[XY_inner; NaN, NaN; XY_boundary];
end

if nargout==0
  plot(XY(:,1),XY(:,2),'k-'); axis image; grid on
end

% Prepare optional output
varns={XY};
varargout=varns(1:nargout);

