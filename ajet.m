function cmap = ajet(n)
%AJET produces a variant of the colormap 'jet'
%
% USAGE:    cmap = ajet 
%           cmap = ajet(n)
%
% produces an alternative colormap to function jet.m. The output cmap
% is a matrix of size [n,3] whose columns contain the R,G,B color levels.
% The integer n defaults to 64.
%
% Jet involves a blue to black transition at one end of the map, and a
% red to black transition at the other end. This maps the two extremes
% look very similar (nerly black). THis function has red transition to 
% magenta instead.
% 
% See also functions jet.m, colormap.m, and caxis.m

% version 2.0         Michael Bevis        22 May 2014
if nargin==0
    n=64;
    cm=jet(n+20);
    cm(1:10,:)=[];
    cm(end-9:end,:)=[];
elseif numel(n)~=1 | n<1 | rem(n,1)~=0
    error('input argument n must be a single positive integer')
else
    m=fix(10*n/6);
    cm=jet(n+2*m);
    cm(1:m,:)=[];
    cm(end-m+1:end,:)=[];
end
cmap=cm;
