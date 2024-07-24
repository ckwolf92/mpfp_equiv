function roots = fastroot(xgrid,fx)
%fast linear interpolation root finding
%(=one Newton step at largest negative function value)
%   stripped down version of interp1 that accepts multiple inputs (max 3)
%   that are interpolated over the same grids x & xi
xgrid=xgrid(:);
fx=reshape(fx,[numel(xgrid),numel(fx)/numel(xgrid)]);

dxgrid=diff(xgrid);
dfx=diff(fx);
idx=ones(1,numel(fx)/numel(xgrid));

% Make use of the fact that the difference equation is monotonically
% increasing in m
idx_min=(fx(1,:)>0); %Corner solutions left (if no solution x* to f(x)=0 exists)
idx_max=(fx(end,:)<0); %Corner solutions right (if no solution x* to f(x)=0 exists)
index=find(and(not(idx_min),not(idx_max))); %interior solutions (if solution x* to f(x)=0 exists)

% Find index of two gridpoints where sign of fx changes from positive to negative,
[~,idx(index)]=max(diff(sign(fx(:,index))));

aux_index  = (0:numel(fx)/numel(xgrid)-1)*numel(xgrid); %aux for linear indexes
aux_index2 = (0:numel(fx)/numel(xgrid)-1)*(numel(xgrid)-1);
fxx  = fx(idx+aux_index);
xl   = xgrid(idx)';
dx   = dxgrid(idx)';
dfxx = dfx(idx+aux_index2);
% Because function is piecewise linear in gridpoints, one newton step is
% enough to find the solution
roots = xl-fxx.*dx./dfxx;

roots(idx_min)=xgrid(1); %constrained choice
roots(idx_max)=xgrid(end); % no-extrapolation
end