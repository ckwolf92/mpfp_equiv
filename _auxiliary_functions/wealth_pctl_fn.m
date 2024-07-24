function [C_pctl,liqwealth_indic,lambda_pctl] = wealth_pctl_fn(pctl_lb,pctl_ub,liqwealth_dist,lambdafull,C_dist,grid_a);

n_y = size(lambdafull,1);

liqwealth_lb = find(cumsum(liqwealth_dist) >= pctl_lb,1); % smallest that is above the lower bound
liqwealth_ub = max(liqwealth_lb,max(1,sum(cumsum(liqwealth_dist) <= pctl_ub))); % largest that is below the upper bound

liqwealth_lb_indic = (repmat(grid_a,n_y,1) >= grid_a(liqwealth_lb));
liqwealth_ub_indic = (repmat(grid_a,n_y,1) <= grid_a(liqwealth_ub));
liqwealth_indic    = min(liqwealth_lb_indic,liqwealth_ub_indic);

lambda_pctl = liqwealth_indic .* lambdafull;
lambda_pctl = lambda_pctl ./ sum(lambda_pctl(:));
 
C_pctl = sum(sum(C_dist(:,:).*lambda_pctl(:,:)));

end