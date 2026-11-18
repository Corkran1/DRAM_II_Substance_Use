function ss = ssfun_mrsa(theta, data, Y0)
% Sum of squared errors for the six observed patient compartments.
% HCW compartments are latent; this is NOT a six-outcome count likelihood.
penalty = 1e12;
ss = penalty;
try
    [~, Y] = ode45(@(t,y) mrsa_ode(t,y,theta), data.time, Y0);
    if size(Y,1) ~= numel(data.time) || any(~isfinite(Y(:)))
        return;
    end
    observed = [data.S_D(:), data.I_D(:), data.R_D(:), ...
                data.S_F(:), data.I_F(:), data.R_F(:)];
    residual = Y(:,1:6) - observed;
    candidate = sum(residual(:).^2);
    if isfinite(candidate)
        ss = candidate;
    end
catch
    % Invalid parameter proposals receive a large finite penalty.
end
end
