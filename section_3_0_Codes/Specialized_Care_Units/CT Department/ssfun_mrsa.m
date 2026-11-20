
function ss = ssfun_mrsa(theta, data, Y0)
    tspan = data.time;
    try
        [~, Y] = ode45(@(t, y) mrsa_ode(t, y, theta), tspan, Y0);
    catch
        ss = 1e12;  % if ODE solver fails, assign large error
        return;
    end

    if size(Y, 1) ~= length(tspan)
        ss = 1e12;  % if output length mismatch, assign large error
        return;
    end

    % Extract model outputs
    S_D = Y(:,1); I_D = Y(:,2); R_D = Y(:,3);
    S_F = Y(:,4); I_F = Y(:,5); R_F = Y(:,6);

    % Sum of squared errors between model and data
    err = (S_D - data.S_D).^2 + (I_D - data.I_D).^2 + (R_D - data.R_D).^2 + ...
          (S_F - data.S_F).^2 + (I_F - data.I_F).^2 + (R_F - data.R_F).^2;

    ss = sum(err);  % total sum of squared errors to minimize
end