function dYdt = mrsa_ode(t, Y, theta)
    % State variables
    S_D = Y(1); I_D = Y(2); R_D = Y(3);
    S_F = Y(4); I_F = Y(5); R_F = Y(6);
    H_S = Y(7); H_C = Y(8);

    % Parameters unpacking
    beta_D     = theta(1);
    beta_F     = theta(2);
    gamma      = theta(3);
    mu_D       = theta(4);
    mu_F       = theta(5);
    delta      = theta(6);
    sigma      = theta(7);
    Lambda_S_D = theta(8);
    Lambda_I_D = theta(9);
    Lambda_S_F = theta(10);
    Lambda_I_F = theta(11);

    % Total population (drug users + drug free)
    N_P = S_D + I_D + R_D + S_F + I_F + R_F;

    % Equations for Drug-Using Patients
    dS_D = Lambda_S_D - beta_D*(H_C/(H_S+H_C)) * S_D - mu_D * S_D;
    dI_D = Lambda_I_D + beta_D*(H_C/(H_S+H_C))*S_D - gamma * I_D - mu_D * I_D;
    dR_D = gamma * I_D - mu_D * R_D;

    % Equations for Drug-Free Patients
    dS_F = Lambda_S_F - beta_F*(H_C/(H_S+H_C))* S_F - mu_F * S_F;
    dI_F = Lambda_I_F + beta_F * (H_C/(H_S+H_C))*S_F - gamma * I_F - mu_F * I_F;
    dR_F = gamma * I_F - mu_F * R_F;

    % Hospital compartments (example transitions)
   

    dH_S = delta * H_C - sigma * H_S * infection_pressure;
    dH_C = sigma * H_S * infection_pressure - delta * H_C;

    % Return derivatives vector
    dYdt = [dS_D; dI_D; dR_D; dS_F; dI_F; dR_F; dH_S; dH_C];
end