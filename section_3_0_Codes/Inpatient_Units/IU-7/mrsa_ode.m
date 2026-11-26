function dYdt = mrsa_ode(~, Y, theta)
% Eight-compartment SU/NS patient and HCW MRSA model.
% theta = [beta_D beta_F gamma mu_D mu_F delta sigma ...
%          Lambda_S_D Lambda_I_D Lambda_S_F Lambda_I_F]
S_D = Y(1); I_D = Y(2); R_D = Y(3);
S_F = Y(4); I_F = Y(5); R_F = Y(6);
H_S = Y(7); H_C = Y(8);

beta_D = theta(1); beta_F = theta(2); gamma = theta(3);
mu_D = theta(4); mu_F = theta(5); delta = theta(6);
sigma = theta(7); Lambda_S_D = theta(8); Lambda_I_D = theta(9);
Lambda_S_F = theta(10); Lambda_I_F = theta(11);

N_P = S_D + I_D + R_D + S_F + I_F + R_F;
N_H = H_S + H_C;
% Guard against zero denominators during numerical integration.
hcw_contaminated_fraction = H_C / max(N_H, eps);
infection_pressure = (I_D + I_F) / max(N_P, eps);

dS_D = Lambda_S_D - beta_D * hcw_contaminated_fraction * S_D - mu_D*S_D;
dI_D = Lambda_I_D + beta_D * hcw_contaminated_fraction * S_D - (gamma+mu_D)*I_D;
dR_D = gamma*I_D - mu_D*R_D;
dS_F = Lambda_S_F - beta_F * hcw_contaminated_fraction * S_F - mu_F*S_F;
dI_F = Lambda_I_F + beta_F * hcw_contaminated_fraction * S_F - (gamma+mu_F)*I_F;
dR_F = gamma*I_F - mu_F*R_F;
dH_S = delta*H_C - sigma*H_S*infection_pressure;
dH_C = sigma*H_S*infection_pressure - delta*H_C;

dYdt = [dS_D; dI_D; dR_D; dS_F; dI_F; dR_F; dH_S; dH_C];
end
