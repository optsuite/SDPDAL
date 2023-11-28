function ATy_h = theta_applyAT(Ref_E, y)

ATy = mex_theta_applyAT(Ref_E, y);
ATy = ATy + ATy';
ATy_h = @(R) R * ATy;

end