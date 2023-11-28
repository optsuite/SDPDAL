function ATy_h = theta_applyAT_full(Ref_E, y)

ATy = full(mex_theta_applyAT(Ref_E, y));
ATy = (ATy + ATy');
ATy_h = @(R) R * ATy;

end
