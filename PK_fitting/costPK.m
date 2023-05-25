function cost = costPK(opt, model, dose_schedule_qsp, tindex, expData)

Qp    = opt(1);
k_cl  = opt(2);
K_P   = opt(3);
Vc    = opt(4);
k_cln = opt(5);
Kcl   = opt(6);

model.parameter(1).Value = exp(Qp);
model.parameter(5).Value = exp(k_cl);
model.parameter(9).Value = exp(K_P);
model.compartment(1).Capacity = exp(Vc);
model.parameter(6).Value = exp(k_cln);
model.parameter(7).Value = exp(Kcl);

MW = 1.463e8; % mg/mole durvalumab

simData = sbiosimulate(model,[],[],dose_schedule_qsp);

K_B = model.parameter(8).Value;
diff = log(expData) - log(simData.Data(:,1).*MW./K_B);
cost = norm(diff)^2;

end
