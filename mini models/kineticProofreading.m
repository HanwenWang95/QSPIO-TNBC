% Mini model of Kinetic proofreading with limited signaling
% from Lever et al. (2014) Analysis

f = figure;
set(f,'units','normalized','outerposition',[0 0 1 1])

N_A = 6.022e23;
N_A = 1;
TCR_tot = 15708/N_A;
k_TCR_p = 1;
phi_TCR = 0.09;
N_TCR = 10;
k_TCR_on = 3.1831e-5*N_A;
P_T = logspace(1,7,5)/N_A;
k_TCR_off = logspace(-1,1,100);
K_D = k_TCR_off./k_TCR_on;

% Varying k_off for different M1p1
subplot(1,3,1)
for i=1:length(P_T)
    C_T = 0.5*(P_T(i) + TCR_tot + K_D - sqrt((P_T(i) + TCR_tot + K_D).^2 -4.*P_T(i).*TCR_tot));
    pTCR = k_TCR_p./(k_TCR_p + k_TCR_off) .* (k_TCR_off./(phi_TCR + k_TCR_off)).^N_TCR.*C_T;
    plot(k_TCR_off,pTCR);hold on; set(gca, 'XScale', 'log'); %set(gca, 'YScale', 'log');
end 
xlabel('$k_{off}$ (1/second)'); ylabel('pTCR');


% Varying M1p1 for different k_off
P_T = logspace(1,7,100)/N_A;
k_TCR_off = logspace(-1,1,6);
K_D = k_TCR_off./k_TCR_on;
subplot(1,3,2)
for i=1:length(k_TCR_off)
    C_T = 0.5*(P_T + TCR_tot + K_D(i) - sqrt((P_T + TCR_tot + K_D(i)).^2 -4.*P_T.*TCR_tot));
    pTCR = k_TCR_p./(k_TCR_p + k_TCR_off(i)) .* (k_TCR_off(i)./(phi_TCR + k_TCR_off(i))).^N_TCR.*C_T;
    plot(P_T,pTCR);hold on; set(gca, 'XScale', 'log');% set(gca, 'YScale', 'log');
end 
% plot(P_T,P_T,'k');
if N_A ==1
    xlabel(['M1p1 ($\#/\mu m^2$)']); ylabel('pTCR');
else
    xlabel(['M1p1 ($mole/\mu m^2$)']); ylabel('pTCR');
end
   
% RO
subplot(2,3,3)
p1_50 = 1e-2/N_A;
for i=1:length(k_TCR_off)
    C_T = 0.5*(P_T + TCR_tot + K_D(i) - sqrt((P_T + TCR_tot + K_D(i)).^2 -4.*P_T.*TCR_tot));
    signal_H = C_T./(C_T+p1_50);
    plot(P_T,signal_H,'k');hold on; set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');ylim([0.5 1]);
end
if N_A ==1
    xlabel(['M1p1 ($\#/\mu m^2$)']); ylabel('Activation signal');
else
    xlabel(['M1p1 ($mole/\mu m^2$)']); ylabel('Activation signal');
end

% KPR
subplot(2,3,6)
p1_50 = 1e-2/N_A;
for i=1:length(k_TCR_off)
    C_T = 0.5*(P_T + TCR_tot + K_D(i) - sqrt((P_T + TCR_tot + K_D(i)).^2 -4.*P_T.*TCR_tot));
    pTCR = k_TCR_p./(k_TCR_p + k_TCR_off(i)) .* (k_TCR_off(i)./(phi_TCR + k_TCR_off(i))).^N_TCR.*C_T;
    signal = pTCR./(pTCR+p1_50);
    plot(P_T,signal);hold on; set(gca, 'XScale', 'log'); set(gca, 'YScale', 'log');ylim([0.5 1]);
end 
if N_A ==1
    xlabel(['M1p1 ($\#/\mu m^2$)']); ylabel('Activation signal');
else
    xlabel(['M1p1 ($mole/\mu m^2$)']); ylabel('Activation signal');
end