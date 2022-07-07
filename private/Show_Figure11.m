function Show_Figure11()
%SHOW_FIGURE11 shows Figure 11 - Evolution of B(p) as a function of different thresholds.

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% init.
K = 1;
Kshift = 9;
ar = 120;
SNR = 100;

% sampling ranges
thAll = [0.01, 0.001, 0.0001];
p = 10.^(-4:0.01:-1);

% layout primitives
t = [SNR,0]';  % like snr in SIC data
phi0 = SIM_Get_phi(ar);

% calculates the desirable opening angle (theta)
[theta, ferr] = fminbnd_inclusive(@(theta)abs(SIM_Calc_Kshift(K,theta,ar)-Kshift), 0, 90);
assert(ferr<0.1,'cannot find theta(Kshift)');

% for each "th" calculates B(p)
Ntests = numel(thAll);
B = cell(1,Ntests);
for n = 1:Ntests
    th = thAll(n);
    [~, AGnom] = SIM_Calc_A(phi0, K, theta, t, th, 'Global');
    [~, ALnom] = SIM_Calc_A(phi0, K, theta, t, th, 'Local');
    B{n} = ALnom(p)./AGnom(p); 
end
B = vertcat(B{:})';

% plots B(p) for every "th"
xAxisRng = [p(1), p(end)];
lineWidth = 1.5;
figure; set(gcf,'Name','B(p;th)');
loglog(p, B, 'LineWidth', lineWidth);
ax = axis; axis([xAxisRng, ax(3:4)])
xlabel('p'); ylabel('B'); grid on;
legend('th=0.01','th=0.001','th=0.0001','location','NorthEast');
set(gcf,'Position',[671,617,888,323])  
