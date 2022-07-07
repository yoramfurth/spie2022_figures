function Show_Figure5()
%SHOW_FIGURE5 shows Figure 5 - Evolution of B(p) as a function of inhomogeneity.

%    Copyright 2017-2022 Yoram Furth (yoram.furth@gmail.com)
%    Dept. Electrical & Computer Engineering, BGU Israel.
%    This code is published under GNU GPLv3 license (see license in "LICENSE." file).

% init.
th = 0.01;
K = 1;
ar = 120;
SNR = 100;

% sampling ranges
Kb = [3, 9, 27];
p = 10.^(-4:0.01:-1);

% layout primitives
t = [SNR,0]';  % like snr in SIC data
phi0 = SIM_Get_phi(ar);

% for each Kb calculates B(p)
Ntests = numel(Kb);
B = cell(1,Ntests);
for n = 1:Ntests    
    Kshift = Kb(n)/K;
    [theta, ferr] = fminbnd_inclusive(@(theta)abs(SIM_Calc_Kshift(K,theta,ar) - Kshift), 0, 90);  % finds the corresponding theta 
    assert(ferr<0.1,'cannot find theta(Kshift)');

    [~, AGnom] = SIM_Calc_A(phi0, K, theta, t, th, 'Global');
    [~, ALnom] = SIM_Calc_A(phi0, K, theta, t, th, 'Local');
    B{n} = ALnom(p)./AGnom(p); 
end
B = vertcat(B{:})';

% plots B(p) for every Kb
xAxisRng = [p(1), p(end)];
lineWidth = 1.5;
figure; set(gcf,'Name','B(p;K)');
loglog(p, B, 'LineWidth', lineWidth);
ax = axis; axis([xAxisRng, ax(3:4)])
xlabel('p'); ylabel('B'); grid on;
legend('Data1', 'Data2', 'Data3', 'location', 'NorthEast');
set(gcf,'Position',[671,617,888,323])  
