% Problem 4
clear all
close all

%part a, for W1
r_hat = 3.003;

%Parameters
kcat = 0.4*3600;
E1 = 0.12;
F6P = 0.1*1000;
ATP = 2.3*1000;
KF6P = 0.11*1000;
KATP = 0.42*1000;

r1 = kcat*E1*(F6P/(KF6P+F6P))*(ATP/(KATP+ATP));

syms W1
eqn1 = r_hat == (W1/(1+W1))*r1;

S1 = double(solve(eqn1,W1))

%part a, for W2
r_hat = 60.306;

syms W2
eqn2 = r_hat == ((S1+W2)/(1+S1*W2))*r1;

S2 = double(solve(eqn2,W2))

% part b
data = table2array(readtable('Data-3-5-AMP.txt'));
act_val = data(:,1)*1000;
overall = data(:,2)/r1;
g = [];
fun =@(z,act_val) (act_val./z(2)).^z(1)./(1+(act_val./z(2)).^z(1));
var = [1.85 .2];
g = lsqcurvefit(fun,var,act_val,overall);
g = real(g);
hill = (act_val./g(2)).^g(1)./(1+(act_val./g(2)).^g(1));

% to smooth...
plzwork = linspace(act_val(1,1),act_val(6,1));
smoothhill = interp1(act_val,hill,plzwork,'pchip');

figure
plot(act_val,overall,'ob')
hold on
plot(plzwork,smoothhill,'-r')
title("f_i function")
hold off

%part c
f = (act_val./g(2)).^g(1)./(1+(act_val./g(2)).^g(1));
u_fun = (S1+S2*f)/(1+S1+S2*f);
r_hat = r1*u_fun;

% to smooth...
plzwork = linspace(act_val(1,1),act_val(6,1));
smoothr_hat = interp1(act_val,r_hat(:,6),plzwork,'pchip');

figure
plot(plzwork/1000,smoothr_hat, '-r')
hold on
errorbar(data(:,1),data(:,2),data(:,3),'ob')
title("cAMP concentration vs. PFK Rate")
xlabel("cAMP concentration (mM)")
ylabel("PFK rate (\muM/hr)")
legend("Estimated Overall","Measured")
hold off
