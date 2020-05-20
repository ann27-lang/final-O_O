% Problem 2
clear all
close all

%Parameters
ke = 1e-4;
kestar = 5e-3;
kf = 5.14e-21;
kr = 2.5e-2;
kdeg = 8e-4;
Vs = 18;
gamma = 1e2;
DL = 1e-10;
q = 1e3;
nc = 3e8;

z = linspace(0,5);

for i = 1:length(z)
    Kss = kestar*kf/(ke*(kr+kestar));

    Shz = (gamma*z(i)^2/DL)^(1/3);

    Lc = (q*nc)/Shz;

    secquan = Kss*Lc/((z(i)/DL)+Kss*Lc);

    firquan = (1/kestar + 1/kdeg);

    Rtot(i) = firquan*secquan*Vs;
end

figure
plot(z,Rtot)
xlabel("z position")
ylabel("R_t_o_t^* concentration")
title("z vs. R_t_o_t^*")
