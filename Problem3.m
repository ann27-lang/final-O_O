% Problem 3

%Parameters
u_value = [0.2 0.204 0.418 0.720 0.965 0.988 0.990];
KX = [0.245 0.270 0.528 0.862 1.107 1.197 1.197];
KL = 1550.388;
tau = 0.0825;
protelnc = 0.055;
ribconc = 18.254;

denominator = 0.464;

for i = 1:7
    numerator = protelnc*ribconc*KX(i)/(tau*KL)*3600;
    
    pistar(i) = numerator/denominator*u_value(i);
end

figure
plot(u_value,pistar)
title("u_i vs p_i^* profile")
xlabel("Promoter model u_i (AU)")
ylabel("p_i^* concentration(nmol/gDW)")

%try part c when Kp > 1, let Kp = 10
Kp = 10;
for i = 1:7
    numeratortry = protelnc*ribconc*KX(i)/(tau*KL/Kp)*3600;
    
    pistartry(i) = numeratortry/denominator*u_value(i);
end

figure
plot(u_value,pistartry)
title("u_i vs p_i^* profile with polysome amplification")
xlabel("Promoter model u_i (AU)")
ylabel("p_i^* concentration(nmol/gDW)")