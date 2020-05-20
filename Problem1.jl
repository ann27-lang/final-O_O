# CHEME 5440
# Problem 1
include("PhasePortraitV2.jl")

function odefun(d1, d2)

    u = @. 1/(1+10*(d2^2/(0.1+d2^2))^2) - d1
    v = @. 1/(1+10*(d1^2/(0.1+d1^2))^2) - d2

    return(u,v)
end


tspan = (0.0, 5.0)

xlimca = (0, 1, 50)
ylimcr = (0, 1, 50)

# To show that in the long-time limit, system will settle into a steady state
xo= ([0.8, 0.7],)

phaseplot(odefun,xlimca,ylimcr,clines=true,xinit=xo,t=tspan,norm=true,scale=0.5)
