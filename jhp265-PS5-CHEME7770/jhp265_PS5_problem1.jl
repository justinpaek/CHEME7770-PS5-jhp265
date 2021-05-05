# Code for Problem Set 5 Problem 1
# CHEME 7770 Spring 2021

# Problem 1
# Construct the phase portrain for the dynamics of Delta in the two cells

using Makie
using AbstractPlotting
using AbstractPlotting.MakieLayout
AbstractPlotting.inline!(true)

# Model for delta dynamics
# D1 - delta concentration in cell 1
# D2 - delta concentration in cell 2
function delta_ODE(D1, D2)
    D1_dt = 1/(1+10*(D2^2/(0.1+D2^2))^2) - D1
    D2_dt = 1/(1+10*(D1^2/(0.1+D1^2))^2) - D2

    return Point(D1_dt, D2_dt)
end

# Construct the streamplot
plt1 = Scene(resolution =(400,400))
streamplot!(plt1, delta_ODE, 0..4, 0..4, colormap = :plasma, 
    gridsize= (32,32), arrow_size = 0.05)

# Display the plot
display(plt1)

# Save the plot
save("delta_phaseplot.png", plt1)


# using Plots, ImplicitEquations
# Plotting the nullclines
# f(x,y) = 1/(1+10*(y^2/(0.1+y^2))^2) - x
# g(x,y) = 1/(1+10*(x^2/(0.1+x^2))^2) - y
# plot(f==0, xlims=(0,1), ylims=(0,1))
# x = [0:0.01:1;]
# d2_nullcline = 1/(1+10*(x.^2/(0.1.+x.^2)).^2)
# d1_nullcline = (((1/(10*x) - 0.1).^(1/2)*0.1)/(1-(1/(10*x) - 0.1).^(1/2)))^(1/2)

# plot(x, d2_nullcline, xaxis = "D1", yaxis = "D2")