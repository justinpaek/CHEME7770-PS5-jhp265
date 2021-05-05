# Code for Problem Set 5 Problem 3
# CHEME 7770 Spring 2021

# Problem 3c
# 
    
using DifferentialEquations
using Plots
gr(show = true)

# Model parameters
L = 0       # start with no ligand to solve for initial conditions
Vmax_R = 1                          # uM/s 
k_bp = 1                            # 1/s
k_b = 0                             # 1/s
alpha1_minus = L/(1+L)              # 1/s
alpha1_plus = 1/(1+L)               # 1/s
beta_1 = (2.5*L)/(1+L)              # 1/s
a_bp = 0.1*1000                     # 1/s*1/uM
a_b = 1*1000                        # 1/s*1/uM
d_bp = 0.01                         # 1/s 
d_b = 1                             # 1/s
k1_minus = 1                        # 1/s
k1_plus = 10                        # 1/s


# Two-state model for adaptation in bacterial chemotaxis
function chemotaxis!(du, u, p, t)
    du[1] = -Vmax_R + k_bp*u[7] + k_b*u[6]
    du[2] = Vmax_R + alpha1_minus*u[3] - alpha1_plus*u[2] +beta_1*(u[7]+u[6])
    du[3] = alpha1_plus*u[2] - alpha1_minus*u[3] - a_bp*u[3]*u[5] - a_b*u[3]*u[4] + d_bp*u[7] + d_b*u[6]
    du[4] = d_b*u[6] - a_b*u[3]*u[4] + k_b*u[6] + beta_1*u[6] + k1_minus*u[5] - k1_plus*u[4]
    du[5] = d_bp*u[7] - a_bp*u[3]*u[5] + k_bp*u[7] + beta_1*u[7] + k1_plus*u[4] - k1_minus*u[5]
    du[6] = a_b*u[3]*u[4] - u[6]*(d_b + k_b + beta_1)
    du[7] = a_bp*u[3]*u[5] - u[7]*(d_bp + k_bp + beta_1)
end

# Solve the steady-state problem to get initial conditions
# Initial guesses
E0₀ = 1.0
E1₀ = 1.0 
E1star₀ = 1.0
B₀ = 1.0
Bp₀ = 1.0
E1starB₀ = 1.0
E1starBp₀ = 1.0
u0 = [E0₀;E1₀;E1star₀;B₀;Bp₀;E1starB₀;E1starBp₀]
tspan = (0.0,100.0)

# Solve the steady state problem
ssprob = SteadyStateProblem(chemotaxis!, u0, tspan)
sssol = solve(ssprob)

# Solve the model for different ligand concentrations

#For L = 0.01:
L = 0.01
alpha1_minus = L/(1+L)              # 1/s
alpha1_plus = 1/(1+L)               # 1/s
beta_1 = 2.5*L/(1+L)                # 1/s
tspan = (0.0,50.0)
prob1 = ODEProblem(chemotaxis!, sssol, tspan)
sol1 = solve(prob1)
activity1 = sol1/1.00285

#For L = 0.1:
L = 0.1
alpha1_minus = L/(1+L)              # 1/s
alpha1_plus = 1/(1+L)               # 1/s
beta_1 = 2.5*L/(1+L)                # 1/s
tspan = (0.0,100.0)
prob2 = ODEProblem(chemotaxis!, sssol, tspan)
sol2 = solve(prob2)
activity2 = sol2/1.00285

#For L = 1:
L = 1
alpha1_minus = L/(1+L)              # 1/s
alpha1_plus = 1/(1+L)               # 1/s
beta_1 = 2.5*L/(1+L)                # 1/s
tspan = (0.0,10000.0)
prob3 = ODEProblem(chemotaxis!, sssol, tspan)
sol3 = solve(prob3)
activity3 = sol3/1.00285    

#For L = 10:
L = 10
alpha1_minus = L/(1+L)              # 1/s
alpha1_plus = 1/(1+L)               # 1/s
beta_1 = 2.5*L/(1+L)                # 1/s
tspan = (0.0,10000.0)
prob4 = ODEProblem(chemotaxis!, sssol, tspan)
sol4 = solve(prob4)
activity4 = sol3/1.00285 

# Plot the results 
using Plots

plt1 = plot(sol1,vars=(0,3), title = "L=0.01", xaxis="time", xlims=(0,50), yaxis = "A/Ast", ylims=(0,1.5), legend = false)
display(plt1)
save("P3_plt1.png", plt1)

plt2 = plot(sol2,vars=(0,3), title = "L=0.1", xaxis="time", xlims=(0,100), yaxis = "A/Ast", ylims=(0,1.5), legend = false)
display(plt2)
save("P3_plt2.png", plt1)

plt3 = plot(sol3,vars=(0,3), title = "L=1", xaxis="time", xlims=(0,10000), yaxis = "A/Ast", ylims=(0,1.5), legend = false)
display(plt3)
save("P3_plt3.png", plt1)

plt4 = plot(sol4,vars=(0,3), title = "L=10", xaxis="time", xlims=(0,10000), yaxis = "A/Ast", ylims=(0,1.5), legend = false)
display(plt4)
save("P3_plt4.png", plt1)