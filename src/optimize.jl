using Random
using DifferentialEquations, SciMLSensitivity
using Zygote, Optimization, OptimizationOptimisers, OptimizationOptimJL

#loads the PK/PD model
include("../assets/model.jl")
#loads the dosing schedule and amount
include("../assets/dosing_unified.jl")
#loads the PK/PD parameters
include("../assets/params.jl")
#loads the utility functions
include("../assets/utils.jl")

#setting up ODE problem
u0 = [17.7,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
tspan = (0,end_time+7).*hours
p = [ode_params;doses]
prob = ODEProblem(pk_pd!,u0,tspan,p)

#sensitivity method
sensealg = ForwardDiffSensitivity()#(convert_tspan=true)
# cb2 = SciMLSensitivity.track_callbacks(hit, prob.tspan[1], prob.u0, prob.p, sensealg)

anim = Animation()

#loss function balancing tumor growth and drug toxicity
function loss(θ)
    q = [ode_params;θ]
    tmp_prob = remake(prob, p=q)
    tmp_sol = solve(tmp_prob, p=q, callback=hit, sensealg=sensealg)#, abstol=1e-12, reltol=1e-12)
    Zygote.@ignore frame(anim, plotter(tmp_sol))
    sols = convert(Array,tmp_sol)
    c = sols[1,:][end]
    return c + sum(abs2, θ)/(dose_amount*dose_amount*(length(θ)))
end

#check:
# loss(doses)
# @time @show Zygote.gradient(loss, doses)

#random initial dose amounts
rng = Random.default_rng()
Random.seed!(rng, 69)
dose_init = rand(length(rg_dosetimes)).*dose_amount

# gradient descent optimization of the dose amounts
adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x,_)->loss(x), adtype)
optprob = Optimization.OptimizationProblem(optf, dose_init)
#callbacks for observing the optimization
callback, loss_values, iter = create_callback()
opt = Optim.LBFGS()#ADAM(0.1)
@time res = Optimization.solve(optprob, opt, callback=callback, maxiters=50)
@show res.u

#plot result
q = [ode_params;res.u]
tmp_prob = remake(prob, p=q)
tmp_sol = solve(tmp_prob, p=q, callback=hit)
display(plotter(tmp_sol))

#save animation
gif(anim, "BFGSoptim_L2.gif", fps = 10)

#save loss values
plot(loss_values, label="loss")
