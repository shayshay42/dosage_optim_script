using DifferentialEquations, Zygote, SciMLSensitivity, Optimization, Optim, OptimizationOptimisers

include("../assets/model.jl")
include("../assets/dosing.jl")
include("../assets/utils.jl")

u0 = [17.7,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
tspan = (0,end_time+7).*hours
p = [ode_params;doses]
prob = ODEProblem(pk_pd!,u0,tspan,p)

sensealg = ForwardDiffSensitivity(convert_tspan=true)
cb2 = SciMLSensitivity.track_callbacks(hit, prob.tspan[1], prob.u0, prob.p, sensealg)

function loss(θ)
    sol = solve(prob, p=[ode_params;θ], callback=cb2, abstol=1e-12, reltol=1e-12, sensealg=sensealg)
    sols = convert(Array,sol)
    c = sols[1,:][end]
    return c + sum(abs, θ)/(1800*(length(θ)))
end

loss(doses)

@time @show Zygote.gradient(loss, doses)

# adtype = Optimization.AutoZygote()
# optf = Optimization.OptimizationFunction((x,_)->loss(x), adtype)
# optprob = Optimization.OptimizationProblem(optf, doses)
# callback, loss_values = create_callback()
# opt = ADAM(0.1)
# @time res = Optimization.solve(optprob, opt, callback = callback, maxiters=50)
# @show res.u
