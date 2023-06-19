using DifferentialEquations, Zygote, SciMLSensitivity, Optimization, OptimizationOptimisers, OptimizationOptimJL

include("../assets/model.jl")
include("../assets/dosing_unified.jl")
include("../assets/params.jl")
include("../assets/utils.jl")

u0 = [17.7,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
tspan = (0,end_time+7).*hours
p = [ode_params;doses]
prob = ODEProblem(pk_pd!,u0,tspan,p)

sensealg = ForwardDiffSensitivity()#(convert_tspan=true)
# cb2 = SciMLSensitivity.track_callbacks(hit, prob.tspan[1], prob.u0, prob.p, sensealg)

function loss(θ)
    q = [ode_params;θ]
    tmp_prob = remake(prob, p=q)
    cb2 = SciMLSensitivity.track_callbacks(hit, tmp_prob.tspan[1], tmp_prob.u0, tmp_prob.p, sensealg)
    tmp_sol = solve(tmp_prob, p=q, callback=cb2, abstol=1e-12, reltol=1e-12, sensealg=sensealg)
    sols = convert(Array,tmp_sol)
    c = sols[1,:][end]
    return c + sum(abs, θ)/(1800.0*avg_huma_surface_area*(length(θ)))
end

loss(doses)

@time @show Zygote.gradient(loss, doses)

adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x,_)->loss(x), adtype)
optprob = Optimization.OptimizationProblem(optf, doses)
callback, loss_values, iter = create_callback()
opt = ADAMW(0.1)#ADAM(0.1)
@time res = Optimization.solve(optprob, opt, callback = callback, maxiters=50)
@show res.u
