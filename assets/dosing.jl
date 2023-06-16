#standard of care (SOC):
#20 - 1800mg /m2
#10 days on 28 day cycle

using Random, DifferentialEquations

include("../assets/utils.jl")
include("../assets/params.jl")

#time frames
hours = 24.0
end_time = 28.0*5.0
end_treat = 42.0

#dose amounts
avg_huma_surface_area = 1.7 #m^2
tmz_treat_dose = 75.0*avg_huma_surface_area
tmz_adjuv_dose = 150.0*avg_huma_surface_area
dose_amount = 1800.0*avg_huma_surface_area

#treatment phase (SOC)
tmz_treat_dosetimes = spaced_list(end_treat,1,0,0).*hours
function tmz_treat!(integrator)
    SciMLBase.set_proposed_dt!( integrator, 0.01)
    integrator.u[3] += tmz_treat_dose
end
tmz_treat_hit = PresetTimeCallback(tmz_treat_dosetimes, tmz_treat!);

#adjuvant/maintenance phase (SOC
tmz_adjuv_dosetimes = spaced_list(end_time,5,23,end_treat+28).*hours
function tmz_adjuv!(integrator)
    SciMLBase.set_proposed_dt!( integrator, 0.01)
    integrator.u[3] += tmz_adjuv_dose
end
tmz_adjuv_hit = PresetTimeCallback(tmz_adjuv_dosetimes, tmz_adjuv!);

#drug2 to optimize
rg_dosetimes = spaced_list(end_time-1,18,10,0).*hours
doses = ones(length(rg_dosetimes)).*dose_amount
function rg_dose!(integrator)
    SciMLBase.set_proposed_dt!( integrator, 0.01)
    hit_rg = integrator.p[length(ode_params)+1:end][findall(x->x==integrator.t,rg_dosetimes)][1]
    integrator.u[6] += relu(hit_rg)
end
hit_rg = PresetTimeCallback(rg_dosetimes, rg_dose!);

#alltogether
hit = CallbackSet(tmz_treat_hit, tmz_adjuv_hit, hit_rg);