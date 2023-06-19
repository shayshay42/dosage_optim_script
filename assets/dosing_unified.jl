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

tmz_treat_dosetimes = spaced_list(end_treat,1.0,0.0,0.0).*hours
tmz_adjuv_dosetimes = spaced_list(end_time,5.0,23.0,end_treat+28.0).*hours
rg_dosetimes = spaced_list(end_time-1.0,18.0,10.0,0.0).*hours

doses = ones(length(rg_dosetimes)).*dose_amount

inject_times = sort(unique([rg_dosetimes;tmz_treat_dosetimes;tmz_adjuv_dosetimes]));

function affect_dose!(integrator)
    SciMLBase.set_proposed_dt!( integrator, 0.01)
    if integrator.t in tmz_treat_dosetimes
        integrator.u[3] += tmz_treat_dose
    end
    if integrator.t in tmz_adjuv_dosetimes
        integrator.u[3] += tmz_adjuv_dose
    end
    if integrator.t in rg_dosetimes
        hit_rg = integrator.p[length(ode_params)+1:end][findall(x->x==integrator.t,rg_dosetimes)][1]
        integrator.u[6] += relu(hit_rg)
    end
end
hit = PresetTimeCallback(inject_times, affect_dose!);