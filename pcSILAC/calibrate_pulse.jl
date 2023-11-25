# Running Bayesian model for experimental design calibration of the PCP dataset
#
# Author: astukalov
###############################################################################
#addprocs(3)
@everywhere (myid()==1) || (ENV["MKL_NUM_THREADS"]=1) # disable MKL threading on workers

using StatsBase, Interpolations, BlackBoxOptim, Distributions, DataFrames, Query,
      ParallelDataTransfer, RData, CodecZlib, JSON

println("Pulse/Chase-SILAC experiments calibration")

job_args = ["phubel_pulse", "20171025", "test", "0"];
isdefined(:job_args) || (job_args = ARGS);

project_id, data_version, calib_version, job_name, job_id = job_args
job_id = parse(Int, job_id)
println("Project '$(project_id)' job: name=$(job_name) dataset version=$(data_version) calibration version=$(calib_version) id=$(job_id)")

@passobj 1 workers() base_scripts_path
include_success = false
while !include_success
    try
        @everywhere include(joinpath(base_scripts_path,"adhoc/phubel_pulse/pulse_dynamics/pulse_dynamics.jl"))
        include_success = true
    catch e
        include_success = false
    end
end

data_path = joinpath(base_data_path, project_id)
scratch_path = joinpath(base_scratch_path, project_id)
result_file = joinpath(scratch_path, "pulse_exp_calib_$(calib_version)_borg.jlser");

import_rdata = false
if import_rdata

const mstags = ["L", "H", "M", "Sum"];
const conditions = ["Mock", "SC35M", "SC35MdelNS1"];

info("Loading $(joinpath(scratch_path, "$(project_id)_$(job_version)")).RData...")
pulse_rdata = load(joinpath(scratch_path, "$(project_id)_$(job_version).RData"), convert=true)
mschannels_df = pulse_rdata["ms_data"]["mschannels"];
mschannels_df[:timepoint] = categorical([float(get(x)) for x in mschannels_df[:timepoint]], true, ordered=true);
mschannels_df[:replicate] = categorical(mschannels_df[:replicate], true, ordered=true);
mschannels_df[:mstag] = levels!(categorical(mschannels_df[:mstag], true, ordered=true), mstags);
mschannels_df[:condition] = levels!(categorical(mschannels_df[:condition], true, ordered=true), conditions);
msrun_shifts_df = pulse_rdata["total_msrun_shifts.df"];
msrun_shifts_df[:msrun] = levels!(categorical(msrun_shifts_df[:msrun], true), levels(mschannels_df[:msrun]))

protein_concentration_df = pulse_rdata["protein_concentration.df"];
protein_concentration_df[:condition] = levels!(categorical(protein_concentration_df[:condition], true, ordered=true),
                                               levels(mschannels_df[:condition]));
protein_concentration_df[:timepoint] = categorical([parse(Float64, get(x)) for x in protein_concentration_df[:timepoint]], true, ordered=true);
protein_concentration_df[:replicate] = categorical([parse(Int, get(x)) for x in protein_concentration_df[:replicate]], true, ordered=true);
protgroups_df = pulse_rdata["ms_data"]["protgroups"];
protgroup_intensities_df = pulse_rdata["ms_data"]["protgroup_intensities"][:, [:mschannel, :protgroup_id, :Intensity, :Intensity_msrun_norm]];
rename!(protgroup_intensities_df, Dict(:Intensity_msrun_norm => :Intensity_norm));
protgroup_intensities_df[:mschannel] = levels!(categorical(protgroup_intensities_df[:mschannel], true),
                                               levels(mschannels_df[:mschannel]));
protgroup_intensities_df = join(protgroup_intensities_df,
                                mschannels_df[:, [:mschannel, :msrun, :mstag, :replicate, :condition, :timepoint, :is_msrun_used]],
                                on=:mschannel, kind=:inner);

open(GzipCompressorStream, joinpath(scratch_path, "$(project_id)_$(job_version).jlser"), "w") do io
    serialize(io, (protgroups_df, protgroup_intensities_df, protein_concentration_df, mschannels_df, msrun_shifts_df))
end;
else

protgroups_df, protgroup_intensities_df, protein_concentration_df, mschannels_df, msrun_shifts_df =
    open(deserialize, GzipDecompressorStream, joinpath(scratch_path, "$(project_id)_$(job_version).jlser"), "r");

end

sort!(protein_concentration_df, cols=[:condition, :timepoint, :replicate]);
concentration_arr = reshape(protein_concentration_df[:concentration],
                            (length(unique(protein_concentration_df[:replicate])),
                             length(unique(protein_concentration_df[:timepoint])),
                             length(unique(protein_concentration_df[:condition]))));

max_timepoint = maximum(levels(mschannels_df[:timepoint]))

@everywhere include(joinpath(base_scripts_path,"adhoc/phubel_pulse/pulse_dynamics/pulse_dynamics.jl"))
ϵ_scale = 1.0
ini_popmatrix = BlackBoxOptim.PopulationMatrix(0,0);
for i in 1:30
    info("Optimization cycle #$i...")
begin
    info("Running Borg...")
    flush(STDERR)
    flush(STDOUT)
    local t_rel = levels(mschannels_df[:timepoint])/max_timepoint
    global fit_problem = PulseSILACDynamics.CurveFitProblem(
        #PulseSILACDynamics.PolyharmonicSplineFactory(2, collect(linspace(0.0, 10.0, 4))),
        PulseSILACDynamics.PolyCurveFactory(3),
        t_rel, concentration_arr)
    fit_bb = PulseSILACDynamics.bbsetup(fit_problem, MaxTime=300.0, MaxSteps=10^8,
                                       ϵ=[10*ϵ_scale,ϵ_scale,10*ϵ_scale], TraceMode=:verbose, TraceInterval=10.0
                                       , Population=ini_popmatrix
                                       )
    global growth_fit_res = BlackBoxOptim.bboptimize(fit_bb)
    global growth_params = best_candidate(growth_fit_res)
    global growth_fitness = best_fitness(growth_fit_res)
    global ini_popmatrix = PulseSILACDynamics.popmatrix(growth_fit_res)

    if startswith(BlackBoxOptim.stop_reason(growth_fit_res), "No epsilon-progress") && length(pareto_frontier(growth_fit_res)) <= 500
        ϵ_scale *= 0.75
        info("New shrinked ϵ_scale=$ϵ_scale")
        if ϵ_scale < 0.01
            info("ϵ_scale limit reached, optimization cycles interrupted")
            break
        end
    end
    flush(STDERR)
    flush(STDOUT)
end
end

exp_calib_data = PulseSILACDynamics.import_pulse_silac_calib_data(protgroups_df, protgroup_intensities_df, protein_concentration_df);
growth_df = vcat([begin
    t = levels(mschannels_df[:timepoint])/max_timepoint
    g_curve = PulseSILACDynamics.curve(growth_params, i, fit_problem)
    logs = growth_params[i]

    cond_df = DataFrame(t=linspace(0.0, max_timepoint, 100),
                        t_rel=linspace(0.0, 1.0, 100))
    cond_df[:condition] = exp_calib_data.conditions[i]
    cond_df[:g] = exp.(g_curve.(cond_df[:t_rel]))
    cond_df[:g_scaled] = exp.(g_curve.(cond_df[:t_rel])+logs)
    cond_df[:g_rate] = map(t -> PulseSILACDynamics.derivative(g_curve, t), cond_df[:t_rel])
    cond_df
end for i in 1:size(concentration_arr, 3)]...);

using PlotlyJS

plot([PlotlyJS.AbstractTrace[scatter(x=trace_df[:t], y=trace_df[:g_scaled],
     name=trace_df[1, :condition],
     mode="lines+markers") for trace_df in DataFrames.groupby(growth_df, :condition)]],
    Layout(width=1600, height=1300, hovermode="closest"))

info("Saving JSON...")
open(joinpath(scratch_path, "protein_concentration_calib.json"), "w") do io
    output = Dict{Symbol,Any}(
        :raw_params => best_candidate(growth_fit_res),
        :fitness => best_fitness(growth_fit_res),
        :growth_dynamics => growth_df
    )
    JSON.print(io, output)
end

MethodTime = 600.0;

# try to load previous results
result_file = joinpath(scratch_path, "pulse_exp_calib_$(job_version)_borg.jlser");
if isfile(result_file)
    info("Loading existing experiments calibration from $result_file...")
    calib_params, exp_calib_data, exp_calib_fit_result_old, exp_calib_fit_old, cond_model_dyn_old, cond_growth_rate_old = open(deserialize, GzipDecompressorStream, result_file, "r");
    ini_popmatrix = Matrix{Float64}(length(BlackBoxOptim.params(first(exp_calib_fit_result_old))), length(exp_calib_fit_result_old));
    for (j, indi) in enumerate(exp_calib_fit_result_old)
        ini_popmatrix[:, j] = BlackBoxOptim.params(indi);
    end
else
    warn("No previous calibration results, preparing data...")
    #msruns, exp_design, sel_pepmods_data = gzopen(joinpath(analysis_path, "$(analysis_prefix)_exp_calib_data.jlser"), "r") do io deserialize(io) end;
    protgroups_df, protgroup_intensities_df, protein_concentration_df, mschannels_df, mschannel_shifts_df =
        open(deserialize, GzipDecompressorStream, joinpath(scratch_path, "$(project_id)_$(job_version).jlser"), "r");
    exp_calib_data = PulseSILACDynamics.import_pulse_silac_calib_data(protgroups_df, protgroup_intensities_df, protein_concentration_df);
    ini_popmatrix = BlackBoxOptim.PopulationMatrix(0,0);
end
flush(STDERR)
flush(STDOUT)

aux_problem = PulseSILACDynamics.bbsetup(exp_calib_data, MaxTime=MethodTime, MaxSteps=10^8, MaxStepsWithoutProgress=20000,
                       TraceMode=:verbose, TraceInterval=10.0,
                       RecoverResults=false, # we are saving intermediate results, so no need to recover
                       ϵ=[1.0, 1.0, 0.05]);
opt_archive = aux_problem.optimizer.evaluator.archive;

function save_results(calib_params, calib_model, calib_model_fitness,
                      problem_frontier, problem,
                      cycle::Int,
                      jlser_file::AbstractString, json_file::AbstractString)
    t_range =  calib_model.ode_solver.ode.t_range
    ts = collect(linspace(t_range[1], t_range[2], 101)[2:end])
    t_rels = collect(linspace(0.0, 1.0, 101)[2:end])

    cond_model_dyn = Vector{Matrix{Float64}}()
    cond_growth_rate = Vector{PulseSILACDynamics.ScaledSigmoid}()
    for i in 1:PulseSILACDynamics.nconditions(calib_model)
        model_dyn = PulseSILACDynamics.solve(calib_model, i, ts)
        push!(cond_model_dyn, model_dyn)
        gs = view(model_dyn, size(model_dyn, 1), :)
        growth_fit_problem = PulseSILACDynamics.ScaledSigmoidFitProblem(t_rels, gs)
        growth_fit_res = BlackBoxOptim.bboptimize(growth_fit_problem, Method=:borg_moea, MaxTime=600.0, MaxSteps=10^8,
                                                  MaxStepsWithoutProgress=30000, ϵ=0.2)
        push!(cond_growth_rate, PulseSILACDynamics.ScaledSigmoid(best_candidate(growth_fit_res)))
    end

    info("Saving optimization results into $jlser_file...")
    open(GzipCompressorStream, jlser_file, "w") do io
        serialize(io, (calib_params, exp_calib_data, problem_frontier,
                       calib_model, calib_model_fitness,
                       cond_model_dyn, cond_growth_rate))
    end

    info("Saving JSON $json_file...")
    open(json_file, "w") do io
        output = Dict{Symbol,Any}(
            :raw_params => calib_params,
            :fitness => calib_model_fitness,
            :exp_calib => calib_model,
            :protgroup_ids => exp_calib_data.protgroup_ids,
            :timepoints => exp_calib_data.timepoints,
            :conditions => exp_calib_data.conditions,
            :pulse_dynamics_per_condition => cond_model_dyn,
            :growth_rate_per_condition => cond_growth_rate,
            :growth_total_per_condition => [exp(PulseSILACDynamics.integrate(grate, 1.0) - PulseSILACDynamics.integrate(grate, 0.0)) for grate in cond_growth_rate],
        )
        JSON.print(io, output)
    end
end

res_jlser_file = joinpath(scratch_path, "pulse_exp_calib_$(calib_version)_borg.jlser");
res_json_file = joinpath(scratch_path, "pulse_exp_calib_$(calib_version)_borg.json");

info("Starting parallel optimization with $(nworkers()) workers...")
best_params = nothing
best_model = nothing
best_model_fitness = nothing
best_model_fitness_agg = NaN
borg_res = nothing # make it global
ϵ_scale = 4.0
last_opt_progress = 0
for i in 1:50
    info("Optimization cycle #$i: ϵ-scale=$(ϵ_scale)...")
    info("Running Borg...")
    flush(STDERR)
    flush(STDOUT)

    # fill the initial population with the current Pareto frontier
    ini_pop = Matrix{Float64}(numdims(aux_problem.problem), length(pareto_frontier(opt_archive)));
    for (j, indi) in enumerate(pareto_frontier(opt_archive))
        ini_pop[:, j] = BlackBoxOptim.params(indi);
    end
    #@everywhere include(joinpath(base_scripts_path,"adhoc/phubel_pulse/pulse_dynamics/pulse_dynamics.jl"))
    #exp_calib_data = PulseSILACDynamics.import_pulse_silac_calib_data(proteins_df, pulse_data_df);
    borg_opt = PulseSILACDynamics.bbsetup(exp_calib_data, MaxTime=MethodTime, MaxSteps=10^8, MaxStepsWithoutProgress=20000,
                           TraceMode=:verbose, TraceInterval=10.0,
                           RecoverResults=false, # we are saving intermediate results, so no need to recover
                           Workers=workers(),
                           ϵ=[0.5*ϵ_scale, ϵ_scale, 0.04*ϵ_scale],
                           Population=ini_pop, PopulationSize=200);
    borg_res = bboptimize(borg_opt);
    info("Done Borg...")
    # update the frontier
    for indi in pareto_frontier(borg_res)
        BlackBoxOptim.add_candidate!(opt_archive, BlackBoxOptim.fitness(BlackBoxOptim.params(indi), aux_problem.problem),
                                     BlackBoxOptim.params(indi), 0, i);
    end

    # don't trust the best_candidate() and check if it's really better than previous iteration;
    # due to discrete nature of ϵ-box frontier the new solution could be worse
    cur_fitness = best_fitness(borg_res)
    cur_fitness_agg = PulseSILACDynamics.posterior_probability(cur_fitness)
    nfrontier = length(pareto_frontier(borg_res))
    #opt_progress = !startswith(BlackBoxOptim.stop_reason(borg_res), "No epsilon-progress")
    if !isfinite(best_model_fitness_agg) || (best_model_fitness_agg < cur_fitness_agg)
        best_model_fitness = cur_fitness
        best_model_fitness_agg = cur_fitness_agg
        best_params = copy(best_candidate(borg_res))
        best_model = PulseSILACDynamics.params2model!(PulseSILACDynamics.empty_model(problem(borg_opt).factory),
                                                      problem(borg_opt).factory, best_params);
        println()
        show(best_model)
        println()
        BlackBoxOptim.show_fitness(STDOUT, best_model_fitness, problem(borg_opt))
        info("\nNew best fitness: ", best_model_fitness)
        opt_progress = true
    else
        opt_progress = false
    end

    info("Saving optimization results...")
    save_results(best_params, best_model, best_model_fitness, pareto_frontier(opt_archive),
                 problem(borg_opt), i, res_jlser_file, res_json_file)

    if !opt_progress && ϵ_scale < 0.1
        info("ϵ_scale limit reached, optimization cycles interrupted")
        break
    elseif !opt_progress && (i > last_opt_progress + 5)
        info("No model improvement for $(i - last_opt_progress) cycle(s), interrupted")
        break
    elseif (!opt_progress && nfrontier <= 200) || (nfrontier < 50)
        ϵ_scale *= 0.75
        info("New shrinked ϵ_scale=$ϵ_scale")
    elseif nfrontier > 1000
        ϵ_scale *= 2.0
    end
    opt_progress && (last_opt_progress = i)

    flush(STDERR)
    flush(STDOUT)
end

info("Done")

pf = pareto_frontier(borg_opt)
weighedfitness(fs) = fs[1]+fs[2]+1000.0*fs[3]
best_wfitness, best_idx = findmax(map(elm -> weighedfitness(fitness(elm)), pf))
opt_model = PulseSILACDynamics.params2model!(PulseSILACDynamics.empty_model(borg_opt.problem.factory), borg_opt.problem.factory, params(pf[best_idx]));

borg_opt = PulseSILACDynamics.bbsetup(exp_calib_data, MaxTime=5.0, MaxSteps=10^8, MaxStepsWithoutProgress=20000,
                       TraceMode=:verbose, TraceInterval=10.0,
                       RecoverResults=false, # we are saving intermediate results, so no need to recover
                       #Workers=workers(),
                       ϵ=[4*ϵ_scale,ϵ_scale,ϵ_scale], Population=ini_popmatrix);
borg_res = bboptimize(borg_opt);
opt_model = PulseSILACDynamics.params2model!(PulseSILACDynamics.empty_model(borg_opt.problem.factory), borg_opt.problem.factory, best_candidate(borg_res));

g_curves_df = vcat([begin
    t_scale = opt_model.ode_solver.ode.t_range[2] - opt_model.ode_solver.ode.t_range[1]
    cond_df = DataFrame(t=linspace(opt_model.ode_solver.ode.t_range[1], opt_model.ode_solver.ode.t_range[2], 100)[2:end],
                        t_rel=linspace(0.0, 1.0, 100)[2:end])
    cond_df[:condition] = borg_opt.problem.data.conditions[i]
    dyn = PulseSILACDynamics.solve(opt_model, i, cond_df[:t])
    cond_df[:g] = dyn[PulseSILACDynamics.nvariables(opt_model.ode_solver), :]
    cond_df[:g_rate] = PulseSILACDynamics.growth_rate(opt_model, i, cond_df[:t], dyn)
    cond_df[:s_rate] = map(opt_model.s_cond[i], cond_df[:t_rel])
    cond_df[:s] = t_scale.*map(PulseSILACDynamics.ScaledSigmoidIntegral(opt_model.s_cond[i]), cond_df[:t_rel])
    cond_df[:d_rate] = map(opt_model.d_cond[i], cond_df[:t_rel])
    cond_df[:d] = t_scale.*map(PulseSILACDynamics.ScaledSigmoidIntegral(opt_model.d_cond[i]), cond_df[:t_rel])
    cond_df
end for i in 1:PulseSILACDynamics.nconditions(opt_model)]...);

g_prot_curves_df = vcat([begin
    t_scale = opt_model.ode_solver.ode.t_range[2] - opt_model.ode_solver.ode.t_range[1]
    timepoints = linspace(opt_model.ode_solver.ode.t_range[1], opt_model.ode_solver.ode.t_range[2], 100)[2:end]
    timepoints_rel=linspace(0.0, 1.0, 100)[2:end]
    dyn = PulseSILACDynamics.solve(opt_model, i, timepoints)
    vcat([DataFrame(t=timepoints, t_rel=timepoints_rel, condition=borg_opt.problem.data.conditions[i],
                    protgroup_id=borg_opt.problem.data.protgroup_ids[j],
                    d_rate=map(opt_model.d[j,i], timepoints_rel),
                    d=t_scale.*map(PulseSILACDynamics.ScaledSigmoidIntegral(opt_model.d[j,i]), timepoints_rel),
                    s_rate=map(opt_model.s[j,i], timepoints_rel),
                    s=t_scale.*map(PulseSILACDynamics.ScaledSigmoidIntegral(opt_model.s[j,i]), timepoints_rel),
                    l=squeeze(dyn[(j-1)*3+1, :], 1),
                    h=squeeze(dyn[(j-1)*3+2, :], 1),
                    m=squeeze(dyn[(j-1)*3+3, :], 1),
                    sum=squeeze(sum(sub(dyn, (j-1)*3+(1:3), :), 1), 1)) for j in 1:PulseSILACDynamics.nproteins(opt_model)]...)
end for i in 1:PulseSILACDynamics.nconditions(opt_model)]...);

using PlotlyJS

plot([PlotlyJS.AbstractTrace[scatter(x=trace_df[:t], y=trace_df[:g_rate],
     name=trace_df[1, :condition],
     mode="lines+markers") for trace_df in DataFrames.groupby(g_curves_df, :condition)]],
    Layout(width=1600, height=1300, hovermode="closest"))

plot([PlotlyJS.AbstractTrace[scatter(x=trace_df[:t], y=trace_df[:g]./trace_df[28,:g],
     name=trace_df[1, :condition],
     mode="lines+markers") for trace_df in DataFrames.groupby(g_curves_df, :condition)]],
    Layout(width=1600, height=1300, hovermode="closest"))

plot([PlotlyJS.AbstractTrace[scatter(x=trace_df[:t], y=trace_df[:sum],
     name=string(trace_df[1, :protgroup_id], " ", trace_df[1, :condition]),
     mode="lines+markers") for trace_df in DataFrames.groupby(g_prot_curves_df, [:protgroup_id, :condition])],
     PlotlyJS.AbstractTrace[x for x in PlotlyJS.AbstractTrace[scatter(x=Float64[12.0,18.0,24.0,30.0,36.0,42.0],
             y=squeeze(sum(exp_calib_data.intensities_raw[(j-1)*3+(1:3),k,:,i],1),(1,2))/median(exp_calib_data.intensities_raw),
          name=string(exp_calib_data.protgroup_ids[j], " ", exp_calib_data.conditions[i], " #", k, " data"),
          mode="markers") for i in 1:PulseSILACDynamics.nconditions(exp_calib_data), j in 1:PulseSILACDynamics.nproteins(exp_calib_data), k in 1:PulseSILACDynamics.nreplicates(exp_calib_data)]]
    ],
    Layout(width=1600, height=1300, hovermode="closest"))

plot([PlotlyJS.AbstractTrace[scatter(x=trace_df[:t], y=trace_df[:d_rate],
     name=string(trace_df[1, :protgroup_id], " ", trace_df[1, :condition]),
     mode="lines+markers") for trace_df in DataFrames.groupby(g_prot_curves_df, [:protgroup_id, :condition])]],
    Layout(width=1600, height=1300, hovermode="closest"))

function condition_curves(channel_ix::Int)
    plot([PlotlyJS.AbstractTrace[scatter(
        x=1:6, y=squeeze(borg_opt.problem.data.intensities_raw[channel_ix,1,:,i,1], (1,2)),
        name=i,
        mode="lines+markers") for i in 1:3]])
end

@everywhere include(joinpath(base_scripts_path,"adhoc/phubel_pulse/pulse_dynamics/pulse_dynamics.jl"))
exp_calib_data = PulseSILACDynamics.import_pulse_silac_calib_data(proteins_df, pulse_data_df);
borg_opt = PulseSILACDynamics.bbsetup(exp_calib_data, MaxTime=3.0, MaxSteps=10^8, MaxStepsWithoutProgress=20000,
                       TraceMode=:verbose, TraceInterval=10.0,
                       RecoverResults=false, # we are saving intermediate results, so no need to recover
                       #Workers=workers(),
                       ϵ=[4*ϵ_scale,ϵ_scale,40*ϵ_scale], Population=ini_popmatrix);
borg_res = bboptimize(borg_opt);
save_results(BlackBoxOptim.problem(borg_opt), borg_res, res_jlser_file, res_json_file)
opt_model = PulseSILACDynamics.params2model!(PulseSILACDynamics.empty_model(borg_opt.problem.factory), borg_opt.problem.factory, best_candidate(borg_res));
growth_rates = PulseSILACDynamics.ScaledSigmoid[begin
    fit_problem = PulseSILACDynamics.ScaledSigmoidFitProblem(convert(Vector, cond_df[:t_rel]), convert(Vector, cond_df[:g]))
    fit_res = BlackBoxOptim.bboptimize(fit_problem, Method=:borg_moea, MaxTime=600.0, MaxSteps=10^8)
    PulseSILACDynamics.ScaledSigmoid(best_candidate(fit_res))
end for cond_df in DataFrames.groupby(g_curves_df, :condition)]
g_curves_df[:g_rate_fit] = 0.0
g_curves_df[:g_fit] = 0.0
for (i, gr) in enumerate(growth_rates)
    g = PulseSILACDynamics.ScaledSigmoidIntegral(gr)
    cond_mask = g_curves_df[:condition] .== exp_calib_data.conditions[i]
    g_curves_df[cond_mask, :g_fit] = exp(map(g, convert(Vector, g_curves_df[cond_mask, :t_rel])))
    g_curves_df[cond_mask, :g_rate_fit] = map(gr, convert(Vector, g_curves_df[cond_mask, :t_rel]))/42.0
end

plot([PlotlyJS.AbstractTrace[scatter(x=trace_df[:t], y=trace_df[:g]./trace_df[28,:g],
     name=trace_df[1, :condition],
     mode="lines+markers") for trace_df in DataFrames.groupby(g_curves_df, :condition)],
     PlotlyJS.AbstractTrace[scatter(x=trace_df[:t], y=trace_df[:g_fit]./trace_df[28,:g_fit],
          name=string(trace_df[1, :condition], " sim"),
          mode="lines+markers", markers=Dict(:line=>Dict(:dash=>"dashed"))) for trace_df in DataFrames.groupby(g_curves_df, :condition)]],
    Layout(width=1600, height=1300, hovermode="closest"))

plot([PlotlyJS.AbstractTrace[scatter(x=trace_df[:t], y=trace_df[:g_rate],
     name=trace_df[1, :condition],
     mode="lines+markers") for trace_df in DataFrames.groupby(g_curves_df, :condition)],
     PlotlyJS.AbstractTrace[scatter(x=trace_df[:t], y=trace_df[:g_rate_fit],
          name=string(trace_df[1, :condition], " sim"),
          mode="lines+markers", markers=Dict(:line=>Dict(:dash=>"dashed"))) for trace_df in DataFrames.groupby(g_curves_df, :condition)]],
    Layout(width=1600, height=1300, hovermode="closest"))

p = [condition_curves(1), condition_curves(2), condition_curves(3)];
p

plot([PlotlyJS.AbstractTrace[scatter(x=1:6, y=squeeze(borg_opt.problem.data.intensities_raw[1,i,:,1,1], (1,2)),
     name=i,
     mode="lines+markers") for i in 1:3]],
    Layout(width=1600, height=1300, hovermode="closest")
)
