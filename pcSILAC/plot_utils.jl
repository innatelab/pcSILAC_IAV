module PulsePlotUtils

using Statistics, DataFrames, PlotlyJS
using Printf: @sprintf

const ConditionsPalette = Dict("Mock" => "#808080", "SC35M" => "#FF0000", "SC35MdelNS1" => "#FFA500")
const RatesPalette = Dict("s" => "MediumSeaGreen", "d0" => "DimGray", "d1" => "LightSteelBlue")

function append_CI_band_traces!(applystyle!::Function, traces::AbstractArray{<:PlotlyJS.AbstractTrace},
                                name::AbstractString,
                                data::AbstractDataFrame, xcol::Symbol,
                                midline_col::Union{Symbol, Nothing} = Symbol("50.0%"),
                                band_cols::Tuple{Symbol,Symbol}=(Symbol("25.0%"), Symbol("75.0%")),
                                fill::Bool = true, hover::Bool=true)
    # CI band
    push!(traces, applystyle!(scatter(data, name=name, x=xcol, y=band_cols[1],
                                      line_width=0.25, mode=:lines, hoverinfo="none"), band_cols[1]))
    push!(traces, applystyle!(scatter(data, name=name, x=xcol, y=band_cols[2],
                                      line_width=0.25, mode=:lines, hoverinfo="none", fill=fill ? "tonexty" : "none"), band_cols[2]))
    # midline
    if midline_col !== nothing
        push!(traces, applystyle!(scatter(data, name=name, x=xcol, y=midline_col,
                                          line_width=2.0, mode=:lines, hoverinfo=hover ? "text" : "none"), midline_col))
    end
    return traces
end

struct SubplotsDef
    subplot_col::Symbol
    trace_col::Symbol
    signal_col::Symbol

    trace2color::Dict{String, String}
    showlegend::Bool
    hover::Bool

    SubplotsDef(; subplot_col::Symbol=:mstag, trace_col::Symbol=:condition,
                signal_col::Symbol=:signal,
                trace2color::Dict{String, String} = ConditionsPalette,
                showlegend::Bool=true, hover::Bool=true) =
        new(subplot_col, trace_col, signal_col,
            trace2color, showlegend, hover)
end

#  PlotlyJS plots for each mstags of Pulse/Chase dynamics timecourse for a given protgroup
# returns array of plots
function data_subplots(intens_df::AbstractDataFrame, def::SubplotsDef)
    qtls = [0.25, 0.5, 0.75]
    qtl_colnames = [Symbol(@sprintf("%.1f%%", qtl*100)) for qtl in qtls]
    has_inference = :is_inferred ∈ names(intens_df)
    res = Vector{PlotlyJS.Plot}()
    showlegend = def.showlegend
    for subplot_df in groupby(intens_df, def.subplot_col)
        subplot_name = subplot_df[1, def.subplot_col]
        subplot_traces = Vector{AbstractTrace}()
        for trace_df in groupby(subplot_df, def.trace_col)
            trace_name = trace_df[1, def.trace_col]
            # individual measurements
            cloud_tr = scatter(trace_df, legendgroup=trace_name, showlegend=showlegend, x=:timepoint, y=def.signal_col,
                               marker_color=def.trace2color[trace_name], marker_line_width=0.0,
                               name=trace_name, mode=:markers, hoverinfo=def.hover ? "text" : "none")
            if def.hover
                cloud_tr[:text]=["$(get(r.timepoint))h #$(get(r.replicate)): $(@sprintf("%.3g", r[def.signal_col]))" for r in eachrow(trace_df)]
            end
            if has_inference
                cloud_tr[:marker_size] = ifelse.(trace_df.is_inferred, 2.0, 5.0)
                if haskey(cloud_tr, :text)
                    for (i, infer) in enumerate(trace_df.is_inferred)
                        infer || continue
                        cloud_tr[:text][i] = cloud_tr[:text][i] * " (inferred)"
                    end
                end
            end
            push!(subplot_traces, cloud_tr)
            trace_stats_df = combine(groupby(trace_df, :timepoint)) do tp_df
                tp_qtls = quantile(tp_df[def.signal_col], qtls)
                DataFrame([[tp_qtl] for tp_qtl in tp_qtls], qtl_colnames)
            end
            append_CI_band_traces!(subplot_traces, trace_name, trace_stats_df,
                                   :timepoint, Symbol("50.0%"), (Symbol("25.0%"), Symbol("75.0%")), true, def.hover) do trace, col
                trace[:legendgroup] = trace_name
                trace[:showlegend] = false
                trace[:marker_color] = def.trace2color[trace_name]
                if trace[:hoverinfo] == "text"
                    trace[:text] = ["$(get(r.timepoint))h median: $(@sprintf("%.3g", r[col]))" for r in eachrow(trace_stats_df)]
                end
                return trace
            end
        end
        showlegend = false # don't duplicate trace legends for the other subplots
        push!(res, Plot(subplot_traces, Layout(hovermode=def.hover ? "closest" : "skip",
                        yaxis_title=subplot_name,
                        xaxis_tickvals=unique(intens_df.timepoint),
                        xaxis_tickmode="array",
                        xaxis_showticklabels=false)))
    end
    isempty(res) || (res[end].layout[:xaxis_showticklabels]=true)
    return res
end

function dynsim_subplots(sim_df::AbstractDataFrame,
                         intens_df::Union{AbstractDataFrame, Nothing},
                         def::SubplotsDef, t_observed::Number = 12.5)
    qtls = [0.025, 0.25, 0.5, 0.75, 0.975]
    qtl_colnames = [Symbol(@sprintf("%.1f%%", qtl*100)) for qtl in qtls]
    qtl_shown_colnames = [Symbol(@sprintf("%.1f%%_shown", qtl*100)) for qtl in qtls]
    signal_shown_col = Symbol(string(def.signal_col)*"_shown")
    res = Vector{PlotlyJS.Plot}()
    if intens_df !== nothing
        intens_df = intens_df[.!ismissing.(intens_df[def.signal_col]), :]
        intens_df[def.signal_col] = convert(Vector{Missings.T(eltype(intens_df[def.signal_col]))}, intens_df[def.signal_col])
        has_inference = :is_inferred ∈ names(intens_df)
    else
        has_inference = false
    end
    showlegend = def.showlegend
    for subplot_df in groupby(sim_df, def.subplot_col)
        subplot_name = subplot_df[1, def.subplot_col]
        # limit shown intensity in case of divergence
        timepoints = subplot_df.timepoint
        if eltype(timepoints) <: CategoricalValue
            timepoints = get.(timepoints)
        end
        observed_mask = timepoints .> t_observed
        max_sim_shown = min(quantile(subplot_df[observed_mask, qtl_colnames[4]], 0.99) * 1.25,
                            maximum(subplot_df[observed_mask, qtl_colnames[5]]))
        subplot_shown_df = convert(DataFrame, subplot_df[[:timepoint, def.trace_col]])
        for (qtl_col, qtl_shown_col) in zip(qtl_colnames, qtl_shown_colnames)
            subplot_shown_df[qtl_col] = subplot_df[qtl_col]
            subplot_shown_df[qtl_shown_col] = min.(subplot_df[qtl_col], max_sim_shown)
        end
        if intens_df !== nothing
            subplot_intens_df = intens_df[intens_df[def.subplot_col] .== subplot_name, :]
            signal_data = subplot_intens_df[def.signal_col]
            max_intens_shown = max(max_sim_shown, min(quantile(signal_data, 0.95)*1.25,
                                   maximum(signal_data), 2*max_sim_shown))
            subplot_intens_df[signal_shown_col] = min.(signal_data, max_intens_shown)
        else
            subplot_intens_df = nothing
        end
        subplot_traces = Vector{AbstractTrace}()
        for trace_df in groupby(subplot_shown_df, def.trace_col)
            trace_name = trace_df[1, def.trace_col]
            #@show subplot_name trace_name
            # individual measurements
            if subplot_intens_df !== nothing
                cloud_tr = scatter(subplot_intens_df[subplot_intens_df[def.trace_col].==trace_name, :],
                        legendgroup=trace_name, showlegend=false,
                        x=:timepoint, y=signal_shown_col,
                        marker_color=def.trace2color[trace_name], marker_line_width=0.0,
                        name=trace_name, mode=:markers, hoverinfo="skip") # don't show labels for measurements
                if def.hover
                    cloud_tr[:text]=["$(get(r.timepoint))h #$(get(r.replicate)): $(@sprintf("%.3g", r[def.signal_col]))" for r in eachrow(subplot_intens_df)]
                end
                if has_inference
                    cloud_tr[:marker_size] = ifelse.(subplot_intens_df.is_inferred .| (subplot_intens_df[def.signal_col] .!= subplot_intens_df[signal_shown_col]), 2.0, 5.0)
                    if def.hover
                        for (i, infer) in enumerate(subplot_intens_df.is_inferred)
                            infer || continue
                            cloud_tr[:text][i] = cloud_tr[:text][i] * " (inferred)"
                        end
                    end
                end
                push!(subplot_traces, cloud_tr)
            end
            append_CI_band_traces!(subplot_traces, trace_name, trace_df,
                                   :timepoint, Symbol("50.0%_shown"), (Symbol("25.0%_shown"), Symbol("75.0%_shown")), true, def.hover) do trace, col
                trace[:legendgroup] = trace_name
                trace[:showlegend] = showlegend && (col == Symbol("50.0%_shown"))
                trace[:marker_color] = def.trace2color[trace_name]
                if trace[:hoverinfo] == "text"
                    real_col = Symbol(replace(string(col), r"_shown$" => ""))
                    trace[:text] = ["$(r.timepoint)h median: $(@sprintf("%.3g", r[real_col]))" for r in eachrow(trace_df)]
                end
                return trace
            end
            append_CI_band_traces!(subplot_traces, trace_name, trace_df, :timepoint,
                                   nothing, (Symbol("2.5%_shown"), Symbol("97.5%_shown")), false, false) do trace, col
                trace[:legendgroup] = trace_name
                trace[:showlegend] = false
                trace[:marker_color] = def.trace2color[trace_name]
                trace[:line_dash] = "dot"
                return trace
            end
        end
        showlegend = false
        subp = Plot(subplot_traces, Layout(hovermode=def.hover ? "x" : "skip",
                    yaxis_title=subplot_name,
                    xaxis_showticklabels=false))
        if intens_df !== nothing
            subp.layout[:xaxis_tickvals]=vcat(minimum(levels(sim_df.timepoint)), unique(intens_df.timepoint))
            subp.layout[:xaxis_tickmode]="array"
        end
        push!(res, subp)
    end
    isempty(res) || (res[end].layout[:xaxis_showticklabels]=true)
    return res
end

function protgroup_label(df::AbstractDataFrame)
    res = "Protgroup #$(df.protgroup_id[1])"
    if :gene_names ∈ names(df)
        res = res * ": $(df.gene_names[1])"
    end
    if :majority_protein_acs ∈ names(df)
        res = res * " ($(df.majority_protein_acs[1]))"
    end
    return res
end

filter_data(df::Union{AbstractDataFrame, Nothing}, protgroup_id::Integer) =
    (df !== nothing) ? df[df.protgroup_id .== protgroup_id, :] : nothing

filter_data(df::Union{AbstractDataFrame, Nothing}, protgroup_id::Integer,
            subplot_col::Union{Symbol, Nothing}, ::Nothing) = filter_data(df, protgroup_id)

function filter_data(df::Union{AbstractDataFrame, Nothing}, protgroup_id::Integer,
                     subplot_col::Symbol,
                     subplots::Union{AbstractSet{T}, AbstractVector{T}}) where T
    pg_df = filter_data(df, protgroup_id)
    subplot_set = subplots isa AbstractSet ? subplots : Set(subplots)
    if isa(pg_df[subplot_col], AbstractCategoricalArray)
        return pg_df[get.(pg_df[subplot_col]) .∈ subplot_set, :]
    else
        return pg_df[pg_df[subplot_col] .∈ subplot_set, :]
    end
end

# PlotlyJS plot of Pulse/Chase dynamics timecourse for a given protgroup
function plot_pulse_data(intens_df::AbstractDataFrame, protgroup_id;
                         signal_col=:signal, mstags = nothing, cond2color=ConditionsPalette)
    sel_intens_df = filter_data(intens_df, protgroup_id, :mstag, mstags)

    mstag_plots = data_subplots(sel_intens_df, SubplotsDef(signal_col=signal_col, trace2color=cond2color))
    res = arrange(mstag_plots)
    res.layout[:hovermode] = "closest"
    res.layout[:title] = "Pulse/Chase data " * protgroup_label(intens_df, protgroup_id)
    return res
end

# PlotlyJS plots of Pulse/Chase dynamics timecourse for a given protgroup
function plot_pulse_data(intens_dfs::AbstractVector{<:AbstractDataFrame}, protgroup_id::Integer;
                         signal_col=:signal, mstags = nothing,
                         cond2color=ConditionsPalette, rate2color=RatesPalette)
    tag_plots = Vector{PlotlyJS.Plot}()
    ntags = 0
    nplots = 0
    for (i, intens_df) in enumerate(intens_dfs)
        sel_intens_df = filter_data(intens_df, protgroup_id, :mstag, mstags)
        df_tag_plots = data_subplots(sel_intens_df,
                SubplotsDef(signal_col=signal_col, showlegend=i==1, trace2color=cond2color))
        isempty(df_tag_plots) && continue # skip empty, FIXME add empty subplot
        if ntags == 0
            ntags = length(df_tag_plots)
        elseif ntags != length(df_tag_plots)
            error("Different number of tags")
        end
        append!(tag_plots, df_tag_plots)
        nplots += 1
    end
    res = arrange(reshape(tag_plots, (ntags, nplots)))
    res.layout[:hovermode] = "closest"
    res.layout[:title] = "Pulse/Chase data " * protgroup_label(filter_data(intens_dfs[1], protgroup_id))
    return res
end

function variable_rows_layout(col2nrow::AbstractVector{<:Integer})
    w, dx = PlotlyJS.subplot_size(1, length(col2nrow), false)[[1, 3]]
    out = Layout()
    x = 0.0
    subplot = 1
    for (col, nr) in enumerate(col2nrow)
        xdom = [x, x + w]::Vector{Float64}

        if nr > 0
            h, dy = PlotlyJS.subplot_size(nr, length(col2nrow), false)[[2, 4]]
            y = 1.0 # start from top
            for row in 1:nr
                out["xaxis$subplot"] = Dict([:domain=>copy(xdom), :anchor=>"y$subplot"])
                out["yaxis$subplot"] = Dict([:domain=>[y - h, y], :anchor=>"x$subplot"])
                #@show row col subplot sub_code x y w h
                subplot += 1
                y -= (dy + h)
            end
        end
        x += dx + w
    end
    out
end

function plot_pulse_dynsim(sim_mstag_df::AbstractDataFrame, protgroup_id::Integer,
                           intens_df::Union{AbstractDataFrame, Nothing} = nothing;
                           signal_col=:signal, mstags = nothing,
                           cond2color=ConditionsPalette)

    sel_sim_mstag_df = filter_data(sim_mstag_df, protgroup_id, :mstag, mstags)
    sel_intens_df = filter_data(intens_df, protgroup_id, :mstag, mstags)

    plots = dynsim_subplots(sel_sim_mstag_df, sel_intens_df,
                            SubplotsDef(signal_col=signal_col, trace2color=cond2color))
    if !isempty(plots)
        res = arrange(plots)
        res.layout[:hovermode] = "x"
        res.layout[:title] = "Pulse/Chase simulation " * protgroup_label(sel_sim_mstag_df)
    else
        res = nothing
    end
    return res
end

function plot_pulse_dynsim(sim_mstag_df::AbstractDataFrame, sim_rate_df::AbstractDataFrame,
                           protgroup_id::Integer, intens_df::Union{AbstractDataFrame, Nothing} = nothing;
                           signal_col=:signal_scaled, mstags = nothing,
                           cond2color=ConditionsPalette, rate2color=RatesPalette)
    sel_sim_mstag_df = filter_data(sim_mstag_df, protgroup_id, :mstag, mstags)
    sel_sim_rate_df = filter_data(sim_rate_df, protgroup_id)
    sel_intens_df = filter_data(intens_df, protgroup_id, :mstag, mstags)

    sim_mstag_plots = dynsim_subplots(sel_sim_mstag_df, sel_intens_df,
            SubplotsDef(signal_col=signal_col, trace2color=cond2color))
    sim_rate_plots = dynsim_subplots(sel_sim_rate_df, nothing,
            SubplotsDef(signal_col=signal_col, subplot_col=:rate, showlegend=false, trace2color=cond2color), 8.0)
    sel_sim_drate_df = sel_sim_rate_df[sel_sim_rate_df.rate .∈ Ref(Set(["d0", "d1"])), :]
    sim_drate_plots = dynsim_subplots(sel_sim_drate_df, nothing,
        SubplotsDef(signal_col=signal_col, subplot_col=:condition, trace_col=:rate, trace2color=rate2color), 8.0)

    res = arrange(variable_rows_layout(filter!(x -> x > 0, [length(sim_mstag_plots), length(sim_rate_plots), length(sim_drate_plots)])),
                  vcat(sim_mstag_plots, sim_rate_plots, sim_drate_plots))
    #res.layout[:hovermode] = "closest"
    res.layout[:title] = "Pulse/Chase simulation " * protgroup_label(sim_mstag_df)
    return res
end

function plot_pulse_ratedyn(sim_rate_df::AbstractDataFrame,
                            protgroup_id::Integer;
                            cond2color=ConditionsPalette, rate2color=RatesPalette)
    sel_sim_rate_df = filter_data(sim_rate_df, protgroup_id)

    plots = dynsim_subplots(sel_sim_rate_df, nothing,
                            SubplotsDef(signal_col=:not_available, subplot_col=:rate,
                                        showlegend=true, trace2color=cond2color), 8.0)
    if !isempty(plots)
        res = arrange(plots)
        res.layout[:hovermode] = "x"
        res.layout[:title] = "Pulse/Chase simulation " * protgroup_label(sel_sim_rate_df)
    else
        res = nothing
    end

    return res
end

function get_protgroup_id(protgroups_df::AbstractDataFrame;
                          protgroup_id::Union{Integer, Nothing} = nothing,
                          gene_name::Union{String, Nothing} = nothing)
    if protgroup_id !== nothing
        gene_name !== nothing && @warn("gene_name=$gene_name keyarg ignored")
        return eltype(protgroups_df.protgroup_id)(protgroup_id)
    end
    gn_rx = gene_name isa String ? Regex("(^|;)"*gene_name*"(;|\$)") : gene_name
    pg_ixs = findall(gn -> !ismissing(gn) && occursin(gn_rx, gn), protgroups_df.gene_names)
    if length(pg_ixs) == 0
        error("Gene name $gene_name not found")
    elseif length(pg_ixs) > 1
        error("Gene name $gene_name matches $(length(pg_ixs)) protein group(s): $(join(protgroups_df.protgroup_id[pg_ixs], ' '))")
    end
    pg_ix = pg_ixs[1]
    return protgroups_df.protgroup_id[pg_ix]
end

function plot_pulse_data(intens_df::AbstractDataFrame, protgroups_df;
                         protgroup_id::Union{Integer, Nothing} = nothing,
                         gene_name::Union{String, Nothing} = nothing,
                         kwargs...)
    pg_id = get_protgroup_id(protgroups_df, protgroup_id=protgroup_id, gene_name=gene_name)
    plot_pulse_data(intens_df, pg_id; kwargs...)
end

function plot_pulse_data(intens_df::AbstractVector{<:AbstractDataFrame}, protgroups_df;
                         protgroup_id::Union{Integer, Nothing} = nothing,
                         gene_name::Union{String, Nothing} = nothing,
                         kwargs...)
    pg_id = get_protgroup_id(protgroups_df, protgroup_id=protgroup_id, gene_name=gene_name)
    plot_pulse_data(intens_df, pg_id; kwargs...)
end

function plot_pulse_dynsim(sim_mstag_df::AbstractDataFrame,
                           protgroups_df::AbstractDataFrame,
                           intens_df::Union{AbstractDataFrame, Nothing} = nothing;
                           protgroup_id::Union{Integer, Nothing} = nothing,
                           gene_name::Union{String, Nothing} = nothing,
                           kwargs...)
    pg_id = get_protgroup_id(protgroups_df, protgroup_id=protgroup_id, gene_name=gene_name)
    res = plot_pulse_dynsim(sim_mstag_df, pg_id, intens_df;
                            kwargs...)
    if res !== nothing
        res.layout[:title] = "Pulse/Chase simulation " * protgroup_label(protgroups_df[protgroups_df.protgroup_id.==pg_id, :])
    end
    return res
end

function plot_pulse_dynsim(sim_mstag_df::AbstractDataFrame, sim_rate_df::AbstractDataFrame,
                           protgroups_df::AbstractDataFrame,
                           intens_df::Union{AbstractDataFrame, Nothing} = nothing;
                           protgroup_id::Union{Integer, Nothing} = nothing,
                           gene_name::Union{String, Nothing} = nothing,
                           kwargs...)
    pg_id = get_protgroup_id(protgroups_df, protgroup_id=protgroup_id, gene_name=gene_name)
    res = plot_pulse_dynsim(sim_mstag_df, sim_rate_df, pg_id, intens_df; kwargs...)
    res.layout[:title] = "Pulse/Chase simulation " * protgroup_label(protgroups_df[protgroups_df.protgroup_id.==pg_id, :])
    return res
end

function embed_hovertext(r)
    res =
    "$(r.gene_names)<br>
     $(r.protein_names)<br>
     $(r.majority_protein_acs)<br>
     $(r.protgroup_id)<br>"
    if :category in names(r)
        res = res * "category: $(r.category)"
        res *= :cluster in names(r) ? ", " : "<br>\n"
    end
    if :cluster in names(r)
        res *= "cluster #$(r.cluster)<br>"
    end
    res *= "N[quant]=$(r.nquant) log(L0)=$(@sprintf("%.2g", r.L0_labu))"
    if :recyc_L in names(r)
        res *= " recyc[L]=$(@sprintf("%.1f%%", r.recyc_L*100)) recyc[H]=$(@sprintf("%.1f%%", r.recyc_H*100))"
    end
    return res
end

end
