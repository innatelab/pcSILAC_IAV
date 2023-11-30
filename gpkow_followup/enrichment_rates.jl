proj_info = (id = "yhuang_iav",
            modelobj = "protgroup",
            data_ver = "20171025",
            oesc_ver = "20230810"
)

const base_scripts_path = "/home/ge54heq/projects"
const base_analysis_path = "/pool/analysis/yhuang"

using Pkg
Pkg.activate(@__DIR__)

using Revise
using RData, DataFrames
using StatsBase
using DelimitedFiles

@info "Project '$(proj_info.id)' dataset version=$(proj_info.data_ver)"

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl");
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const results_path = joinpath(analysis_path, "results")
const reports_path = joinpath(analysis_path, "reports")
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")
const party3rd_data_path = "/pool/analysis/yhuang/pub3rdparty" #use my own pub3rdparty

includet(joinpath(misc_scripts_path, "frame_utils.jl"));
includet(joinpath(misc_scripts_path, "msglm_utils.jl"));
includet(joinpath(misc_scripts_path, "delimdata_utils.jl"));

objid_col = Symbol(string(proj_info.modelobj, "_id"));

full_rdata = load(joinpath(scratch_path, "phubel_pulse_$(proj_info.data_ver).RData"), convert=true)
data, header = readdlm(joinpath(data_path, "data.hits.intersected.relaxed.txt"), header = true)
contrasts_df = DataFrame(data, vec(header))|> MSGLMUtils.fix_object_id! 
contrasts_df.protgroup_id = Int64.(contrasts_df.protgroup_id);
contrasts_df[!,[:p_value, :median]] = Float64.(contrasts_df[!,[:p_value, :median]]);
contrasts_df[!, [:gene_name, :contrasts, :rate]] = string.(contrasts_df[!,[:gene_name, :contrasts, :rate]]);
contrasts_df.is_hit = contrasts_df[!, "is.hit"] .== "TRUE";

data, header = readdlm(joinpath(data_path, "phubel_pulse_clusters_20190131.txt"), header = true)
clusters_df = DataFrame(data, vec(header), makeunique = true)|> MSGLMUtils.fix_object_id! 
clusters_df.cluster= Int64.(clusters_df.cluster);

obj2protac_df = select!(
         full_rdata["ms_data"][string("protein2", proj_info.modelobj,"s")],
         [objid_col, :majority_protein_ac]) |> unique! |> MSGLMUtils.fix_object_id! 
rename!(obj2protac_df, :majority_protein_ac => :protein_ac)
objects_df = copy(full_rdata["ms_data"][string(proj_info.modelobj, "s")]) |> MSGLMUtils.fix_object_id!;

includet(joinpath(misc_scripts_path, "gmt_reader.jl"));
#add https://github.com/alyst/OptEnrichedSetCover.jl
include(joinpath(misc_scripts_path, "optcover_utils.jl"));
includet(joinpath(misc_scripts_path, "omics_collections.jl"));

@info "Loading Human annotations..."
# human mappings from http://download.baderlab.org/EM_Genesets/October_01_2020/Human/UniProt/
# using this older version to match the database used for the paper
# FIXME using all evidence codes
genesets_df, genesets_coll = GMT.read(String,
        #joinpath(party3rd_data_path, "Human_GO_AllPathways_with_GO_iea_October_01_2020_UniProt.gmt"),
        joinpath(party3rd_data_path, "Human_GO_AllPathways_with_GO_iea_August_08_2023_UniProt.gmt"), #use the latest one!
        id_col = :term_id, src_col = :term_src);
# strip Reactome version
genesets_df.term_id = [ifelse(r.term_src == "Reactome", replace(r.term_id, r"\.\d+$" => ""), r.term_id)
                       for r in eachrow(genesets_df)]
genesets_coll = Dict(ifelse(contains(k, r"^R-HSA-"), replace(k, r"\.\d+$" => ""), k) => v
                     for (k, v) in pairs(genesets_coll))

pcomplexes_df, pcomplex_iactors_df, pcomplex_iactor2ac_df =
    OmicsCollections.ppicollection(joinpath("/pool/pub3rdparty", "complexes_20191217.RData"), seqdb=:uniprot);
pcomplexes_df[!, :coll_id] .= "protein_complexes";
# make complexes collections, keep complexes with at least 2 participants
uprot_pcomplex_coll = FrameUtils.frame2collection(innerjoin(pcomplex_iactors_df, pcomplex_iactor2ac_df,
            on=[:file, :entry_index, :interaction_id, :interactor_id]),
            set_col=:complex_id, obj_col=:protein_ac, min_size=2)

protac_sets = merge!(genesets_coll, uprot_pcomplex_coll)

terms_df = vcat(rename(genesets_df[!, [:term_src, :term_id, :name, :descr]],
                       :term_src => :coll_id, :name=>:term_name, :descr=>:term_descr),
                rename(pcomplexes_df[!, [:coll_id, :complex_id, :interaction_label, :interaction_name]],
                       :complex_id=>:term_id, :interaction_label=>:term_name, :interaction_name=>:term_descr));
protac2term_df = FrameUtils.collection2frame(protac_sets, terms_df,
                                             setid_col=:term_id, objid_col=:protein_ac)

# link protein group IDs to annots and create protgroup collections
obj2term_df = select!(innerjoin(obj2protac_df, protac2term_df, on = :protein_ac),
                      Not([:protein_ac])) |> unique!
protac_colls = FrameUtils.frame2collections(protac2term_df, obj_col=:protein_ac,
                                            set_col=:term_id, coll_col=:coll_id)
obj_colls = FrameUtils.frame2collections(obj2term_df, obj_col=objid_col,
                                         set_col=:term_id, coll_col=:coll_id)

@info "Preparing mosaics..."
observed_protacs = Set(obj2protac_df.protein_ac) # all annotation ACs observed in the data
obj_mosaics = OptCoverUtils.collections2mosaics(obj_colls, 
                                    protac_colls, observed_protacs,
                                      setXset_frac_extra_elms=0.05,
                                      verbose=true);
                                      
ENV["MKL_NUM_THREADS"] = 1
using OptEnrichedSetCover
cover_params = CoverParams(setXset_factor=0.4,
                           uncovered_factor=0.5, covered_factor=0.0) #to match what Alexey used in his version

# GO enrichment for all the vs_Mock comparisons
obj_contrast_hit_sets = OmicsCollections.effects_collection(
    #contrasts_df,
    filter(r -> r.contrasts in (("SC35M_vs_Mock", "SC35MdelNS1_vs_Mock")), contrasts_df),
    obj_col=:object_id, change_col=:median, sel_col=:is_hit,
    group_cols=[:rate], #don't skip the group_cols even though it's not always meaningful
    effect_col=:contrasts)    

obj_contrast_hit_mosaics = Dict(begin
    @info "Masking $mosaic_name dataset by contrast hits..."
    mosaic_name => OptCoverUtils.automask(mosaic, obj_contrast_hit_sets,
                                          max_sets=2000, min_nmasked=2, max_setsize=200)
    end for (mosaic_name, mosaic) in pairs(obj_mosaics));

obj_contrast_hit_mosaics_v = collect(pairs(obj_contrast_hit_mosaics))
obj_contrast_hit_covers_v = similar(obj_contrast_hit_mosaics_v, Pair)
Threads.@threads for i in eachindex(obj_contrast_hit_mosaics_v)
    mosaic_name, masked_mosaic = obj_contrast_hit_mosaics_v[i]
    @info "Covering $mosaic_name by contrast hits..."
    obj_contrast_hit_covers_v[i] =
        mosaic_name => collect(masked_mosaic, cover_params,
            CoverEnumerationParams(max_set_score=0.0, max_covers=1),
            MultiobjOptimizerParams(ϵ=[0.02, 0.2], MaxSteps=2_000_000, WeightDigits=2, NWorkers=1),
            true)
end
obj_contrast_hit_covers = Dict(k => v for (k, v) in obj_contrast_hit_covers_v)

# GO enrichment for all the d1_vs_d0 comparisons
obj_contrast_d_hit_sets = OmicsCollections.effects_collection(
    filter(r -> r.rate == "d1_vs_d0", contrasts_df),
    obj_col=:object_id, change_col=:median, sel_col=:is_hit,
    group_cols=[:rate], 
    effect_col=:contrasts)    

obj_contrast_d_hit_mosaics = Dict(begin
    @info "Masking $mosaic_name dataset by contrast hits..."
    mosaic_name => OptCoverUtils.automask(mosaic, obj_contrast_d_hit_sets,
                                          max_sets=2000, min_nmasked=2, max_setsize=200)
    end for (mosaic_name, mosaic) in pairs(obj_mosaics));

obj_contrast_d_hit_mosaics_v = collect(pairs(obj_contrast_d_hit_mosaics))
obj_contrast_d_hit_covers_v = similar(obj_contrast_d_hit_mosaics_v, Pair)
Threads.@threads for i in eachindex(obj_contrast_d_hit_mosaics_v)
    mosaic_name, masked_mosaic = obj_contrast_d_hit_mosaics_v[i]
    @info "Covering $mosaic_name by contrast hits..."
    obj_contrast_d_hit_covers_v[i] =
        mosaic_name => collect(masked_mosaic, cover_params,
            CoverEnumerationParams(max_set_score=0.0, max_covers=1),
            MultiobjOptimizerParams(ϵ=[0.02, 0.2], MaxSteps=2_000_000, WeightDigits=2, NWorkers=1),
            true)
end
obj_contrast_d_hit_covers = Dict(k => v for (k, v) in obj_contrast_d_hit_covers_v)

# GO enrichment for each cluster
obj_clusters = FrameUtils.frame2collection(clusters_df, set_col=:cluster, obj_col=:object_id);
obj_clusters_sets = [obj_clusters[cix] for cix in 1:length(obj_clusters)];

obj_cluster_mosaics = Dict(begin
    @info "Masking $mosaic_name by clusters..."
    mosaic_name => OptCoverUtils.automask(mosaic, obj_clusters_sets,
                                          max_sets=2000, min_nmasked=2, max_setsize=200, max_pvalue=1E-3)
    end for (mosaic_name, mosaic) in pairs(obj_mosaics));

obj_cluster_mosaics_v = collect(pairs(obj_cluster_mosaics))
obj_cluster_covers_v = similar(obj_cluster_mosaics_v, Pair)
Threads.@threads for i in eachindex(obj_cluster_mosaics_v)
    mosaic_name, masked_mosaic = obj_cluster_mosaics_v[i]
    @info "Covering $mosaic_name by clusters..."
    obj_cluster_covers_v[i] =
        mosaic_name => collect(masked_mosaic, cover_params,
            CoverEnumerationParams(max_set_score=0.0, max_covers=1),
            MultiobjOptimizerParams(ϵ=[0.02, 0.2], MaxSteps=2_000_000, WeightDigits=2, NWorkers=1),
            true)
end
obj_cluster_covers = Dict(k => v for (k, v) in obj_cluster_covers_v)

using JLD2

@info "Saving data and analysis results"
hit_covers_filename = joinpath(scratch_path, "$(proj_info.id)_contrasts_$(proj_info.oesc_ver)_covers.jld2")
@save(hit_covers_filename,
      proj_info, protac_colls, obj_colls, obj_mosaics,
      obj2term_df, terms_df,
      objects_df, #obj_contrasts_df,
      #contrasts_df,
      obj_contrast_hit_sets, obj_contrast_d_hit_sets, obj_clusters_sets, 
      obj_contrast_hit_mosaics, obj_contrast_d_hit_mosaics, obj_cluster_mosaics,
      cover_params, obj_contrast_hit_covers, obj_contrast_d_hit_covers, obj_cluster_covers)
if !@isdefined(obj_hit_covers)
using JLD2, CSV, DataFrames, OptEnrichedSetCover
@load(hit_covers_filename,
      proj_info, protac_colls, obj_colls, obj_mosaics,
      obj2term_df, terms_df,
      objects_df, #obj_contrasts_df,
      #contrasts_df,
      obj_contrast_hit_sets, obj_contrast_d_hit_sets, obj_clusters_sets, 
      obj_contrast_hit_mosaics, obj_contrast_d_hit_mosaics, obj_cluster_mosaics,
      cover_params, obj_contrast_hit_covers, obj_contrast_d_hit_covers, obj_cluster_covers)
end

include(joinpath(misc_scripts_path, "optcover_utils.jl"));
@info "Preparing protgroup->gene_name map..."
obj_id2name = Dict(r.object_id => r[Symbol("gene_names")]
                   for r in eachrow(objects_df))

obj_contrast_hit_covers_df = 
    OptCoverUtils.covers_report(
    obj_contrast_hit_covers, obj_contrast_hit_sets, 
    obj_colls, 
    obj_id2name, terms_df,
    cover_params = cover_params,
    experimentid_col=[:rate, :contrast], weightedset_col_prefix="hit")
obj_contrast_hit_covers_df.intersect_genes = [join(unique(vcat(split.(split(genes, ' '), ';')...)), ' ') for genes in obj_contrast_hit_covers_df.intersect_genes]

# don't remove the sets 
obj_contrast_hit_covers_signif_df = combine(groupby(obj_contrast_hit_covers_df, :term_collection)) do coll_df
    @info "Processing $(coll_df.term_collection[1])..."
    df = OptCoverUtils.filter_multicover(coll_df, set_cols=[:contrast],
                                         max_term_pvalue=1E-3, max_set_pvalue=nothing, min_set_overlap=nothing)
    return select!(df, Not(:term_collection))
end

obj_contrast_d_hit_covers_df = 
    OptCoverUtils.covers_report(
    obj_contrast_d_hit_covers, obj_contrast_d_hit_sets, 
    obj_colls, 
    obj_id2name, terms_df,
    cover_params = cover_params,
    experimentid_col=[:rate, :contrast], weightedset_col_prefix="hit")
obj_contrast_d_hit_covers_df.intersect_genes = [join(unique(vcat(split.(split(genes, ' '), ';')...)), ' ') for genes in obj_contrast_d_hit_covers_df.intersect_genes]

# don't remove the sets 
obj_contrast_d_hit_covers_signif_df = combine(groupby(obj_contrast_d_hit_covers_df, :term_collection)) do coll_df
    @info "Processing $(coll_df.term_collection[1])..."
    df = OptCoverUtils.filter_multicover(coll_df, set_cols=[:contrast],
                                         max_term_pvalue=1E-3, max_set_pvalue=nothing, min_set_overlap=nothing)
    return select!(df, Not(:term_collection))
end

obj_cluster_covers_df = OptCoverUtils.covers_report(
    obj_cluster_covers, obj_clusters, 
    obj_colls, 
    obj_id2name, terms_df,
    cover_params = cover_params,
    experimentid_col=[:cluster], weightedset_col_prefix="cluster");
obj_cluster_covers_df.intersect_genes = [join(unique(vcat(split.(split(genes, ' '), ';')...)), ' ') for genes in obj_cluster_covers_df.intersect_genes]

obj_cluster_covers_signif_df = combine(groupby(obj_cluster_covers_df, :term_collection)) do coll_df
    @info "Processing $(coll_df.term_collection[1])..."
    df = OptCoverUtils.filter_multicover(coll_df, set_cols=[:cluster],
                                         max_term_pvalue=1E-4, max_set_pvalue=1E-3, min_set_overlap=nothing)
    return select!(df, Not(:term_collection))
end

using CSV
CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_contrasts_hit_covers_$(proj_info.oesc_ver).txt"),
          obj_contrast_hit_covers_df[obj_contrast_hit_covers_df.nmasked .> 0, :],
          missingstring="", delim='\t');
CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_contrasts_hit_covers_signif_$(proj_info.oesc_ver).txt"),
          obj_contrast_hit_covers_signif_df[obj_contrast_hit_covers_signif_df.nmasked .> 0, :],
          missingstring="", delim='\t');

CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_d1_vs_d0_hit_covers_$(proj_info.oesc_ver).txt"),
          obj_contrast_d_hit_covers_df[obj_contrast_d_hit_covers_df.nmasked .> 0, :],
          missingstring="", delim='\t');
CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_d1_vs_d0_hit_covers_signif_$(proj_info.oesc_ver).txt"),
          obj_contrast_d_hit_covers_signif_df[obj_contrast_d_hit_covers_signif_df.nmasked .> 0, :],
          missingstring="", delim='\t');

CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_cluster_covers_$(proj_info.oesc_ver).txt"),
          obj_cluster_covers_df[obj_cluster_covers_df.nmasked .> 0, :],
          missingstring="", delim='\t');
CSV.write(joinpath(analysis_path, "reports", "$(proj_info.id)_cluster_covers_signif_$(proj_info.oesc_ver).txt"),
          obj_cluster_covers_signif_df[obj_cluster_covers_signif_df.nmasked .> 0, :],
          missingstring="", delim='\t');

Revise.includet(joinpath(misc_scripts_path, "frame_utils.jl"))
Revise.includet(joinpath(misc_scripts_path, "optcover_plots.jl"))
include(joinpath(misc_scripts_path, "optcover_heatmap.jl"))

oesc_plots_path = joinpath(plots_path, "oesc_contrasts_$(proj_info.oesc_ver)")
isdir(oesc_plots_path) || mkdir(oesc_plots_path)

using PlotlyJS, TextWrap#, ORCA

heatmap_layout_attrs = Dict(
    ("SigDB_C2", true) => Dict(:margin_l => 500),
    ("SigDB_C2", false) => Dict(:margin_l => 500),
    #("Reactome", true) => Dict(:margin_l => 500),
    #("Reactome", false) => Dict(:margin_l => 500),
    ("Reactome", true) => Dict(:margin_l => 600),
    ("Reactome", false) => Dict(:margin_l => 600),
    ("GO_CC", true) => Dict(:margin_l => 200),
    ("GO_CC", false) => Dict(:margin_l => 200),
)

treatment_info = Dict("SC35M" => (color="#ee3a39", label = "IAV"),
                      "SC35MdelNS1"  => (color="#f79837", label = "IAVΔNS1"),
                      "Mock" => (color="gray", label = "Mock"))

function stylize_contrast_multi(contrast::AbstractString, condition_info::AbstractDict{<:AbstractString})
    contrast_match = match(r"^(?<lhs>.+)_vs_(?<rhs>[^+-]+)(?<sign>[+-])?$", contrast)
    if (isnothing(contrast_match))
        @warn "No format match for contrast '$(contrast)'"
        return contrast
    else
        res = OptCoverHeatmap.stylize_treatment(contrast_match[:lhs], condition_info) * " <b>vs</b> " *
        OptCoverHeatmap.stylize_treatment(contrast_match[:rhs], condition_info);
        if (!isnothing(contrast_match[:sign]))
            res *= OptCoverHeatmap.stylize_change(contrast_match[:sign]);
        end    
        return res
    end
end

function process_contrast_axis(contrasts_df)
    contrast_labels = contrasts_df.rate .* stylize_contrast_multi.(contrasts_df.contrast, Ref(treatment_info))

    contrasts_df,
    contrast_labels,
    string.("contrast: ", contrast_labels)
end


#treatment_order = Dict(:MPXV => 1, :VACV_WR => 2, :CVA_152 => 3, :MVA_F => 4, :MVA_dE9L => 5)

for term_coll in unique(obj_contrast_hit_covers_df.term_collection), signif in (false, true)
    @info "Plotting $(signif ? "signif " : "")contrast heatmap for $term_coll..."
    layout_attrs = get(heatmap_layout_attrs, (term_coll, signif), Dict())
    df = filter(r -> r.term_collection == term_coll, signif ? obj_contrast_hit_covers_signif_df : obj_contrast_hit_covers_df)
    if nrow(df) == 0
        @warn "No term_collection=$term_coll rows"
        continue
    end
    coll_heatmap = OptCoverHeatmap.oesc_heatmap(df,
            elements_label="protein",
            experiment_axis_title = "contrast",
            experiment_cols = [:rate, :contrast],
            process_experiment_axis=process_contrast_axis, 
            process_term_axis=OptCoverHeatmap.process_term_axis,
            margin_l=get(layout_attrs, :margin_l, 400),
            margin_b=get(layout_attrs, :margin_b, 80),
            cell_width=40, cell_height=30,
            colorscale = "Hot", reversescale=false,
            plot_bgcolor="#FFF", gridcolor="#DDD",
            #experiment_order=contrasts -> begin
            #    contrasts.treatment_order = [treatment_order[Symbol(r.treatment_lhs)]
            #                                for r in eachrow(contrasts)]
            #    return sortperm(contrasts, [:change, :treatment_order, :timepoint_lhs])
            #end
            )
    (coll_heatmap === nothing) && continue
    for (k, v) in [#:width=>800, :height=>400,
                   :margin_r=>80,
                   :yaxis_tickfont_size=>12, :xaxis_tickangle=>45]
        coll_heatmap.plot.layout[k] = v
    end
    for outformat in ["svg", "pdf","html"]
        plot_fname = joinpath(oesc_plots_path,
            "$(proj_info.id)_$(proj_info.oesc_ver)_$(term_coll)_contrasts$(signif ? "_signif" : "")_heatmap.$(outformat)")
        try
                savefig(coll_heatmap, plot_fname, format=outformat, width=coll_heatmap.plot.layout[:width], height=coll_heatmap.plot.layout[:height]);
        catch e
            if e isa InterruptException
                rethrow(e)
            else
                @warn "$term_coll generation failed: $e"
            end
        end
    end
end


for term_coll in unique(obj_contrast_d_hit_covers_df.term_collection), signif in (false, true)
    @info "Plotting $(signif ? "signif " : "")contrast heatmap for $term_coll..."
    layout_attrs = get(heatmap_layout_attrs, (term_coll, signif), Dict())
    df = filter(r -> r.term_collection == term_coll, signif ? obj_contrast_d_hit_covers_signif_df : obj_contrast_d_hit_covers_df)
    if nrow(df) == 0
        @warn "No term_collection=$term_coll rows"
        continue
    end
    coll_heatmap = OptCoverHeatmap.oesc_heatmap(df,
            elements_label="protein",
            experiment_axis_title = "contrast",
            experiment_cols = [:contrast, :rate, :nhit],
            #process_experiment_axis=process_contrast_axis,
            process_term_axis=OptCoverHeatmap.process_term_axis,
            margin_l=get(layout_attrs, :margin_l, 400),
            margin_b=get(layout_attrs, :margin_b, 80),
            cell_width=40, cell_height=30,
            colorscale = "Hot", reversescale=false,
            plot_bgcolor="#FFF", gridcolor="#DDD",
            #experiment_order=contrasts -> begin
            #    contrasts.treatment_order = [treatment_order[Symbol(r.treatment_lhs)]
            #                                for r in eachrow(contrasts)]
            #    return sortperm(contrasts, [:change, :treatment_order, :timepoint_lhs])
            #end
            )
    (coll_heatmap === nothing) && continue
    for (k, v) in [#:width=>800, :height=>400,
                   :margin_r=>80,
                   :yaxis_tickfont_size=>12, :xaxis_tickangle=>45]
        coll_heatmap.plot.layout[k] = v
    end
    for outformat in ["svg", "pdf","html"]
        plot_fname = joinpath(oesc_plots_path,
            "$(proj_info.id)_$(proj_info.oesc_ver)_$(term_coll)_d1_vs_d0_contrasts$(signif ? "_signif" : "")_heatmap.$(outformat)")
        try
                savefig(coll_heatmap, plot_fname, format=outformat, width=coll_heatmap.plot.layout[:width], height=coll_heatmap.plot.layout[:height]);
        catch e
            if e isa InterruptException
                rethrow(e)
            else
                @warn "$term_coll generation failed: $e"
            end
        end
    end
end

for term_coll in unique(obj_cluster_covers_df.term_collection), signif in (false, true)
    @info "Plotting $(signif ? "signif " : "")cluster heatmap for $term_coll..."
    layout_attrs = get(heatmap_layout_attrs, (term_coll, signif), Dict())
    df = filter(r -> r.term_collection == term_coll, signif ? obj_cluster_covers_signif_df : obj_cluster_covers_df)
    if nrow(df) == 0
        @warn "No term_collection=$term_coll rows"
        continue
    end
    coll_heatmap = OptCoverHeatmap.oesc_heatmap(df,
            elements_label="protein",
            experiment_axis_title = "cluster",
            experiment_cols = [:cluster, :ncluster],
            process_experiment_axis=OptCoverHeatmap.process_cluster_axis,
            process_term_axis=OptCoverHeatmap.process_term_axis,
            margin_l=get(layout_attrs, :margin_l, 400),
            margin_b=get(layout_attrs, :margin_b, 80),
            cell_width=25, cell_height=25,
            colorscale = "Hot", reversescale=false,
            plot_bgcolor="#FFF", gridcolor="#DDD")
    (coll_heatmap === nothing) && continue
    for (k, v) in [#:width=>800, :height=>400,
                   :margin_r=>80,
                   :yaxis_tickfont_size=>12, :xaxis_tickangle=>45]
        coll_heatmap.plot.layout[k] = v
    end
    for outformat in [#="svg",=# "pdf", "html"]
        plot_fname = joinpath(oesc_plots_path,
            "$(proj_info.id)_$(proj_info.oesc_ver)_$(term_coll)_clusters$(signif ? "_signif" : "")_heatmap.$(outformat)")
        try
                savefig(coll_heatmap, plot_fname, format=outformat, width=coll_heatmap.plot.layout[:width], height=coll_heatmap.plot.layout[:height]);
        catch e
            if e isa InterruptException
                rethrow(e)
            else
                @warn "$term_coll generation failed: $e"
            end
        end
    end
end

#only for debugging

modelobjid_col = Symbol(proj_info.modelobj, "_id")
ModelobjType = Int32

function add_modelobjid!(df::AbstractDataFrame)
    df.modelobj_id = copy(df[!, modelobjid_col])
    return df
end

function deduplicate_ids!(df::AbstractDataFrame; cols=[:gene_names, :majority_protein_acs, :protein_acs])
    for col in cols
        if !(col in propertynames(df))
            @warn "Column $col not found"
            continue
        end
        df[!, col] = [ismissing(r[col]) ? missing : DelimDataUtils.rejoin_unique_substrings([r[col]], delim=";", joindelim=";")
                      for r in eachrow(df)]
    end
    return df
end

function fix_msglm_columns!(df::AbstractDataFrame)
    df |> MSGLMUtils.fix_quantile_columns! |> MSGLMUtils.delete_service_cols! |> add_modelobjid! |> deduplicate_ids!
    if hasproperty(df, :p_value)
        rename!(df, :p_value => :pvalue)
    end
    if hasproperty(df, Symbol("50.0%"))
        df.median = copy(df[!, "50.0%"])
    end
    return df
end


fit_rdata = load(joinpath("/pool/analysis/astukalov/phubel_pulse/data", "phubel_pulse_$(proj_info.data_ver)_fit_all.RData"), convert=true)
fit_status_df = copy(fit_rdata["fit_status.df"]) |> add_modelobjid! |> deduplicate_ids!;
modelobj_global_df = copy(fit_rdata["fit_stats"]["global"]) |> fix_msglm_columns!;
#modelobj_recyc_df = copy(fit_rdata["fit_stats"]["recyc"]) |> fix_msglm_columns!; # may not exist
modelobj_scales_df = select!(filter(r -> r.var == "L0_labu", modelobj_global_df), [:modelobj_id, :median]);
modelobj_scales_df.modelobj_scale = exp.(max.(-5.0, modelobj_scales_df.median));
select!(modelobj_scales_df, Not(:median))

modelobj_rate_metrics_df = copy(fit_rdata["fit_stats"]["rate_params"]) |> fix_msglm_columns!;
rename!(modelobj_rate_metrics_df, :param => :metric, :condition => :treatment);
modelobj_rate_metrics_df.metric = replace.(modelobj_rate_metrics_df.metric, Ref("cum" => "sum"));
modelobj_rate_metrics_df.var = replace.(modelobj_rate_metrics_df.var, Ref("cum" => "sum"));
FrameUtils.matchcategoricals!(modelobj_rate_metrics_df, mschannels_df, cols=[:treatment]);
modelobj_rate_sum_df = filter(r -> startswith(r.metric, "sum"), modelobj_rate_metrics_df)
# add s_sum_rel
modelobj_s_rate_sum_df = innerjoin(filter(r -> startswith(r.var, "s_sum"), modelobj_rate_sum_df),
                                   modelobj_scales_df, on=:modelobj_id);
modelobj_s_rate_sum_df.var .= modelobj_s_rate_sum_df.var .* "_rel";
modelobj_s_rate_sum_df.metric .= modelobj_s_rate_sum_df.metric .* "_rel";
for col in MSGLMUtils.DefStat_cols
    (col in propertynames(modelobj_s_rate_sum_df)) || continue
    modelobj_s_rate_sum_df[!, col] ./= modelobj_s_rate_sum_df.modelobj_scale
end
select!(modelobj_s_rate_sum_df, Not(:modelobj_scale));
modelobj_rate_sum_df = vcat(modelobj_rate_sum_df, modelobj_s_rate_sum_df)
modelobj_rate_sum_df.var_treatment = modelobj_rate_sum_df.var .* "_" .* String.(modelobj_rate_sum_df.treatment);
categorical!(modelobj_rate_sum_df, [:metric, :var, :var_treatment])
categorical!(modelobj_rate_metrics_df, [:metric, :var])


modelobj_rate_metric_contrast_df = copy(fit_rdata["fit_cond_contrast_stats"]["rate_params"]) |> fix_msglm_columns!;
rename!(modelobj_rate_metric_contrast_df, :param => :metric, :contrasts => :contrast);
modelobj_rate_metric_contrast_df.rate = categorical(replace.(modelobj_rate_metric_contrast_df.var, Ref(r"_(.+)$" => "")));
modelobj_rate_metric_contrast_df.metric = replace.(modelobj_rate_metric_contrast_df.metric, Ref("cum" => "sum"));
modelobj_rate_metric_contrast_df.var = replace.(modelobj_rate_metric_contrast_df.var, Ref("cum" => "sum"));
contrast_matches = match.(Ref(r"(.+)_vs_(.+)"), modelobj_rate_metric_contrast_df.contrast);
modelobj_rate_metric_contrast_df.treatment_lhs = FrameUtils.matchcategorical([m[1] for m in contrast_matches], mschannels_df.treatment)
modelobj_rate_metric_contrast_df.treatment_rhs = FrameUtils.matchcategorical([m[2] for m in contrast_matches], mschannels_df.treatment)
categorical!(modelobj_rate_metric_contrast_df, [:metric, :rate, :var, :contrast])
modelobj_rate_sum_contrast_df = filter(r -> startswith(string(r.metric), "sum"), modelobj_rate_metric_contrast_df)
modelobj_s_rate_sum_contrast_df = innerjoin(filter(r -> r.rate == "s" && startswith(string(r.metric), "sum"), modelobj_rate_sum_contrast_df),
                                        modelobj_scales_df, on=:modelobj_id);
modelobj_s_rate_sum_contrast_df.var .= string.(modelobj_s_rate_sum_contrast_df.var) .* "_rel";
modelobj_s_rate_sum_contrast_df.metric = string.(modelobj_s_rate_sum_contrast_df.metric) .* "_rel";
for col in MSGLMUtils.DefStat_cols
    (col in propertynames(modelobj_s_rate_sum_contrast_df)) || continue
    modelobj_s_rate_sum_contrast_df[!, col] ./= modelobj_s_rate_sum_contrast_df.modelobj_scale
end
select!(modelobj_s_rate_sum_contrast_df, Not(:modelobj_scale));

modelobj_rate_sum_contrast_df = vcat(modelobj_rate_sum_contrast_df, modelobj_s_rate_sum_contrast_df);
modelobj_rate_sum_contrast_df.rate_metric = String.(modelobj_rate_sum_contrast_df.rate) .* "_" .* String.(modelobj_rate_sum_contrast_df.metric);
modelobj_rate_sum_contrast_df.rate_metric_contrast = modelobj_rate_sum_contrast_df.rate_metric .* "_" .* String.(modelobj_rate_sum_contrast_df.contrast);
modelobj_rate_sum_contrast_df.mlog10_pvalue = -log10.(modelobj_rate_sum_contrast_df.pvalue);
modelobj_rate_sum_contrast_df.is_hit_strong = [(r.rate_metric == "d1_sum" || r.rate_metric == "d0_sum") ?
                                                (r.pvalue <= 1E-4) && (abs(r.median) >= 0.25) :
                                               (r.rate_metric == "s_sum_scaled_rel") ?
                                                (r.pvalue <= 1E-4) && (abs(r.median) >= 2) : missing
                                               for r in eachrow(modelobj_rate_sum_contrast_df)]
modelobj_rate_sum_contrast_df.is_hit_relaxed = [(r.rate_metric == "d1_sum" || r.rate_metric == "d0_sum") ?
                                                (r.pvalue <= 1E-3) && (abs(r.median) >= 0.25) :
                                               (r.rate_metric == "s_sum_scaled_rel") ?
                                                (r.pvalue <= 1E-3) && (abs(r.median) >= 1) : missing
                                               for r in eachrow(modelobj_rate_sum_contrast_df)]
categorical!(modelobj_rate_sum_contrast_df, [:var, :rate, :metric, :rate_metric, :rate_metric_contrast])
