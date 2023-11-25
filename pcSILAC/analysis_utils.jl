module PulseAnalysisUtils

using DataFrames, Statistics

"""
Convert data frame into 3D tensor with `data_col` values and
return it along with the axes information (as a 4-tuple):
  * fractions (`fraction_col` values)
  * objects (`obj_cols` subframe)
  * experiments (`exp_cols` subframe).
"""
function expdata_tensor(data_df::AbstractDataFrame;
                        data_col=:Intensity_norm,
                        exp_cols=[:condition, :replicate],
                        tag_cols=[:mstag],
                        timepoint_cols=[:timepoint],
                        obj_cols=[:protgroup_id],
                        verbose::Bool=false)
    verbose && @info("Sorting the data...")
    sort!(data_df, vcat(obj_cols, exp_cols, tag_cols, timepoint_cols));
    verbose && @info("Reshaping into 4-tensor...")
    res = reshape(coalesce.(data_df[!, data_col], 0.0),
                  size(unique(data_df[!, timepoint_cols]), 1),
                  size(unique(data_df[!, tag_cols]), 1),
                  size(unique(data_df[!, exp_cols]), 1),
                  size(unique(data_df[!, obj_cols]), 1))
    verbose && @info("Returning the 3D-tensor and the axes labels...")
    res,
    # axes info
    data_df[1:size(res, 1), timepoint_cols],
    data_df[(0:(size(res, 2)-1))*stride(res, 2).+1, tag_cols],
    data_df[(0:(size(res, 3)-1))*stride(res, 3).+1, exp_cols],
    data_df[(0:(size(res, 4)-1))*stride(res, 4).+1, obj_cols]
end

# add sum signal for each msrun
function append_sumtag(df; data_col=:signal)
    grp_cols = names(df)
    deleteat!(grp_cols, findfirst(grp_cols, data_col))
    deleteat!(grp_cols, findfirst(grp_cols, :mstag))
    sum_df = map(groupby(df, grp_cols)) do df
        sum_df = df[1:1, :]
        sum_df[!, :mstag] .= "Sum"
        sum_df[!, data_col] .= sum(df[!, data_col])
        sum_df
    end |> combine
    vcat(df, sum_df)
end

end
