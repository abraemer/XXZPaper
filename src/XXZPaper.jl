module XXZPaper

# Write your package code here.

using Statistics
using SimLib
using JLD2: jldopen

function entropy_to_average(dir)
    dir = joinpath(dir, "entropy")
    files = readdir(dir)
    for file in files
        data = jldopen(joinpath(dir, file))["data"]
        newdata = HCEData(data.descriptor, mean(data.data; dims=1))
        newdata.descriptor.derivedfrom.pathdata.suffix = "avg"
        save(newdata)
    end
end

function merge_runs(dir)
    run_folders = filter(isdir, joinpath.(abspath(dir), readdir(dir)))
    outdir = joinpath(abspath(dir), "merged")

    subdirs = ["entropy", "locality", "lsr", "pr"]
    for subdir in subdirs
        subfiles = readdir(joinpath(run_folders[1], subdir))
        for subfile in subfiles
            data = load(joinpath(run_folders[1], subdir, subfile))
            for run_folder in run_folders[2:end]
                path = joinpath(run_folder, subdir, subfile)
                isfile(path) || continue
                data = merge_data(data, load(path))
            end
            save(data; prefix=outdir)
        end
    end
end

function merge_data(data1::T, data2::T) where {T<:SimLib.AbstractSimpleData}
    return T(data1.descriptor, cat(data1.data,data2.data; dims=_mergeaxis(data1)))
end

function _mergeaxis end
_mergeaxis(::Union{LSRData,ELData,PRData}) = 3
_mergeaxis(::HCEData) = 4
end
