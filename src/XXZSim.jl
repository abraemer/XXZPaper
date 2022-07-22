module XXZSim

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

end
