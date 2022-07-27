import Pkg
Pkg.activate(".")
using XXZPaper
using SimLib: path_prefix

XXZPaper.merge_runs(joinpath(path_prefix("xxzpaper"), "high_density-avg"))
XXZPaper.merge_runs(joinpath(path_prefix("xxzpaper"), "low_density-avg"))
