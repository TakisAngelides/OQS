using ITensors
using LinearAlgebra
using HDF5
using SparseArrays
using Arpack
using KrylovKit
using TupleTools
using OpenQuantumTools
using Statistics
using Dates
using BenchmarkTools
using Plots
include("/lustre/fs24/group/cqta/tangelides/OQS/Closed/Utilities.jl")

# ---------------------------------------------------------------------------------------------------------------

let
    f = h5open("inputs.h5", "r")
    ff = open("scrap.txt", "w")
    counter = 0
    for job_id in collect(keys(f))
        group = f["$(job_id)"]
        group_attributes = attributes(group)
        inputs = Dict()
        for key in keys(group_attributes)
            inputs[key] = read(group_attributes[key])
        end
        N, ma, cutoff, wis, l_0_1, x = inputs["N"], inputs["ma"], inputs["cutoff"], inputs["wis"], inputs["l_0_1"], inputs["x"]
        if wis == "dirac_vacuum" # this is always 1 before
            continue   
        end
        if N == 12 && cutoff == 1e-9 && ma == 0.0 && x == 4.0
            counter += 1
            println(counter, ": ", job_id)
        end
    end
end

# ---------------------------------------------------------------------------------------------------------------

