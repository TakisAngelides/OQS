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
include("/lustre/fs24/group/cqta/tangelides/OQS/OQS_Purified/utilities.jl")
ITensors.disable_warn_order()

let

    t = time()

    h5open("/lustre/fs24/group/cqta/tangelides/OQS/OQS_Purified/DAGS/46/inputs.h5", "r") do inputs_file
        
    for j in 1:2:16

        h5open("/lustre/fs24/group/cqta/tangelides/OQS/OQS_Purified/DAGS/46/HDF5/$(j).h5", "r+") do f1
        h5open("/lustre/fs24/group/cqta/tangelides/OQS/OQS_Purified/DAGS/46/HDF5/$(j+1).h5", "r+") do f2
    
        group = inputs_file["$(j)"]
        group_attributes = attributes(group)
        inputs = Dict()
        for key in keys(group_attributes)
            inputs[key] = read(group_attributes[key])
        end
        D, m, l_0 = inputs["aD"], inputs["ma"], inputs["l_0_1"]
        println(D, " ", m, " ", l_0)
        flush(stdout)
        
        steps = 2000
        N = 12
        mutual_info_list = Float64[]
        for i in 1:steps
            println(i, ": ", time() - t)
            flush(stdout)
            mps1 = read(f1, "$(i)", MPS)
            mps2 = read(f2, "$(i)", MPS)
            push!(mutual_info_list, real(get_mutual_info(mps2, [div(N, 2) - 1], [div(N, 2) + 2]) - get_mutual_info(mps1, [div(N, 2) - 1], [div(N, 2) + 2])))
        end

        h5open("Plots/mutual_info_data/$(D)_$(m)_$(l_0).h5", "w") do f
            write(f, "mutual_info_list", mutual_info_list)
        end

        end
        end
    
    end

    end

end
