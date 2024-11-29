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

    h5open("/lustre/fs24/group/cqta/tangelides/OQS/OQS_Purified/DAGS/50/inputs.h5", "r") do inputs_file
        
    counter = 0
    for j in 1:2:320

        group = inputs_file["$(j)"]
        group_attributes = attributes(group)
        inputs = Dict()
        for key in keys(group_attributes)
            inputs[key] = read(group_attributes[key])
        end
        D, m, l_0, T = inputs["aD"], inputs["ma"], inputs["l_0_1"], inputs["aT"]
        if !(T == 7.0 || T == 31.473684210526315 || T == 55.94736842105263 || T == 80.42105263157895 || T == 100.0) || D != 5.0 || m != 1.0 || l_0 != 0.5

            continue
        
        else

            counter += 1
            println(counter, ": D = ", D, ", m = ", m, ", l_0 = ", l_0, ", T = ", T)
            flush(stdout)
        
            h5open("/lustre/fs24/group/cqta/tangelides/OQS/OQS_Purified/DAGS/50/HDF5/$(j).h5", "r") do f1
            h5open("/lustre/fs24/group/cqta/tangelides/OQS/OQS_Purified/DAGS/50/HDF5/$(j+1).h5", "r") do f2

            steps = 2000
            N = 12
            mutual_info_list = Float64[]
            for i in 1:steps
                t = time()
                mps1 = read(f1, "$(i)", MPS)
                mps2 = read(f2, "$(i)", MPS)
                push!(mutual_info_list, real(get_mutual_info(mps2, [div(N, 2) - 2, div(N, 2) - 1, div(N, 2)], [div(N, 2) + 1, div(N, 2) + 2, div(N, 2) + 3]) - get_mutual_info(mps1, [div(N, 2) - 2, div(N, 2) - 1, div(N, 2)], [div(N, 2) + 1, div(N, 2) + 2, div(N, 2) + 3])))
                println(i, ": ", time() - t)
                flush(stdout)
            end

            h5open("Plots/mutual_info_data_50_extra/$(D)_$(m)_$(l_0)_$(T)_2_sites.h5", "w") do f
                write(f, "mutual_info_list", mutual_info_list)
            end

            end
            end

        end

    end

    end

end
