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

# ---

# let
#     f = h5open("inputs.h5", "r")
#     ff = open("scrap.txt", "w")
#     counter = 0
#     for job_id in collect(keys(f))
#         group = f["$(job_id)"]
#         group_attributes = attributes(group)
#         inputs = Dict()
#         for key in keys(group_attributes)
#             inputs[key] = read(group_attributes[key])
#         end
#         N, ma, aD, cutoff, wis, l_0_1 = inputs["N"], inputs["ma"], inputs["aD"], inputs["cutoff"], inputs["wis"], inputs["l_0_1"]
#         if wis == "dirac_vacuum" # this is always 1 before
#             continue   
#         end
#         if N == 12 && cutoff == 1e-11 && ma == 0.0
#             counter += 1
#             println(counter, ": ", job_id)
#         end
#     end
# end

# Checking a saved state from a simulation how its energy is changing with time

let

ff = h5open("inputs.h5", "r")

for file_name in readdir("/lustre/fs24/group/cqta/tangelides/OQS/OQS_Purified/DAGS/18/specific_states")
    
    if occursin("h5", file_name) && !occursin("inputs.h5", file_name)

        group = ff["$(file_name[1:end-3])"]
        group_attributes = attributes(group)
        inputs = Dict()
        for key in keys(group_attributes)
            inputs[key] = read(group_attributes[key])
        end   

        wis = inputs["wis"]
        if wis == "dirac_vacuum_with_string"

            f = h5open("/lustre/fs24/group/cqta/tangelides/OQS/OQS_Purified/DAGS/18/specific_states/$(file_name)", "r+")
            steps1 = []
            states = []
            for k in collect(keys(f))
                if occursin("n", k)
                    continue
                else
                    push!(steps1, parse(Int, k))
                    push!(states, read(f, k, MPS))
                end
            end
            x = inputs["x"]
            l_0_1 = inputs["l_0_1"]
            ma = inputs["ma"]
            lambda = inputs["lambda"]
            aD = inputs["aD"]
            side = "left"
            energy1 = []
            for (state_idx, state) in enumerate(states)
                sites = siteinds(state)
                H = get_double_aH_Hamiltonian(sites, x, l_0_1, ma, lambda, side)
                H = MPO(H, sites)
                energy_value = measure_mpo(state, H)
                push!(energy1, energy_value)
            end
            # write(f, "energy", real(energy1))

            tmp = parse(Int, file_name[1:end-3])
            tmp -= 1
            f = h5open("/lustre/fs24/group/cqta/tangelides/OQS/OQS_Purified/DAGS/18/specific_states/$(tmp).h5", "r+")
            steps2 = []
            states = []
            for k in collect(keys(f))
                if occursin("n", k)
                    continue
                else
                    push!(steps2, parse(Int, k))
                    push!(states, read(f, k, MPS))
                end
            end
            x = inputs["x"]
            l_0_1 = inputs["l_0_1"]
            ma = inputs["ma"]
            lambda = inputs["lambda"]
            aD = inputs["aD"]
            side = "left"
            energy2 = []
            for (state_idx, state) in enumerate(states)
                sites = siteinds(state)
                H = get_double_aH_Hamiltonian(sites, x, l_0_1, ma, lambda, side)
                H = MPO(H, sites)
                energy_value = measure_mpo(state, H)
                push!(energy2, energy_value)
            end
            # write(f, "energy", real(energy2))

            perm1 = sortperm(steps1)
            perm2 = sortperm(steps2)
            steps_energy1, energy1 = steps1[perm1], energy1[perm1]
            steps_energy2, energy2 = steps1[perm2], energy2[perm2]
            p = plot()
            plot!(steps_energy1, real(energy1), label = "With string")
            plot!(steps_energy2, real(energy2), label = "No string")
            title!("aD = $(aD), l_0 = $(l_0_1)")
            savefig("/lustre/fs24/group/cqta/tangelides/OQS/OQS_Purified/DAGS/18/plots_specific_states/$(file_name).png")

        end

    end

end

end

# ------------------------------------------------------------------------------------------------------------------------------------