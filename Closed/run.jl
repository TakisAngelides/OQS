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
include("Utilities.jl")

# Input arguments for file and opening the results h5 file
println("Now getting the input arguments for the file and opening the results h5 file ", now())
flush(stdout)
path_to_project_number = ARGS[1] 
path_to_inputs_h5 = "$(path_to_project_number)/inputs.h5" # Where the h5 file containing the inputs is
job_id = ARGS[2] # The job id to get the inputs within this h5 file for this specific job
path_to_results = "$(path_to_project_number)/HDF5/$(job_id).h5" # The path to the h5 file to save all results of states
results_file = h5open(path_to_results, "w")
println("Finished getting the input arguments for the file and opening the results h5 file ", now())
flush(stdout)

# Build a dictionary of the inputs
println("Now building a dictionary of the inputs", now())
flush(stdout)
inputs_file = h5open(path_to_inputs_h5, "r")
group = inputs_file["$(job_id)"]
group_attributes = attributes(group)
inputs = Dict()
for key in keys(group_attributes)
    inputs[key] = read(group_attributes[key])
end
close(inputs_file)
println("The inputs are ", inputs)
flush(stdout)
println("Finished building a dictionary of the inputs ", now())
flush(stdout)

# Prepare the initial state
println("Now getting the initial state ", now())
flush(stdout)
function get_initial_state(which_initial_state, conserve_qns, N, inputs, results_file, sites)

    if which_initial_state == "dirac_vacuum"
    
        return get_dirac_vacuum_mps(sites)

    elseif which_initial_state == "gs_naive" # This is naive because it takes the outer product instead of doing optimization so bond dimension will be D^2

        x = inputs["x"]
        l_0_initial_state = inputs["l_0_initial_state"]
        ma = inputs["ma"]
        lambda = inputs["lambda"]
        max_sweeps_dmrg = inputs["msdmrg"]
        maxdim_dmrg = inputs["mddmrg"]
        energy_tol_dmrg = inputs["etdmrg"]
        cutoff_dmrg = inputs["cdmrg"]
        state = [isodd(n) ? "0" : "1" for n = 1:N]
        psi = MPS(sites, state)
        H = get_aH_Hamiltonian(sites, x, l_0_initial_state, ma, lambda)
        sweeps = Sweeps(max_sweeps_dmrg, maxdim = maxdim_dmrg, cutoff = cutoff_dmrg)
        observer = DMRGObserver(;energy_tol = energy_tol_dmrg)
        gs_energy, gs = dmrg(H, psi, sweeps; outputlevel = 1, observer = observer, ishermitian = true)
        write(results_file, "gs_energy", gs_energy)
        write(results_file, "gs", gs)
        
        return gs

    else # case of which_initial_state = "dirac_vacuum_with_string"

        flip_sites = inputs["fs"]
        psi = get_dirac_vacuum_mps(sites; flip_sites) # This is a normal non-purified MPS
        
        return psi

    end

end
which_initial_state = inputs["wis"]
conserve_qns = parse(Bool, inputs["cqns"])
N = inputs["N"]
sites = siteinds("S=1/2", N, conserve_qns = conserve_qns)
mps = get_initial_state(which_initial_state, conserve_qns, N, inputs, results_file, sites)
orthogonalize!(mps, 1) # put the MPS in right canonical form to start the ATDDMRG as expected by the functions
println("Finished getting the initial state ", now())
flush(stdout)

# Get the taylor, odd and even opsum groups
x = inputs["x"]
t_over_a = 0 # starting the time variable
which_applied_field = inputs["waf"]
l_0 = get_applied_field(which_applied_field, inputs, t_over_a)
ma = inputs["ma"]
tau = inputs["tau"]
Ho_mpo_list = get_exp_Ho(sites, -1im*tau/2, x, ma, l_0) # odd/2
He_mpo_list = get_exp_He(sites, -1im*tau, x, ma, l_0) # even
Ho_mpo_list_2 = get_exp_Ho(sites, -1im*tau, x, ma, l_0) # odd
taylor_mpo_tmp = get_exp_Hz(sites, -1im*tau/2, x) # H_z/2
taylor_expansion_order = inputs["teo"]
taylor_expansion_cutoff_1 = inputs["tec_1"]
truncate!(taylor_mpo_tmp; cutoff = taylor_expansion_cutoff_1)
taylor_expansion_cutoff_2 = inputs["tec_2"]
taylor_mpo = get_mpo_taylor_expansion(taylor_mpo_tmp, taylor_expansion_order, taylor_expansion_cutoff_2, sites)

# Starting the lists of the observables we want to keep track of
println("Now getting the lists for the tracked observables ", now())
flush(stdout)
number_of_time_steps = inputs["nots"]
z_configs = zeros(ComplexF64, number_of_time_steps+1, N)
z_configs[1, :] = expect(mps, "Z")
link_dims = zeros(Int64, number_of_time_steps+1, N-1)
link_dims[1, :] = linkdims(mps)
println("Finished getting the lists for the tracked observables ", now())
flush(stdout)

# Perform the time evolution
function evolve(which_applied_field, Ho_mpo_list, He_mpo_list, Ho_mpo_list_2, taylor_mpo, results_file, inputs, z_configs, mps)

    cutoff = inputs["cutoff"]
    maxdim = inputs["md"]
    time_varying_applied_field_flag = parse(Bool, inputs["tvaff"])
    which_steps_to_save_state = inputs["wstss"]

    if time_varying_applied_field_flag

        t_over_a = 0
        l_0_list = [get_applied_field(which_applied_field, inputs, t_over_a)]

        for step in 1:number_of_time_steps

            t = time() # Starting the time for the step
    
            # Recalculate only the gates which are affected by l_0 and save the new l_0 to a list
            t_over_a += tau # Incrementing the t_over_a time variable as it was initiated to 0 and observables were measured
            l_0 = get_applied_field(which_applied_field, inputs, t_over_a)
            push!(l_0_list, l_0)
            Ho_mpo_list = get_exp_Ho(sites, -1im*tau/2, x, ma, l_0) # odd/2
            He_mpo_list = get_exp_He(sites, -1im*tau, x, ma, l_0) # even
            Ho_mpo_list_2 = get_exp_Ho(sites, -1im*tau, x, ma, l_0) # odd

            # One time step with ATDDMRG
            if step == 1
                apply_Ho_mpo_list!(Ho_mpo_list, mps; cutoff = cutoff, maxdim = maxdim)
                mps = apply(taylor_mpo, mps; cutoff = cutoff, normalize = true, maxdim = maxdim)
                apply_He_mpo_list!(He_mpo_list, mps; cutoff = cutoff, maxdim = maxdim)
                mps = apply(taylor_mpo, mps; cutoff = cutoff, normalize = true, maxdim = maxdim)
            elseif step == number_of_time_steps
                apply_Ho_mpo_list!(Ho_mpo_list_2, mps; cutoff = cutoff, maxdim = maxdim)
                mps = apply(taylor_mpo, mps; cutoff = cutoff, normalize = true, maxdim = maxdim)
                apply_He_mpo_list!(He_mpo_list, mps; cutoff = cutoff, maxdim = maxdim)
                mps = apply(taylor_mpo, mps; cutoff = cutoff, normalize = true, maxdim = maxdim)
                apply_Ho_mpo_list!(Ho_mpo_list, mps; cutoff = cutoff, maxdim = maxdim)    
            else
                apply_Ho_mpo_list!(Ho_mpo_list_2, mps; cutoff = cutoff, maxdim = maxdim)
                mps = apply(taylor_mpo, mps; cutoff = cutoff, normalize = true, maxdim = maxdim)
                apply_He_mpo_list!(He_mpo_list, mps; cutoff = cutoff, maxdim = maxdim)
                mps = apply(taylor_mpo, mps; cutoff = cutoff, normalize = true, maxdim = maxdim)
            end

            # Compute the tracked observables
            z_configs[step+1, :] = expect(mps, "Z")
            linkdims_of_step = linkdims(mps)
            link_dims[step+1, :] = linkdims_of_step

            # Save state to file
            if step in which_steps_to_save_state
                write(results_file, "$(step)", mps)
            end

            println("Step = $(step), Time = $(time() - t), Links = $(linkdims_of_step)")
            flush(stdout)

        end

        # Write tracked observables to results h5 file
        println("Now writing the observables to results h5 file ", now())
        flush(stdout)
        write(results_file, "z_configs", z_configs)
        write(results_file, "link_dims", link_dims)
        write(results_file, "l_0_list", l_0_list)
        println("Finished writing the observables to results h5 file ", now())
        flush(stdout)

    else

        for step in 1:number_of_time_steps

            t = time() # Starting the time for the step
    
            # One time step with ATDDMRG
            if step == 1
                apply_Ho_mpo_list!(Ho_mpo_list, mps; cutoff = cutoff, maxdim = maxdim)
                mps = apply(taylor_mpo, mps; cutoff = cutoff, normalize = true, maxdim = maxdim)
                apply_He_mpo_list!(He_mpo_list, mps; cutoff = cutoff, maxdim = maxdim)
                mps = apply(taylor_mpo, mps; cutoff = cutoff, normalize = true, maxdim = maxdim)
            elseif step == number_of_time_steps
                apply_Ho_mpo_list!(Ho_mpo_list_2, mps; cutoff = cutoff, maxdim = maxdim)
                mps = apply(taylor_mpo, mps; cutoff = cutoff, normalize = true, maxdim = maxdim)
                apply_He_mpo_list!(He_mpo_list, mps; cutoff = cutoff, maxdim = maxdim)
                mps = apply(taylor_mpo, mps; cutoff = cutoff, normalize = true, maxdim = maxdim)
                apply_Ho_mpo_list!(Ho_mpo_list, mps; cutoff = cutoff, maxdim = maxdim)    
            else
                apply_Ho_mpo_list!(Ho_mpo_list_2, mps; cutoff = cutoff, maxdim = maxdim)
                mps = apply(taylor_mpo, mps; cutoff = cutoff, normalize = true, maxdim = maxdim)
                apply_He_mpo_list!(He_mpo_list, mps; cutoff = cutoff, maxdim = maxdim)
                mps = apply(taylor_mpo, mps; cutoff = cutoff, normalize = true, maxdim = maxdim)
            end

            # Compute the tracked observables
            z_configs[step+1, :] = expect(mps, "Z")
            linkdims_of_step = linkdims(mps)
            link_dims[step+1, :] = linkdims_of_step

            # Save state to file
            if step in which_steps_to_save_state
                write(results_file, "$(step)", mps)
            end

            println("Step = $(step), Time = $(time() - t), Links = $(linkdims(mps))")
            flush(stdout)

        end

        # Write tracked observables to results h5 file
        println("Now writing the observables to results h5 file ", now())
        flush(stdout)
        write(results_file, "z_configs", z_configs)
        write(results_file, "link_dims", link_dims)
        println("Finished writing the observables to results h5 file ", now())
        flush(stdout)

    end

end
evolve(which_applied_field, Ho_mpo_list, He_mpo_list, Ho_mpo_list_2, taylor_mpo, results_file, inputs, z_configs, mps)

close(results_file)

println("Finished ", now())
flush(stdout)
