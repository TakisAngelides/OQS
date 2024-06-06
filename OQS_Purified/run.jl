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
include("utilities.jl")

# Input arguments for file and opening the results h5 file
println("Now getting the input arguments for the file and opening the results h5 file ", now())
flush(stdout)
path_to_inputs_h5 = ARGS[1] # Where the h5 file containing the inputs is
job_id = ARGS[2] # The job id to get the inputs within this h5 file for this specific job
path_to_save_results_h5 = ARGS[3] # The path to the h5 file to save all results
results_file = h5open(path_to_save_results_h5, "w")
write(results_file, "path_to_inputs_h5", path_to_inputs_h5) # This can be used later to get the constants such as N etc
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
function get_initial_state(which_initial_state, conserve_qns, N, results_file, inputs)

    sites_initial_state = siteinds("S=1/2", N, conserve_qns = conserve_qns)

    if which_initial_state == "dirac_vacuum"

        psi = get_dirac_vacuum_mps(sites_initial_state) # This is a normal non-purified MPS
        rho = outer(psi', psi) # Get the density matrix
        rho_vec = convert(MPS, rho) # Convert the density matrix to a purified MPS
        mps = rho_vec_to_mps(rho_vec) # Split the sites so that each site has one physical index of dimension 2
        orthogonalize!(mps, 1) # Bring the center of orthogonalization to the very left
        
        return mps

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
        psi = MPS(sites_initial_state, state)
        H = get_aH_Hamiltonian(sites_initial_state, x, l_0_initial_state, ma, lambda)
        sweeps = Sweeps(max_sweeps_dmrg, maxdim = maxdim_dmrg, cutoff = cutoff_dmrg)
        observer = DMRGObserver(;energy_tol = energy_tol_dmrg)
        gs_energy, gs = dmrg(H, psi, sweeps; outputlevel = 1, observer = observer, ishermitian = true)
        write(results_file, "gs_energy", gs_energy)
        write(results_file, "gs", gs)
        rho = outer(gs', gs) # Get the density matrix
        rho_vec = convert(MPS, rho) # Convert the density matrix to a purified MPS
        mps = rho_vec_to_mps(rho_vec) # Split the sites so that each site has one physical index of dimension 2
        orthogonalize!(mps, 1) # Bring the center of orthogonalization to the very left

        return mps

    else # case of which_initial_state = "dirac_vacuum_with_string"

        flip_sites = inputs["fs"]
        psi = get_dirac_vacuum_mps(sites_initial_state; flip_sites) # This is a normal non-purified MPS
        rho = outer(psi', psi) # Get the density matrix
        rho_vec = convert(MPS, rho) # Convert the density matrix to a purified MPS
        mps = rho_vec_to_mps(rho_vec) # Split the sites so that each site has one physical index of dimension 2
        orthogonalize!(mps, 1) # Bring the center of orthogonalization to the very left
        
        return mps

    end

end
which_initial_state = inputs["wis"]
conserve_qns = parse(Bool, inputs["cqns"])
N = inputs["N"]
mps = get_initial_state(which_initial_state, conserve_qns, N, results_file, inputs)
println("Finished getting the initial state ", now())
flush(stdout)

# Get the taylor, odd and even opsum groups without the l0 terms
println("Now getting the odd, even and taylor gates without the l0 terms ", now())
flush(stdout)
x = inputs["x"]
ma = inputs["ma"]
lambda = inputs["lambda"]
aT = inputs["aT"]
aD = inputs["aD"]
env_corr_type = inputs["env_corr_type"]
sites = siteinds(mps) # This is done so that the odd, even gates and taylor MPO have physical legs matching the purified MPS and combining this with the swapprime done on the operators later the transpose is taken on the operators acting on the even sites which correspond to operators acting on the right of the density matrix
for i in 2:2:length(sites)
    sites[i] = dag(sites[i])
end
tau = inputs["tau"]
dissipator_sites = inputs["ds"]
opsum_without_l0_terms = get_Lindblad_opsum_without_l0_terms(sites, x, ma, lambda, aT, aD, env_corr_type, inputs, dissipator_sites)
nn_odd_without_l0_terms, nn_even_without_l0_terms, taylor = get_odd_even_taylor_groups(opsum_without_l0_terms, sites)
println("Finished getting the odd, even and taylor gates without the l0 terms ", now())
flush(stdout)

# Get the odd and even opsum groups with just the l0 terms
println("Now getting the odd, even and taylor gates with just the l0 terms ", now())
flush(stdout)
which_applied_field = inputs["waf"]
t_over_a = 0 # starting the time variable
l_0 = get_applied_field(which_applied_field, inputs, t_over_a)
opsum_just_l0_terms = get_Lindblad_opsum_just_l0_terms(sites, x, l_0, lambda)
nn_odd_just_l0_terms, nn_even_just_l0_terms, _ = get_odd_even_taylor_groups(opsum_just_l0_terms, sites)
println("Finished getting the odd, even and taylor gates with just the l0 terms ", now())
flush(stdout)

# Gather the two odd and even gates
println("Now putting all the even and odd together ", now())
flush(stdout)
odd = get_odd(sites, tau/2, nn_odd_without_l0_terms .+ nn_odd_just_l0_terms)
even = get_even(sites, tau, nn_even_without_l0_terms .+ nn_even_just_l0_terms)
println("Finished putting all the even and odd together ", now())
flush(stdout)

# Get the MPO for the Taylor expansion
println("Now getting the MPO for the taylor expansion ", now())
flush(stdout)
taylor_expansion_order = inputs["teo"]
taylor_expansion_cutoff_1 = inputs["tec_1"]
taylor_mpo_tmp = 0.5*tau*MPO(taylor, sites)
truncate!(taylor_mpo_tmp; cutoff = taylor_expansion_cutoff_1)
for i in 2:2:length(taylor_mpo_tmp) # Transpose the MPO on the even sites which would correspond to the bottom legs of the MPO
    taylor_mpo_tmp[i] = swapprime(taylor_mpo_tmp[i], 0, 1; :tags => "Site")
end
taylor_expansion_cutoff_2 = inputs["tec_2"]
taylor_mpo = get_mpo_taylor_expansion(taylor_mpo_tmp, taylor_expansion_order, taylor_expansion_cutoff_2, sites)
println("The taylor_mpo with taylor order $(taylor_expansion_order) and cutoffs $(taylor_expansion_cutoff_1), $(taylor_expansion_cutoff_2) has bond dimensions ", linkdims(taylor_mpo))
flush(stdout)
println("Finished getting the MPO for the taylor expansion ", now())
flush(stdout)

# Starting the lists of the observables we want to keep track of
println("Now getting the lists for the tracked observables ", now())
flush(stdout)
number_of_time_steps = inputs["nots"]
z_configs = [[] for _ in 0:number_of_time_steps]
z_configs[1] = measure_z_config(mps)
println("Finished getting the lists for the tracked observables ", now())
flush(stdout)

# Perform the time evolution
function evolve(which_applied_field, odd, even, taylor_mpo, nn_odd_without_l0_terms, nn_even_without_l0_terms, results_file, inputs, z_configs, mps)

    cutoff = inputs["cutoff"]
    maxdim = inputs["md"]
    time_varying_applied_field_flag = parse(Bool, inputs["tvaff"])

    if time_varying_applied_field_flag

        t_over_a = 0
        l_0_list = [get_applied_field(which_applied_field, inputs, t_over_a)]

        for step in 1:number_of_time_steps

            t = time() # Starting the time for the step
    
            # Recalculate only the gates which are affected by l_0 and save the new l_0 to a list
            t_over_a += tau # Incrementing the t_over_a time variable as it was initiated to 0 and observables were measured
            l_0 = get_applied_field(which_applied_field, inputs, t_over_a)
            push!(l_0_list, l_0)
            opsum_just_l0_terms = get_Lindblad_opsum_just_l0_terms(sites, x, l_0, lambda)
            nn_odd_just_l0_terms, nn_even_just_l0_terms, _ = get_odd_even_taylor_groups(opsum_just_l0_terms, sites)
            odd = get_odd(sites, tau/2, nn_odd_without_l0_terms .+ nn_odd_just_l0_terms)
            even = get_even(sites, tau, nn_even_without_l0_terms .+ nn_even_just_l0_terms)

            # One time step with ATDDMRG
            apply_odd!(odd, mps, cutoff, maxdim)
            mps = apply(taylor_mpo, mps; cutoff = cutoff, maxdim = maxdim)
            apply_even!(even, mps, cutoff, maxdim)
            mps = apply(taylor_mpo, mps; cutoff = cutoff, maxdim = maxdim)
            apply_odd!(odd, mps, cutoff, maxdim)

            # Compute the tracked observables
            z_configs[step+1] = measure_z_config(mps)

            println("Step = $(step), Time = $(time() - t), Links = $(linkdims(mps))")
            flush(stdout)

        end

        write(results_file, "l_0_list", l_0_list)

    else

        for step in 1:number_of_time_steps

            t = time() # Starting the time for the step
    
            # One time step with ATDDMRG
            apply_odd!(odd, mps, cutoff, maxdim)
            mps = apply(taylor_mpo, mps; cutoff = cutoff, maxdim = maxdim)
            apply_even!(even, mps, cutoff, maxdim)
            mps = apply(taylor_mpo, mps; cutoff = cutoff, maxdim = maxdim)
            apply_odd!(odd, mps, cutoff, maxdim)

            # Compute the tracked observables
            z_configs[step+1] = measure_z_config(mps)

            println("Step = $(step), Time = $(time() - t), Links = $(linkdims(mps))")
            flush(stdout)

        end

    end

    # Write tracked observables to results h5 file
    println("Now writing the observables to results h5 file ", now())
    flush(stdout)
    z_configs_group = create_group(results_file, "z_configs_group")
    for step in 1:number_of_time_steps+1
        write_attribute(z_configs_group, "$(step)", z_configs[step])
    end
    println("Finished writing the observables to results h5 file ", now())
    flush(stdout)

end
evolve(which_applied_field, odd, even, taylor_mpo, nn_odd_without_l0_terms, nn_even_without_l0_terms, results_file, inputs, z_configs, mps)

close(results_file)

println("Finished ", now())
flush(stdout)
