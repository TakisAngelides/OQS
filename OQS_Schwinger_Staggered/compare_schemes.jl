using ITensors
using LinearAlgebra
using HDF5
using Plots
using SparseArrays
using Arpack
using KrylovKit
using TupleTools
using OpenQuantumTools
using Statistics
include("Utilities.jl")
ITensors.disable_warn_order()

function run_ATDDMRG(scheme)

    t = time()
    
    sites = siteinds("S=1/2", N, conserve_qns = true)
    rho = get_dirac_vacuum_density_matrix(sites)
    
    # Get the exponential of the Lindblad terms for evolution 
    odd_gates_2 = get_exp_Ho_list(sites, -1im*tau/2) # odd/2
    odd_gates = get_exp_Ho_list(sites, -1im*tau) # odd
    even_gates = get_exp_He_list(sites, -1im*tau) # even

    # Prepare the lists for the tracked observables and the associated MPO
    max_bond_list = Int[]
    avg_bond_list = Float64[]
    ee_list = Float64[]
    step_num_list = Int[]
    z_mpo = [MPO(get_Z_site_operator(idx), sites) for idx in 1:N]
    z_list = [[] for _ in 1:N]
    particle_number_mpo = get_particle_number_MPO(sites)
    particle_number = [real(tr(apply(rho, particle_number_mpo)))]
    if get_state_diff_norm
        state_list = [project_zeroq(mpo_to_matrix(rho))]
    end

    # Push into the lists the initial state observables
    push!(step_num_list, 0)
    push!(max_bond_list, maxlinkdim(rho))
    push!(avg_bond_list, mean(linkdims(rho)))
    for idx in 1:N
        push!(z_list[idx], real(tr(apply(rho, z_mpo[idx]))))
    end
    rhodim = Int(binomial(N, div(N, 2)))
    # push!(ee_list, get_entanglement_entropy_reduced_matrix(N, reshape(project_zeroq(mpo_to_matrix(rho)), rhodim, rhodim)))
    if get_entanglement
        push!(ee_list, get_entanglement_entropy_mpo(rho, 1:div(N, 2), sites))
    end
    
    # Print statements
    println("The time to get the initial rho and MPO lists is: $(time() - t)\n")
    println("Now starting the ATDDMRG algorithm\n")
    flush(stdout)

    # test
    # rho = contract(rho)

    t = time()
    for step in 1:max_steps

        println("The step is ", step)

        if scheme == "simple"
            apply_odd!(odd_gates, rho; cutoff = cutoff, max_rho_D = max_rho_D)
            rho = apply_taylor_part(rho, tau, sites, x, l_0, ma, aD_0, sigma_over_a, env_corr_type, aT; cutoff = cutoff, max_rho_D = max_rho_D)
            apply_even!(even_gates, rho; cutoff = cutoff, max_rho_D = max_rho_D)
        else
            if step == 1
                apply_odd!(odd_gates_2, rho; cutoff = cutoff, max_rho_D = max_rho_D)
                rho = apply_taylor_part(rho, tau, sites, x, l_0, ma, aD_0, sigma_over_a, env_corr_type, aT; cutoff = cutoff, max_rho_D = max_rho_D)
                apply_even!(even_gates, rho; cutoff = cutoff, max_rho_D = max_rho_D)
                rho = apply_taylor_part(rho, tau, sites, x, l_0, ma, aD_0, sigma_over_a, env_corr_type, aT; cutoff = cutoff, max_rho_D = max_rho_D)           
            elseif step == max_steps            
                apply_odd!(odd_gates, rho; cutoff = cutoff, max_rho_D = max_rho_D)
                rho = apply_taylor_part(rho, tau, sites, x, l_0, ma, aD_0, sigma_over_a, env_corr_type, aT; cutoff = cutoff, max_rho_D = max_rho_D)
                apply_even!(even_gates, rho; cutoff = cutoff, max_rho_D = max_rho_D)
                rho = apply_taylor_part(rho, tau, sites, x, l_0, ma, aD_0, sigma_over_a, env_corr_type, aT; cutoff = cutoff, max_rho_D = max_rho_D)
                apply_odd!(odd_gates_2, rho; cutoff = cutoff, max_rho_D = max_rho_D)                 
            else    
                apply_odd!(odd_gates, rho; cutoff = cutoff, max_rho_D = max_rho_D)
                rho = apply_taylor_part(rho, tau, sites, x, l_0, ma, aD_0, sigma_over_a, env_corr_type, aT; cutoff = cutoff, max_rho_D = max_rho_D)
                apply_even!(even_gates, rho; cutoff = cutoff, max_rho_D = max_rho_D)
                rho = apply_taylor_part(rho, tau, sites, x, l_0, ma, aD_0, sigma_over_a, env_corr_type, aT; cutoff = cutoff, max_rho_D = max_rho_D)
            end

        end
            
        # Take care of hermiticity and positivity of the density matrix
        rho = add(dag(swapprime(rho, 0, 1)), rho; cutoff = cutoff, maxdim = max_rho_D)/2 # fix hermiticity with rho -> rho dagger + rho over 2
        rho = rho/tr(rho)
        
        if (step % measure_every == 0) || (step == max_steps)

            println("Now measuring observables")

            println("The link dimensions are: $(linkdims(rho))")

            # println("The L_taylor*dt/2 should be much less than one: $(real(inner(L_taylor_expanded_part_tmp, rho))*tau/2)\n")

            # println("The trace should always be 1: ", tr(rho))

            # Measure the observables
            push!(step_num_list, step)
            push!(avg_bond_list, mean(linkdims(rho)))
            push!(particle_number, real(tr(apply(rho, particle_number_mpo))))
            push!(max_bond_list, maxlinkdim(rho))
            if get_state_diff_norm
                push!(state_list, project_zeroq(mpo_to_matrix(rho)))
            end

            # println("The trace of rho should be 1: ", tr(rho))
            
            # test
            # push!(state_list, project_zeroq(reshape(Array(rho, inds(rho; :plev => 1)..., inds(rho; :plev => 0)...), 2^N, 2^N)))

            # sum_for_z = 0

            for idx in 1:N
                exp_val_z = real(tr(apply(rho, z_mpo[idx])))
                # exp_val_z = real(tr(project_zeroq(mpo_to_matrix(rho))*project_zeroq(mpo_to_matrix(z_mpo[idx]))))
                # sum_for_z += exp_val_z
                push!(z_list[idx], exp_val_z)
            end

            # println("The total charge should be 0: ", sum(sum_for_z))
            # push!(ee_list, get_entanglement_entropy_mpo(rho, div(N, 2)+1:N, sites))

            # push!(ee_list, get_entanglement_entropy_reduced_matrix(N, reshape(project_zeroq(mpo_to_matrix(rho)), rhodim, rhodim)))
            if get_entanglement
                push!(ee_list, get_entanglement_entropy_mpo(rho, 1:div(N, 2), sites))
            end

            # test
            # rhodim = Int(binomial(N, div(N, 2)))
            # push!(ee_list, get_entanglement_entropy_reduced_matrix(N, reshape(project_zeroq(reshape(Array(rho, inds(rho; :plev => 1)..., inds(rho; :plev => 0)...), 2^N, 2^N)), rhodim, rhodim)))

            # Refresh time variable and print statement
            println("Step: $step finished, Time = $(time()-t), Average Step Time = $((time() - t)/measure_every)\n")            
            flush(stdout)
            t = time()

        end

    end

    # Prepare the initial state and the Z observable at the middle of the lattice to be tracked
    if get_sparse
        particle_number_sparse_operator = project_zeroq(get_particle_number_operator_sparse(N))
        z_op = [project_zeroq(get_op(["Z"], [idx], N)) for idx in 1:N]
        
        # rho = outer(gs', gs; cutoff = 0)
        # rho_m = mpo_to_matrix(rho)
        # rho_m = project_zeroq(rho_m)
    
        rho_m = get_dirac_vacuum_zeroq_density_matrix_sparse(N)
        
        rho_v = reshape(rho_m, length(rho_m))
        if get_state_diff_norm
            state_sparse_list = [rho_m]
        end

        # Get the Lindblad operator and its exponential which is the evolution operator
        L = get_Lindblad_reduced_sparse_matrix(N, x, ma, l_0, lambda, aD_0, sigma_over_a, aT, env_corr_type)
        evolution_operator = exp(Matrix(L)*tau)
        
        # Prepare the list to store the tracked observables and get the initial state values
        rhodim = Int(binomial(N, div(N, 2)))
        step_num_sparse_list = Int[]
        push!(step_num_sparse_list, 0)
        ee_list_sparse = [get_entanglement_entropy_reduced_matrix(N, reshape(rho_v, rhodim, rhodim))]
        z_list_sparse = [[] for _ in 1:N]
        for idx in 1:N
            push!(z_list_sparse[idx], real(tr(rho_m*z_op[idx])))
        end
        particle_number_sparse = [real(tr(particle_number_sparse_operator*rho_m))]
        
        # Do the evolution
        for step in 1:max_steps

            # One time step evolution
            rho_v = evolution_operator*rho_v

            # Measure tracked observables
            push!(ee_list_sparse, get_entanglement_entropy_reduced_matrix(N, reshape(rho_v, rhodim, rhodim)))
            push!(step_num_sparse_list, step)
            rho_m = reshape(rho_v, (rhodim, rhodim))
            if get_state_diff_norm
                push!(state_sparse_list, rho_m)
            end
            for idx in 1:N
                push!(z_list_sparse[idx], real(tr(rho_m*z_op[idx])))
            end
            push!(particle_number_sparse, real(tr(particle_number_sparse_operator*rho_m)))

        end

        return max_bond_list, avg_bond_list, ee_list, step_num_list, particle_number, particle_number_sparse

    else 

        return max_bond_list, avg_bond_list, ee_list, step_num_list, particle_number

    end

end

N = 6
tau = 0.01 # 1/N^2 # time step in time evolution rho -> exp(-tau L) after one step
cutoff = 1e-16 # cutoff for SVD
max_rho_D = 1000
max_steps = 100
tol = 1e-16 # tolerance for DMRG convergence
e = 0.8
x = 1/(e)^2
ma = 0.5
l_0 = 0.2
lambda = 0.0
aD_0 = 1
beta = 0.1
aT = 1/beta
sigma_over_a = 3.0
env_corr_type = "delta"
max_sweeps = 1000
l_0_initial = 0.0
measure_every = 1 # this determines how often to save rho and measure the energy in ATDDMRG
get_entanglement = true
get_state_diff_norm = false
scheme_list = ["simple", "2nd order"]

p1 = plot(title = "Particle Number Density vs time step number\nN=$(N),tau=$(tau),cutoff=$(cutoff),aD=$(aD_0)\nbeta=$(beta),e=$(e),l_0=$(l_0),ma=$(ma)", legend = true)

let

for scheme in scheme_list

    if scheme == "simple"
        get_sparse = true
    else
        get_sparse = false
    end

    results = run_ATDDMRG(scheme)

    if get_sparse
        steps, pnd_sp, pnd = results[4], results[6], results[5]
        plot!(p1, steps, pnd_sp, label = "sparse", linestyle = :auto)
    else
        steps, pnd = results[4], results[5]
    end

    plot!(p1, steps, pnd, label = scheme, linestyle = :auto)

end

end

display(p1)


