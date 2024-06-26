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

N = 14
tau = 0.1 # 1/N^2 # time step in time evolution rho -> exp(-tau L) after one step
cutoff = 1e-16 # cutoff for SVD
max_rho_D = 100
max_steps = 3
tol = 1e-9 # tolerance for DMRG convergence and ATDDMRG convergence
e = 0.8
x = 1/(e)^2
ma = 0.5
l_0 = 0.0
lambda = 0.0
aD_0 = 1
beta = 0.001
aT = 1/beta
sigma_over_a = 3.0
env_corr_type = "delta"
max_sweeps = 1000
l_0_initial = 0.0
measure_every = 1 # this determines how often to save rho and measure the energy in ATDDMRG
get_entanglement = false
get_state_diff_norm = false
get_sparse = false

function run_ATDDMRG()

    t = time()
    
    # Prepare initial rho
    # println("Initializing with the MPO corresponding to the ground state of H_system\n")
    sites = siteinds("S=1/2", N, conserve_qns = true)
    # state = [isodd(n) ? "0" : "1" for n = 1:N]
    # mps = randomMPS(sites, state)
    # H = get_aH_Hamiltonian(sites, x, l_0_initial, ma, lambda)
    # sweeps = Sweeps(max_steps; maxdim = D)
    # observer = DMRGObserver(;energy_tol = tol)
    # gs_energy, gs = dmrg(H, mps, sweeps; outputlevel = 1, observer = observer, ishermitian = true)
    # println("The ground state energy was found to be $(gs_energy)\n")
    # rho = outer(gs', gs; cutoff = cutoff)

    # Test Dirac vacuum as initial state
    # state = [isodd(n) ? "1" : "0" for n = 1:N]
    # mps = MPS(sites, state)
    # rho = outer(mps', mps; cutoff = cutoff)
    rho = get_dirac_vacuum_density_matrix(sites)

    L_taylor_expanded_part_tmp = get_L_taylor(sites, x, l_0, ma, aD_0, sigma_over_a, env_corr_type, aT)
    println("The expectation value of the taylor expanded part of the Lindblad operator is $(real(inner(L_taylor_expanded_part_tmp, rho)))\n")
    println("The L_taylor*dt/2 should be much less than one: $(real(inner(L_taylor_expanded_part_tmp, rho))*tau/2)\n")
    println("The trace of initial rho should be 1: $(tr(rho))\n")
    flush(stdout)
    
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

        # Take care of hermiticity and positivity of the density matrix
        rho = add(dag(swapprime(rho, 0, 1)), rho; cutoff = cutoff, maxdim = max_rho_D)/2 # fix hermiticity with rho -> rho dagger + rho over 2
        rho = rho/tr(rho)
        
        if (step % measure_every == 0) || (step == max_steps)

            println("Now measuring observables")

            println("The link dimensions are: $(linkdims(rho))")

            # println("The L_taylor*dt/2 should be much less than one: $(real(inner(L_taylor_expanded_part_tmp, rho))*tau/2)\n")

            # println("The trace should always be 1: ", tr(rho))

            # Measure the observables
            # push!(step_num_list, step)
            # push!(avg_bond_list, mean(linkdims(rho)))
            push!(particle_number, real(tr(apply(rho, particle_number_mpo))))
            # push!(max_bond_list, maxlinkdim(rho))
            if get_state_diff_norm
                push!(state_list, project_zeroq(mpo_to_matrix(rho)))
            end

            # println("The trace of rho should be 1: ", tr(rho))
            
            # test
            # push!(state_list, project_zeroq(reshape(Array(rho, inds(rho; :plev => 1)..., inds(rho; :plev => 0)...), 2^N, 2^N)))

            # sum_for_z = 0

            # for idx in 1:N
            #     exp_val_z = real(tr(apply(rho, z_mpo[idx])))
            #     # exp_val_z = real(tr(project_zeroq(mpo_to_matrix(rho))*project_zeroq(mpo_to_matrix(z_mpo[idx]))))
            #     # sum_for_z += exp_val_z
            #     push!(z_list[idx], exp_val_z)
            # end

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

    println("The final ee is: ", get_entanglement_entropy_mpo(rho, 1:div(N, 2), sites))

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
    end
    
    p1 = plot()
    plot!(max_bond_list, label = "Max")
    plot!(avg_bond_list, label = "Average")
    title!("Maximum/Average bond dimension vs step number")
    display(p1)

    p5 = plot()
    # particle_number_from_z = []
    # for t_idx in 1:max_steps+1
    #     push!(particle_number_from_z, 0.5*N + sum([z_list[i][t_idx]*(-1)^(i-1)*0.5 for i in 1:N]))
    # end
    title!("Particle Number vs step number\nN = $(N), tau = $(tau), cutoff = $(cutoff)")
    # plot!(p5, particle_number_from_z, label = "from Z")
    plot!(p5, particle_number, label = "direct", linestyle = :dash)
    if get_sparse
        plot!(p5, particle_number_sparse, label = "sparse")
    end
    display(p5)

    if get_entanglement
        p2 = plot()
        plot!(step_num_list[1:end-1], abs.(ee_list - ee_list_sparse)[1:end-1])
        # plot!(step_num_sparse_list, ee_list_sparse, label = "Sparse")
        # plot!(step_num_list, ee_list, label = "MPO")
        println("Final entropy: ", ee_list[end], " ", ee_list_sparse[end])
        title!("Difference in Entanglement Entropy MPO vs Sparse")
        display(p2)
    end

    # p3 = plot()
    # for idx in 1:N
    #     # plot!(step_num_list, z_list[idx], label = "MPO, $(idx)")
    #     # plot!(step_num_sparse_list, z_list_sparse[idx], label = "Sparse, $(idx)", linestyle = :dash)
    #     plot!(step_num_list[1:end-1], abs.((z_list[idx] - z_list_sparse[idx])[1:end-1]), label = "Site: $(idx)")
    # end
    # title!("Difference of Middle Z Operator\nExpectation Value MPO vs Sparse")
    # display(p3)

    if get_state_diff_norm
        norm_diff = []
        for idx in 1:max_steps+1
            e1 = sort(real(eigen(state_list[idx]).values))
            e2 = sort(real(eigen(state_sparse_list[idx]).values))
            # push!(norm_diff, norm(state_list[idx] - state_sparse_list[idx]))
            push!(norm_diff, norm(e1 - e2))
        end
        p4 = plot()
        plot!(step_num_list, norm_diff)
        title!("Norm of difference of MPO and sparse states")
        display(p4)
    end

end

run_ATDDMRG()
