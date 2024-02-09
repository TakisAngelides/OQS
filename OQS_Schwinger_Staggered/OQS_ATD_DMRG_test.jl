using ITensors
using LinearAlgebra
using HDF5
using Plots
using SparseArrays
using Arpack
using KrylovKit
using TupleTools
using OpenQuantumTools
include("Utilities.jl")

N = 4
tau = 0.001 # time step in time evolution rho -> exp(-tau L) after one step
cutoff = 1e-20 # cutoff for SVD
tol = 1e-9 # tolerance for DMRG convergence and ATDDMRG convergence
e = 0.8
x = 1/e^2
l_0 = 0.45
D = 1000
lambda = 0.0
ma = 0.5
aD_0 = 0.0
aT = 10.0
beta = 1/aT
sigma_over_a = 3.0
env_corr_type = "delta"
max_sweeps = 1000
max_steps = 3
measure_every = 1 # this determines how often to save rho and measure the energy in ATDDMRG

function run_ATDDMRG()

    t = time()
    
    # Prepare initial rho
    println("Initializing with the MPO corresponding to the ground state of H_system\n")
    sites = siteinds("S=1/2", N, conserve_qns = false)
    state = [isodd(n) ? "0" : "1" for n = 1:N]
    mps = randomMPS(sites, state)
    H = get_aH_Hamiltonian(sites, x, l_0, ma, 100.0)
    sweeps = Sweeps(max_steps; maxdim = D)
    observer = DMRGObserver(;energy_tol = tol)
    gs_energy, gs = dmrg(H, mps, sweeps; outputlevel = 1, observer = observer, ishermitian = true)
        
    rho = outer(gs', gs)
    println("The ground state energy was found to be $(gs_energy)\n")
    L_taylor_expanded_part_tmp = get_L_taylor(sites, x, l_0, ma, aD_0, sigma_over_a, env_corr_type, aT)
    L_taylor_expectation_value = real(inner(L_taylor_expanded_part_tmp, rho))
    println("The expectation value of the taylor expanded part of the Lindblad operator is $(L_taylor_expectation_value)\n")
    println("The L_taylor*dt/2 should be much less than one: $(L_taylor_expectation_value*tau/2)\n")
    println("The trace of initial rho should be 1: $(tr(rho))\n")
    flush(stdout)
    
    # Get the exponential of the Lindblad terms for evolution 
    odd_gates_2 = get_exp_Ho_list(sites, -1im*tau/2) # odd/2
    odd_gates = get_exp_Ho_list(sites, -1im*tau) # odd
    even_gates = get_exp_He_list(sites, -1im*tau) # even

    # Prepare the lists for the tracked observables and the associated MPO
    max_bond_list = Int[]
    ee_list = Float64[]
    step_num_list = Int[]
    z_mpo = [MPO(get_Z_site_operator(idx), sites) for idx in 1:N]
    z_list = [[] for _ in 1:N]

    # Push into the lists the initial state observables
    push!(step_num_list, 0)
    push!(max_bond_list, maxlinkdim(rho))
    for idx in 1:N
        push!(z_list[idx], real(tr(apply(rho, z_mpo[idx]))))
    end
    push!(ee_list, get_entanglement_entropy_mpo(rho, div(N, 2)+1:N, sites))
    
    # Print statements
    println("The time to get the initial rho and MPO lists is: $(time() - t)\n")
    println("Now starting the ATDDMRG algorithm\n")
    flush(stdout)

    t = time()
    for step in 1:max_steps
    
        if step == 1

            apply_odd!(odd_gates_2, rho; cutoff = cutoff)
            rho = rho/tr(rho)
        
            rho = apply_taylor_part(rho, cutoff, tau, sites, x, l_0, ma, aD_0, sigma_over_a, env_corr_type, aT)
            rho = rho/tr(rho)

            apply_even!(even_gates, rho; cutoff = cutoff)
            rho = rho/tr(rho)

            rho = apply_taylor_part(rho, cutoff, tau, sites, x, l_0, ma, aD_0, sigma_over_a, env_corr_type, aT)
            rho = rho/tr(rho)
            
        elseif step == max_steps
        
            apply_odd!(odd_gates, rho; cutoff = cutoff)
            rho = rho/tr(rho)

            rho = apply_taylor_part(rho, cutoff, tau, sites, x, l_0, ma, aD_0, sigma_over_a, env_corr_type, aT)
            rho = rho/tr(rho)

            apply_even!(even_gates, rho; cutoff = cutoff)
            rho = rho/tr(rho)

            rho = apply_taylor_part(rho, cutoff, tau, sites, x, l_0, ma, aD_0, sigma_over_a, env_corr_type, aT)
            rho = rho/tr(rho)

            apply_odd!(odd_gates_2, rho; cutoff = cutoff)
            rho = rho/tr(rho)
                
        else

            apply_odd!(odd_gates, rho; cutoff = cutoff)
            rho = rho/tr(rho)

            rho = apply_taylor_part(rho, cutoff, tau, sites, x, l_0, ma, aD_0, sigma_over_a, env_corr_type, aT)
            rho = rho/tr(rho)

            apply_even!(even_gates, rho; cutoff = cutoff)
            rho = rho/tr(rho)

            rho = apply_taylor_part(rho, cutoff, tau, sites, x, l_0, ma, aD_0, sigma_over_a, env_corr_type, aT)
            rho = rho/tr(rho)

        end

        # Take care of hermiticity and positivity of the density matrix
        rho = add(dag(swapprime(rho, 0, 1)), rho; cutoff = cutoff)/2 # fix hermiticity with rho -> rho dagger + rho over 2
        rho = apply(dag(swapprime(rho, 0, 1)), rho; cutoff = cutoff) # fix positivity with rho beta = rho beta over 2 dagger times rho beta over 2 - this results in a right canonical form MPO
        rho = rho/tr(rho)
        
        if (step % measure_every == 0) || (step == max_steps)

            # Measure the observables
            push!(step_num_list, step)
            push!(max_bond_list, maxlinkdim(rho))
            for idx in 1:N
                push!(z_list[idx], real(tr(apply(rho, z_mpo[idx]))))
            end
            push!(ee_list, get_entanglement_entropy_mpo(rho, div(N, 2)+1:N, sites))

            # Refresh time variable and print statement
            println("Step: $step, Time = $(time()-t), Average Step Time = $((time() - t)/measure_every)\n")            
            flush(stdout)
            t = time()

        end

    end

    # Prepare the initial state and the Z observable at the middle of the lattice to be tracked
    z_op = [project_zeroq(get_op(["Z"], [idx], N)) for idx in 1:N]
    rho = outer(gs', gs)
    rho_m = mpo_to_matrix(rho)
    rho_v = reshape(rho_m, length(rho_m))

    # Get the Lindblad operator and its exponential which is the evolution operator
    L = get_Lindblad_reduced_sparse_matrix(N, x, ma, l_0, lambda, aD_0, sigma_over_a, aT, env_corr_type)
    evolution_operator = exp(Matrix(L)*tau)
    
    # Prepare the list to store the tracked observables and get the initial state values
    ee_list_sparse = [get_entanglement_entropy_reduced_matrix(N, reshape(rho_v, rhodim, rhodim))]
    z_list_sparse = [[] for _ in 1:N]
    for idx in 1:N
        push!(z_list_sparse[idx], real(tr(rho_m*z_op[idx])))
    end
    
    # Do the evolution
    for _ in 1:max_steps

        # One time step evolution
        rho_v = evolution_operator*rho_v

        # Measure tracked observables
        push!(ee_list_sparse, get_entanglement_entropy_reduced_matrix(N, reshape(rho_v, rhodim, rhodim)))
        for idx in 1:N
            push!(z_list_sparse[idx], real(tr(reshape(rho_v, (rhodim, rhodim))*z_op[idx])))
        end
        
    end

    # p1 = plot(step_num_list, max_bond_list)
    # title!("Maximum Bond Dimension")
    # display(p1)
    
    p2 = plot(step_num_list, ee_list, label = "MPO")
    plot!(step_num_list, ee_list_sparse, label = "Sparse")
    title!("Entanglement Entropy")
    display(p2)

    p3 = plot(step_num_list, z_list[div(N, 2)], label = "MPO, $(div(N, 2))")
    plot!(step_num_list, z_list_sparse[div(N, 2)], label = "Sparse, $(div(N, 2))", linestyle = :dash)
    title!("Middle Z Operator Expectation Value")
    display(p3)

    # norm_diff = []
    # for (idx, element) in enumerate(max_imag)
    #     push!(norm_diff, norm(element - max_imag_sp[idx]))
    #     # println(idx, ": ", norm(element - max_imag_sp[idx]))
    # end

    # p4 = plot(norm_diff)
    # display(p4)

end

run_ATDDMRG()
