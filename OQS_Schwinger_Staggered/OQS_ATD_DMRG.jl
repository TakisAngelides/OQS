using ITensors
using LinearAlgebra
using HDF5
using SparseArrays
using Arpack
using KrylovKit
using TupleTools
using OpenQuantumTools
using Statistics
include("Utilities.jl")

# Inputs
N = parse(Int, ARGS[1])
tau = parse(Float64, ARGS[2]) # time step in time evolution rho -> exp(-tau L) after one step
cutoff = parse(Float64, ARGS[3]) # cutoff for SVD
tol = parse(Float64, ARGS[4]) # tolerance for DMRG convergence and ATDDMRG convergence
x = parse(Float64, ARGS[5])
l_0 = parse(Float64, ARGS[6])
ma = parse(Float64, ARGS[7])
max_steps = parse(Int, ARGS[8]) # another stopping condition for DMRG and ATDDMRG
project_number = parse(Int, ARGS[9])
h5_path = ARGS[10]
measure_every = parse(Int, ARGS[11]) # this determines how often to save rho and measure the energy in ATDDMRG
h5_previous_path = ARGS[12]
file = h5open(h5_path, "w")
D = parse(Int64, ARGS[13])
lambda = parse(Float64, ARGS[14])
aD_0 = parse(Float64, ARGS[15])
aT = parse(Float64, ARGS[16])
sigma_over_a = parse(Float64, ARGS[17])
env_corr_type = ARGS[18]
max_sweeps = parse(Int, ARGS[19])
sparse_evol = parse(Bool, ARGS[20])
l_0_initial_state = parse(Float64, ARGS[21])
dirac_vacuum_initial_state = parse(Bool, ARGS[22])
max_rho_D = parse(Int, ARGS[23])
l_0_small = parse(Float64, ARGS[24])
Omega = parse(Float64, ARGS[25])
omega = parse(Float64, ARGS[26])

function get_applied_field(at, l_0, l_0_small, Omega, omega)

    # we consider a strong and slow electric field pulse superimposed by a weak and fast pulse
    # E_S >> E >> epsilon > 0 and m >> omega >> Omega > 0 which implies l_0 >> l_0_small
    return l_0/cosh(Omega*at)^2 + l_0_small/cosh(omega*at)^2

end

function run_ATDDMRG()

    t = time()
    
    # Prepare initial rho
    if h5_previous_path == "None" 
        
        # Initialize rho from the ground state of the Hamiltonian or from the Dirac vacuum
        sites = siteinds("S=1/2", N, conserve_qns = true)
        
        if dirac_vacuum_initial_state
            println("Initializing with the MPO corresponding to the Dirac vacuum\n")
            rho_initial = get_dirac_vacuum_density_matrix(sites)
        else
            println("Initializing with the MPO corresponding to the ground state of H_system\n")
            state = [isodd(n) ? "0" : "1" for n = 1:N]
            mps = randomMPS(sites, state)
            H = get_aH_Hamiltonian(sites, x, l_0_initial_state, ma, lambda) 
            sweeps = Sweeps(max_sweeps; maxdim = D)
            observer = DMRGObserver(;energy_tol = tol)
            gs_energy, gs = dmrg(H, mps, sweeps; outputlevel = 1, observer = observer, ishermitian = true) 
            rho_initial = outer(gs', gs; cutoff = cutoff)
            println("The ground state energy was found to be $(gs_energy)\n")
        end
        rho = rho_initial
        # L_taylor_expanded_part_tmp = get_L_taylor(sites, x, l_0, ma, aD_0, sigma_over_a, env_corr_type, aT)
        # L_taylor_expectation_value = real(inner(L_taylor_expanded_part_tmp, rho))
        # println("The expectation value of the taylor expanded part of the Lindblad operator is $(L_taylor_expectation_value)\n")
        # println("The L_taylor*tau/2 should be much less than one: $(L_taylor_expectation_value*tau/2)\n")
        # println("The trace of initial rho should be 1: $(tr(rho))\n")
        flush(stdout)

    else 

        # This is the case when the time step has decreased but continues evolving a given state from a larger time step
        println("Initializing the MPO from h5_previous_path = $(h5_previous_path)\n")
        flush(stdout)
        previous_file = h5open(h5_previous_path, "r")
        keys_previous_file = keys(previous_file)
        rho_keys = filter(key -> occursin("rho_", key), keys_previous_file)
        max_rho_key_num = rho_keys[argmax(parse.(Int, [split(item, "_")[2] for item in rho_keys]))]
        max_rho_key = "$(max_rho_key_num)"
        rho = read(previous_file, max_rho_key, MPO)
        ITensors.orthogonalize!(rho, 1) # put the MPO in right canonical form
        rho = rho/tr(rho)
        sites = dag(reduce(vcat, siteinds(rho; :plev => 0)))
    
    end

    # Get the exponential of the Lindblad terms for evolution 
    odd_gates_2 = get_exp_Ho_list(sites, -1im*tau/2) # odd/2
    odd_gates = get_exp_Ho_list(sites, -1im*tau) # odd
    even_gates = get_exp_He_list(sites, -1im*tau) # even

    # Prepare the lists for the tracked observables and the associated MPO
    max_bond_list = Int[]
    avg_bond_list = Float64[]
    avg_step_time = Float64[]
    if N <= 8
        ee_list = Float64[]
    end
    at_list = Float64[]
    z_mpo = [MPO(get_Z_site_operator(idx), sites) for idx in 1:N]
    z_list = [Float64[] for _ in 1:N]

    # Push into the lists the initial state observables
    push!(at_list, 0.0)
    push!(max_bond_list, maxlinkdim(rho))
    push!(avg_bond_list, mean(linkdims(rho)))
    for idx in 1:N
        push!(z_list[idx], real(tr(apply(rho, z_mpo[idx]))))
    end
    if N <= 8
        push!(ee_list, get_entanglement_entropy_mpo(rho, div(N, 2)+1:N, sites))
    end

    # Write initial state to file and print statements
    println("The time to get the initial rho and MPO lists is: $(time() - t)\n")
    println("Now starting the ATDDMRG algorithm\n")
    flush(stdout)
    # write(file, "rho_0", rho)

    t = time()
    at = 0
    for step in 1:max_steps

        println("Step $(step) begins")
        flush(stdout)

        at += tau
        applied_field = get_applied_field(at, l_0, l_0_small, Omega, omega)
        
        if step == 1

            apply_odd!(odd_gates_2, rho; cutoff = cutoff, max_rho_D = max_rho_D)        
            rho = apply_taylor_part(rho, tau, sites, x, applied_field, ma, aD_0, sigma_over_a, env_corr_type, aT; cutoff = cutoff, max_rho_D = max_rho_D)
            apply_even!(even_gates, rho; cutoff = cutoff, max_rho_D = max_rho_D)
            rho = apply_taylor_part(rho, tau, sites, x, applied_field, ma, aD_0, sigma_over_a, env_corr_type, aT; cutoff = cutoff, max_rho_D = max_rho_D)            
        
        elseif step == max_steps
        
            apply_odd!(odd_gates, rho; cutoff = cutoff, max_rho_D = max_rho_D)
            rho = apply_taylor_part(rho, tau, sites, x, applied_field, ma, aD_0, sigma_over_a, env_corr_type, aT; cutoff = cutoff, max_rho_D = max_rho_D)
            apply_even!(even_gates, rho; cutoff = cutoff, max_rho_D = max_rho_D)
            rho = apply_taylor_part(rho, tau, sites, x, applied_field, ma, aD_0, sigma_over_a, env_corr_type, aT; cutoff = cutoff, max_rho_D = max_rho_D)
            apply_odd!(odd_gates_2, rho; cutoff = cutoff, max_rho_D = max_rho_D)
                
        else

            apply_odd!(odd_gates, rho; cutoff = cutoff, max_rho_D = max_rho_D)
            rho = apply_taylor_part(rho, tau, sites, x, applied_field, ma, aD_0, sigma_over_a, env_corr_type, aT; cutoff = cutoff, max_rho_D = max_rho_D)
            apply_even!(even_gates, rho; cutoff = cutoff, max_rho_D = max_rho_D)
            rho = apply_taylor_part(rho, tau, sites, x, applied_field, ma, aD_0, sigma_over_a, env_corr_type, aT; cutoff = cutoff, max_rho_D = max_rho_D)

        end

        # Take care of hermiticity and positivity of the density matrix
        rho = add(dag(swapprime(rho, 0, 1)), rho; cutoff = cutoff, maxdim = max_rho_D)/2 # fix hermiticity with rho -> rho dagger + rho over 2
        rho = rho/tr(rho) # fix the trace to 1 again
        
        if (step % measure_every == 0) || (step == max_steps)

            # Measure the observables
            push!(at_list, at_list[end] + tau*measure_every)
            push!(max_bond_list, maxlinkdim(rho))
            push!(avg_bond_list, mean(linkdims(rho)))
            for idx in 1:N
                push!(z_list[idx], real(tr(apply(rho, z_mpo[idx]))))
            end
            if N <= 8
                push!(ee_list, get_entanglement_entropy_mpo(rho, div(N, 2)+1:N, sites))
            end

            # Write the state to file
            # write(file, "rho_$(step)", rho)

            # Refresh time variable and print statement
            tmp_time_measurement = (time() - t)/measure_every
            push!(avg_step_time, tmp_time_measurement)
            println("Step $(step) finished, Time = $(time()-t), Average Step Time = $(tmp_time_measurement), Linkdims = $(linkdims(rho)), applied field = $(applied field)\n")            
            flush(stdout)
            t = time()

        end

    end

    # Write observable lists to file and the last density matrix
    write(file, "at_list", at_list)
    write(file, "max_bond_list", max_bond_list)
    write(file, "avg_bond_list", avg_bond_list)
    write(file, "avg_step_time", avg_step_time)
    for idx in 1:N
        write(file, "z_list_$(idx)", z_list[idx])
    end
    if N <= 8
        write(file, "ee_list", ee_list)
    end
    write(file, "rho", rho)

    # Sparse matrix evolution if required
    if sparse_evol

        # Print statement
        println("Now starting sparse evolution\n")
        flush(stdout)

        # Prepare initial rho
        if h5_previous_path == "None" 
            
            # Initialize rho from the ground state of the Hamiltonian
            println("Initializing with the MPO corresponding to the ground state of H_system for the sparse evolution\n")
            flush(stdout)
            rho_sparse_inital = rho_initial
            rho_matrix = mpo_to_matrix(rho_sparse_inital)
            rho_zeroq_subspace_matrix = project_zeroq(rho_matrix)
            rho_v = reshape(rho_zeroq_subspace_matrix, length(rho_zeroq_subspace_matrix))
            
        else 

            # This is the case when the time step has decreased but continues evolving a given state from a larger time step
            println("Initializing the sparse initial rho_v from h5_previous_path = $(h5_previous_path)\n")
            flush(stdout)
            previous_file = h5open(h5_previous_path, "r")
            rho_v = read(previous_file, "last_sparse_v")
            rho_zeroq_subspace_matrix = reshape(rho_v, rhodim, rhodim)
        
        end

        # Prepare the initial state and the Z observable at the middle of the lattice to be tracked
        z_op = [project_zeroq(get_op(["Z"], [idx], N)) for idx in 1:N]
        
        # Prepare the list to store the tracked observables and get the initial state values
        rhodim = Int(binomial(N, div(N, 2)))
        if N <= 8
            ee_list_sparse = [get_entanglement_entropy_reduced_matrix(N, rho_zeroq_subspace_matrix)]
        end
        z_list_sparse = [Float64[] for _ in 1:N]
        for idx in 1:N
            push!(z_list_sparse[idx], real(tr(rho_zeroq_subspace_matrix*z_op[idx])))
        end
        
        # Do the evolution
        at = 0
        for step in 1:max_steps

            # Get the Lindblad operator and its exponential which is the evolution operator
            at += tau
            applied_field = get_applied_field(at, l_0, l_0_small, Omega, omega)
            L = get_Lindblad_reduced_sparse_matrix(N, x, ma, applied_field, lambda, aD_0, sigma_over_a, aT, env_corr_type)
            evolution_operator = exp(Matrix(L)*tau)

            # One time step evolution
            rho_v = evolution_operator*rho_v

            # Measure tracked observables
            rho_zeroq_subspace_matrix = reshape(rho_v, rhodim, rhodim)
            if N <= 8
                push!(ee_list_sparse, get_entanglement_entropy_reduced_matrix(N, rho_zeroq_subspace_matrix))
            end
            for idx in 1:N
                push!(z_list_sparse[idx], real(tr(rho_zeroq_subspace_matrix*z_op[idx])))
            end
            
        end

        # Write the sparse observable lists to file
        if N <= 8
            write(file, "ee_list_sparse", ee_list_sparse)
        end
        for idx in 1:N
            write(file, "z_list_sparse_$(idx)", z_list_sparse[idx])
        end
        write(file, "last_sparse_v", rho_v)

    end

end

run_ATDDMRG()
close(file)

println("Finished")
flush(stdout)
