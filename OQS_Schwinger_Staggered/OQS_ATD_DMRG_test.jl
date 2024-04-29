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

N = 60
tau = 0.001 # 1/N^2 # time step in time evolution rho -> exp(-tau L) after one step
cutoff = 1e-9 # cutoff for SVD
max_rho_D = 100
max_steps = 100
tol = 1e-9 # tolerance for DMRG convergence and ATDDMRG convergence
e = 1
x = 1/(e)^2
ma = 1
l_0 = 0.02
l_0_small = 0.5*l_0
type = "static"
omega = 1.91*ma
lambda = 0.0
aD_0 = 0.1
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

function get_applied_field(at, l_0, l_0_small, type, omega)

    if type == "sauter"
        return l_0 + l_0_small/cosh(omega*at)^2
    elseif type == "gaussian"
        return l_0 + l_0_small*exp(-(omega*at)^2)
    elseif type == "static"
        return l_0
    else
        return l_0 + l_0_small*cos(omega*at)
    end

end

function run_ATDDMRG()

    t = time()
    
    # Prepare initial rho
    # println("Initializing with the MPO corresponding to the ground state of H_system\n")
    sites = siteinds("S=1/2", N, conserve_qns = true)
    state = [isodd(n) ? "0" : "1" for n = 1:N]
    mps = randomMPS(sites, state)
    H = get_aH_Hamiltonian(sites, x, l_0_initial, ma, lambda)
    sweeps = Sweeps(max_steps; maxdim = 100)
    observer = DMRGObserver(;energy_tol = tol)
    gs_energy, gs = dmrg(H, mps, sweeps; outputlevel = 1, observer = observer, ishermitian = true)
    # println("The ground state energy was found to be $(gs_energy)\n")
    rho = outer(gs', gs; cutoff = cutoff)

    # Test Dirac vacuum as initial state
    # state = [isodd(n) ? "1" : "0" for n = 1:N]
    # mps = MPS(sites, state)
    # rho = outer(mps', mps; cutoff = cutoff)
    # rho = get_dirac_vacuum_density_matrix(sites)
    
    # Get the exponential of the Lindblad terms for evolution 
    odd_gates_2 = get_exp_Ho_list(sites, -1im*tau/2) # odd/2
    odd_gates = get_exp_Ho_list(sites, -1im*tau) # odd
    even_gates = get_exp_He_list(sites, -1im*tau) # even

    # Prepare the lists for the tracked observables and the associated MPO
    at_list = [0.0]
    particle_number_mpo = get_particle_number_MPO(sites)
    pnd = [real(tr(apply(rho, particle_number_mpo)))/N]
    rhodim = Int(binomial(N, div(N, 2)))

    # Print statements
    println("The time to get the initial rho and MPO lists is: $(time() - t)\n")
    println("Now starting the ATDDMRG algorithm\n")
    flush(stdout)

    t = time()
    at = 0
    for step in 1:max_steps

        println("Step $(step) begins")
        flush(stdout)

        at += tau
        applied_field = get_applied_field(at, l_0, l_0_small, type, omega)
        
        if step == 1

            t1 = time()
            apply_odd!(odd_gates_2, rho; cutoff = cutoff, max_rho_D = max_rho_D)
            println("In step $step the first odd gates took ", time() - t1)   
            t1 = time()     
            rho = apply_taylor_part(rho, tau, sites, x, applied_field, ma, aD_0, sigma_over_a, env_corr_type, aT; cutoff = cutoff, max_rho_D = max_rho_D)
            println("In step $step the first taylor gates took ", time() - t1) 
            t1 = time()
            apply_even!(even_gates, rho; cutoff = cutoff, max_rho_D = max_rho_D)
            println("In step $step the first even gates took ", time() - t1)   
            t1 = time()
            rho = apply_taylor_part(rho, tau, sites, x, applied_field, ma, aD_0, sigma_over_a, env_corr_type, aT; cutoff = cutoff, max_rho_D = max_rho_D)      
            println("In step $step the second taylor gates took ", time() - t1)       
        
        elseif step == max_steps
        
            t1 = time()
            apply_odd!(odd_gates, rho; cutoff = cutoff, max_rho_D = max_rho_D)
            println("In step $step the first odd gates took ", time() - t1)
            t1 = time()   
            rho = apply_taylor_part(rho, tau, sites, x, applied_field, ma, aD_0, sigma_over_a, env_corr_type, aT; cutoff = cutoff, max_rho_D = max_rho_D)
            println("In step $step the first taylor gates took ", time() - t1)
            t1 = time()   
            apply_even!(even_gates, rho; cutoff = cutoff, max_rho_D = max_rho_D)
            println("In step $step the first even gates took ", time() - t1)
            t1 = time()   
            rho = apply_taylor_part(rho, tau, sites, x, applied_field, ma, aD_0, sigma_over_a, env_corr_type, aT; cutoff = cutoff, max_rho_D = max_rho_D)
            println("In step $step the second taylor gates took ", time() - t1)
            t1 = time()   
            apply_odd!(odd_gates_2, rho; cutoff = cutoff, max_rho_D = max_rho_D)
            println("In step $step the second odd gates took ", time() - t1)   

                
        else

            t1 = time()
            apply_odd!(odd_gates, rho; cutoff = cutoff, max_rho_D = max_rho_D)
            println("In step $step the first odd gates took ", time() - t1)
            t1 = time()
            rho = apply_taylor_part(rho, tau, sites, x, applied_field, ma, aD_0, sigma_over_a, env_corr_type, aT; cutoff = cutoff, max_rho_D = max_rho_D)
            println("In step $step the first taylor gates took ", time() - t1)
            t1 = time()
            apply_even!(even_gates, rho; cutoff = cutoff, max_rho_D = max_rho_D)
            println("In step $step the first even gates took ", time() - t1)
            t1 = time()
            rho = apply_taylor_part(rho, tau, sites, x, applied_field, ma, aD_0, sigma_over_a, env_corr_type, aT; cutoff = cutoff, max_rho_D = max_rho_D)
            println("In step $step the second taylor gates took ", time() - t1)

        end

        # Take care of hermiticity and positivity of the density matrix
        rho = add(dag(swapprime(rho, 0, 1)), rho; cutoff = cutoff, maxdim = max_rho_D)/2 # fix hermiticity with rho -> rho dagger + rho over 2
        rho = rho/tr(rho) # fix the trace to 1 again
        
        if (step % measure_every == 0) || (step == max_steps)

            # Measure the observables
            push!(at_list, at_list[end] + tau*measure_every)
            push!(pnd, real(tr(apply(rho, particle_number_mpo)))/N)
    
            # Refresh time variable and print statement
            tmp_time_measurement = (time() - t)/measure_every
            println("Step $(step) finished, Time = $(time()-t), Average Step Time = $(tmp_time_measurement), Linkdims = $(linkdims(rho)), applied field = $(applied_field), SPND = $(pnd[end]-pnd[1])\n")            
            flush(stdout)
            t = time()

        end

    end
    
    spnd = pnd .- pnd[1]
    display(plot(at_list, spnd))

end

run_ATDDMRG()
