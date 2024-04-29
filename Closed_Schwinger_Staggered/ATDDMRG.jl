using ITensors
using LinearAlgebra
using HDF5
using Statistics
using Plots
include("Utilities.jl")

N = parse(Int, ARGS[1])
tau = parse(Float64, ARGS[2])
cutoff = parse(Float64, ARGS[3])
tol = parse(Float64, ARGS[4])
x = parse(Float64, ARGS[5])
l_0 = parse(Float64, ARGS[6])
ma = parse(Float64, ARGS[7])
max_steps = parse(Int, ARGS[8])
project_number = parse(Int, ARGS[9])
l_0_initial = parse(Bool, ARGS[10])
h5_path = ARGS[11]
measure_every = parse(Int, ARGS[12])
h5_previous_path = ARGS[13]
lambda = parse(Float64, ARGS[14])
D = parse(Int, ARGS[15])
max_sweeps = parse(Int, ARGS[16])
l_0_small = parse(Float64, ARGS[17])
omega = parse(Float64, ARGS[18])
type = ARGS[19]
file = h5open(h5_path, "w")

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

function get_dmrg_results(N, x, l_0_initial, ma, lambda, tol, max_sweeps, D, file)

    sites = siteinds("S=1/2", N, conserve_qns = true) 
    H = get_aH_Hamiltonian(sites, x, l_0_initial, ma, lambda)
    state = [isodd(n) ? "0" : "1" for n = 1:N]
    psi0 = randomMPS(sites, state, linkdims = 2)
    observer = DMRGObserver(;energy_tol = tol)
    sweeps = Sweeps(max_sweeps; maxdim = D)

    dmrg_ground_energy, dmrg_ground_state = dmrg(H, psi0, sweeps; outputlevel = 1, observer = observer, ishermitian = true) 

    write(file, "dmrg_ground_energy", dmrg_ground_energy)
    write(file, "dmrg_ground_state", dmrg_ground_state)

    return sites, dmrg_ground_state, H

end

function run_attDMRG(N, tau, cutoff, tol, x, l_0, l_0_small, type, omega, ma, max_steps, l_0_initial, measure_every, h5_previous_path, lambda, D, max_sweeps, file)

    t = time()

    energy_list = Float64[]
    max_bond_list = Int[]
    avg_bond_list = Float64[]
    entanglement_entropy_list = Float64[]
    z_config_list = [Float64[] for _ in 1:N]
    at_list = Float64[]
    
    # Prepare the initial state for evolution
    if h5_previous_path == "None"
        println("Initializing the MPS with the ground state of the hamiltonian\n")
        flush(stdout)
        sites, mps, H = get_dmrg_results(N, x, l_0_initial, ma, lambda, tol, max_sweeps, D, file)
    else
        println("Initializing the MPS from h5_previous_path = $(h5_previous_path)\n")
        flush(stdout)
        previous_file = h5open(h5_previous_path, "r")
        mps = read(previous_file, "last_mps", MPS)
        orthogonalize!(mps, 1) # put the MPS in right canonical form as it is saved in left canonical form
        sites = dag(reduce(vcat, siteinds(mps; :plev => 0)))
        H = get_aH_Hamiltonian(sites, x, l_0, ma, lambda)
    end

    Ho_mpo_list = get_exp_Ho(sites, -1im*tau/2) # odd/2
    He_mpo_list = get_exp_He(sites, -1im*tau) # even
    Ho_mpo_list_2 = get_exp_Ho(sites, -1im*tau) # odd
    Hz_mpo = get_exp_Hz(sites, -1im*tau/2, x, l_0, ma) # 1+aH_z/2

    at = 0.0
    push!(at_list, 0.0)
    push!(energy_list, real(inner(mps', H, mps)))
    push!(max_bond_list, maxlinkdim(mps))
    push!(avg_bond_list, mean(linkdims(mps)))
    push!(entanglement_entropy_list, get_entanglement_entropy(mps, div(N, 2)))
    z_config = real(get_Z_configuration(mps))
    for i in 1:N
        push!(z_config_list[i], z_config[i])
    end

    println("The time to get the Hamiltonian, initial MPS and MPO lists is: $(time() - t)\n")
    flush(stdout)
    
    println("Now starting the attDMRG algorithm\n")
    flush(stdout)

    t = time()
    for step in 1:max_steps

        at += tau
        if type != "static"
            applied_field = get_applied_field(at, l_0, l_0_small, type, omega)
            Hz_mpo = get_exp_Hz(sites, -1im*tau/2, x, applied_field, ma) # 1+aH_z/2
        end

        if step == 1
            apply_Ho_mpo_list!(Ho_mpo_list, mps; cutoff = cutoff)
            mps = apply(Hz_mpo, mps; cutoff = cutoff, normalize = true)
            apply_He_mpo_list!(He_mpo_list, mps; cutoff = cutoff)
            mps = apply(Hz_mpo, mps; cutoff = cutoff, normalize = true)
        elseif step == max_steps
            apply_Ho_mpo_list!(Ho_mpo_list_2, mps; cutoff = cutoff)
            mps = apply(Hz_mpo, mps; cutoff = cutoff, normalize = true)
            apply_He_mpo_list!(He_mpo_list, mps; cutoff = cutoff)
            mps = apply(Hz_mpo, mps; cutoff = cutoff, normalize = true)
            apply_Ho_mpo_list!(Ho_mpo_list, mps; cutoff = cutoff)    
        else
            apply_Ho_mpo_list!(Ho_mpo_list_2, mps; cutoff = cutoff)
            mps = apply(Hz_mpo, mps; cutoff = cutoff, normalize = true)
            apply_He_mpo_list!(He_mpo_list, mps; cutoff = cutoff)
            mps = apply(Hz_mpo, mps; cutoff = cutoff, normalize = true)
        end
        
        if step % measure_every == 0

            push!(at_list, at)
            push!(energy_list, real(inner(mps', H, mps)))
            push!(max_bond_list, maxlinkdim(mps))
            push!(avg_bond_list, mean(linkdims(mps)))
            push!(entanglement_entropy_list, get_entanglement_entropy(mps, div(N, 2)))
            z_config = real(get_Z_configuration(mps))
            for i in 1:N
                push!(z_config_list[i], z_config[i])
            end

            println("Step: $step, Time = $(time()-t)\n")
            flush(stdout)
            t = time()

        end
    end

    write(file, "at_list", at_list)
    write(file, "energy_list", energy_list)
    write(file, "max_bond_list", max_bond_list)
    write(file, "avg_bond_list", avg_bond_list)
    write(file, "entanglement_entropy_list", entanglement_entropy_list)
    for idx in 1:N
        write(file, "z_list_$(idx)", z_config_list[idx])
    end
    close(file)

end

run_attDMRG(N, tau, cutoff, tol, x, l_0_val, l_0_small, type, omega, ma, max_steps, l_0_initial, measure_every, h5_previous_path, lambda, D, max_sweeps, file)
close(file)
println("Finished.")
flush(stdout)
