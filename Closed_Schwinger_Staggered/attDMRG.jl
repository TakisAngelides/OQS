using ITensors
using LinearAlgebra
using HDF5
include("Utilities.jl")

N = parse(Int, ARGS[1])
tau = parse(Float64, ARGS[2])
cutoff = parse(Float64, ARGS[3])
tol = parse(Float64, ARGS[4])
x = parse(Float64, ARGS[5])
l_0 = parse(Float64, ARGS[6])
mg = parse(Float64, ARGS[7])
max_steps = parse(Int, ARGS[8])
project_number = parse(Int, ARGS[9])
get_dmrg = parse(Bool, ARGS[10])
h5_path = ARGS[11]
measure_every = parse(Int, ARGS[12])
h5_previous_path = ARGS[13]
file = h5open(h5_path, "w")

function get_dmrg_results()

    sites = siteinds("S=1/2", N, conserve_qns = true) 
    write(file, "sites_dmrg", sites)
    H = get_Hamiltonian(sites, x, l_0, mg)
    state = [isodd(n) ? "0" : "1" for n = 1:N]
    psi0 = randomMPS(sites, state, linkdims = 2)
    observer = DMRGObserver(;energy_tol = tol)
    sweeps = Sweeps(max_steps)

    dmrg_energy, dmrg_state = dmrg(H, psi0, sweeps; outputlevel = 1, observer = observer, cutoff = cutoff)

    write(file, "dmrg_energy", dmrg_energy)
    write(file, "dmrg_state", dmrg_state)

end

function run_iattDMRG()

    t = time()
    
    if h5_previous_path == "None" 
        sites = siteinds("S=1/2", N, conserve_qns = true) 
        println("Initializing with a random MPS\n")
        flush(stdout)
        state = [isodd(n) ? "0" : "1" for n = 1:N]
        mps = randomMPS(sites, state, linkdims = 2) # Starts in right canonical form
    else
        println("Initializing the MPS from h5_previous_path = $(h5_previous_path)\n")
        flush(stdout)
        previous_file = h5open(h5_previous_path, "r")
        keys_previous_file = keys(previous_file)
        mps_keys = filter(key -> occursin("mps_", key), keys_previous_file)
        max_mps_key_num = mps_keys[argmax(parse.(Int, [split(item, "_")[2] for item in mps_keys]))]
        max_mps_key = "$(max_mps_key_num)"
        mps = read(previous_file, max_mps_key, MPS)
        orthogonalize!(mps, 1) # put the MPS in right canonical form as it is saved in left canonical form
        sites = siteinds(mps)
    end
    H = get_Hamiltonian(sites, x, l_0, mg)

    Ho_mpo_list = get_exp_Ho(sites, -tau/2, x) # odd/2
    Hz_mpo = get_exp_Hz(sites, -tau/2, x, l_0, mg) # 1+aH_z/2
    He_mpo_list = get_exp_He(sites, -tau, x) # even
    Ho_mpo_list_2 = get_exp_Ho(sites, -tau, x) # odd

    energy_list = Float64[]
    max_bond_list = Int[]
    step_num_list = Int[]
    E_previous = real(inner(mps', H, mps))
    push!(step_num_list, 0)
    push!(energy_list, E_previous) # the 0 here is the step number
    push!(max_bond_list, maxlinkdim(mps))
    E_current = 0
    
    println("The time to get the Hamiltonian, initial MPS, MPO lists and the first E_previous is: $(time() - t)\n")
    flush(stdout)
    
    println("Now starting the iattDMRG algorithm\n")
    flush(stdout)

    println("Step: 0, E = $(E_previous)\n")
    flush(stdout)
    write(file, "mps_0", mps)

    t = time()
    no_convergence = true # flag for some print statement after the for loop
    for step in 1:max_steps
    
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

            E_current = real(inner(mps', H, mps))
            push!(step_num_list, step)
            push!(energy_list, E_current)
            push!(max_bond_list, maxlinkdim(mps))
            write(file, "mps_$(step)", mps)

            println("Step: $step, E = $E_current, Time = $(time()-t), Average Step Time = $((time() - t)/measure_every)\n")
            flush(stdout)
            t = time()
            
            e = abs(E_current-E_previous)/N 
            if e < tol
                no_convergence = false
                println("The absolute value of the difference in energy at time step $step was found to be $e which is less than tol = $tol, hence the while loop breaks here.\n")
                flush(stdout)
                write(file, "energy_list", energy_list)
                write(file, "max_bond_list", max_bond_list)
                write(file, "step_num_list", step_num_list)
                break
            end

            E_previous = E_current
            
        end

        if step == max_steps 
            write(file, "energy_list", energy_list)
            write(file, "max_bond_list", max_bond_list)
            write(file, "step_num_list", step_num_list)
            if step % measure_every != 0
                write(file, "mps_$(step)", mps)
            end
        end

    end
    
    if no_convergence
        println("The absolute value of the difference in energy after $max_steps steps did not reach the desired tol = $tol, hence the function stops here.\n")
        flush(stdout)
    end

end

function run_attDMRG()

    t = time()
    energy_list = Float64[]
    max_bond_list = Int[]
    particle_number_list = Float64[]
    entanglement_entropy_list = Float64[]
    z_config_list = Float64[]
    q_config_list = Float64[]
    ef_config_list = Float64[]
    step_num_list = Int[]
    
    H = get_Hamiltonian(sites, x, l_0, mg)
    PN = get_particle_number_MPO(sites)
    
    if h5_previous_path == "None"
        sites = siteinds("S=1/2", N, conserve_qns = true) 
        println("Initializing with a random MPS\n")
        flush(stdout)
        state = [isodd(n) ? "0" : "1" for n = 1:N]
        mps = randomMPS(sites, state, linkdims = 2) # Starts in right canonical form
    else
        println("Initializing the MPS from h5_previous_path = $(h5_previous_path)\n")
        flush(stdout)
        previous_file = h5open(h5_previous_path, "r")
        keys_previous_file = keys(previous_file)
        mps_keys = filter(key -> occursin("mps_", key), keys_previous_file)
        max_mps_key_num = argmax(mps_keys)
        max_mps_key = "$(max_mps_key_num)"
        mps = read(previous_file, max_mps_key, MPS)
        orthogonalize!(mps, 1) # put the MPS in right canonical form as it is saved in left canonical form
        sites = siteinds(mps)
    end

    Ho_mpo_list = get_exp_Ho(sites, -1im*tau/2, x) # odd/2
    Hz_mpo = get_exp_Hz(sites, -1im*tau/2, x, l_0, mg) # 1+aH_z/2
    He_mpo_list = get_exp_He(sites, -1im*tau, x) # even
    Ho_mpo_list_2 = get_exp_Ho(sites, -1im*tau, x) # odd

    push!(step_num_list, 0)
    push!(energy_list, real(inner(mps', H, mps)))
    push!(max_bond_list, maxlinkdim(mps))
    push!(particle_number_list, real(inner(mps', PN, mps)))
    push!(entanglement_entropy_list, get_entanglement_entropy(mps, div(N, 2)))
    push!(z_config_list, get_Z_configuration(mps))
    q, ef = get_charge_and_electric_field_configurations(mps, l_0)
    push!(q_config_list, q)
    push!(ef_config_list, ef)
    write(file, "mps_0", mps)
    
    println("The time to get the Hamiltonian, initial mps, MPO lists and the first E_previous is: $(time() - t)\n")
    flush(stdout)
    
    println("Now starting the attDMRG algorithm\n")
    flush(stdout)

    t = time()
    for step in 1:max_steps
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
            push!(step_num_list, step)
            push!(energy_list, real(inner(mps', H, mps)))
            push!(max_bond_list, maxlinkdim(mps))
            push!(particle_number_list, real(inner(mps', PN, mps)))
            push!(entanglement_entropy_list, get_entanglement_entropy(mps, div(N, 2)))
            push!(z_config_list, get_Z_configuration(mps))
            q, ef = get_charge_and_electric_field_configurations(mps, l_0)
            push!(q_config_list, q)
            push!(ef_config_list, ef)
            write(file, "mps_$(step)", mps)

            println("Step: $step, Time = $(time()-t)\n")
            flush(stdout)
            t = time()
        end

        if (step == max_steps) && (step % measure_every == 0)
            push!(step_num_list, step)
            push!(energy_list, real(inner(mps', H, mps)))
            push!(max_bond_list, maxlinkdim(mps))
            push!(particle_number_list, real(inner(mps', PN, mps)))
            push!(entanglement_entropy_list, get_entanglement_entropy(mps, div(N, 2)))
            push!(z_config_list, get_Z_configuration(mps))
            q, ef = get_charge_and_electric_field_configurations(mps, l_0)
            push!(q_config_list, q)
            push!(ef_config_list, ef)
            write(file, "mps_$(step)", mps)
        end

    end
end

if get_dmrg
    println("DMRG starting now\n")
    get_dmrg_results()
    println("DMRG is finished, starting iattDMRG now\n")
    run_iattDMRG()
else
    run_attDMRG() 
end

close(file)
