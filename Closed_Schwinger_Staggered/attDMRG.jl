using ITensors
using LinearAlgebra
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
h5_path = parse(String, ARGS[11])
measure_every = parse(Int, ARGS[12])
h5_previous_path = parse(String, ARGS[13])
sites = siteinds("S=1/2", N, conserve_qns = true) 
file = h5open(h5_path, "w")
write(file, "sites", sites)

function get_dmrg_results()

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
    H = get_Hamiltonian(sites, x, l_0, mg)
    
    if h5_previous_path[end] != "5" 
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
        max_mps_key = "mps_$(max_mps_key_num)"
        mps = read(previous_file, max_mps_key, MPS)
        sites = read(previous_file, "sites", Vector{Index{Int64}})
    end
    println("The initial MPS needs to start from right canonical form: $(get_which_canonical_form(mps))\n")
    flush(stdout)

    Ho_mpo_list = get_exp_Ho(sites, -tau/2, x) # odd/2
    Hz_mpo = get_exp_Hz(sites, -tau/2, x, l_0, mg) # 1+aH_z/2
    He_mpo_list = get_exp_He(sites, -tau, x, l_0) # even
    Ho_mpo_list_2 = get_exp_Ho(sites, -tau, x) # odd

    E_previous = real(inner(mps', H, mps))
    push!(energy_list, (0, E_previous)) # the 0 here is the step number
    push!(max_bond_list, (0, maxlinkdim(mps)))
    E_current = 0
    
    println("The time to get the Hamiltonian, initial MPS, MPO lists and the first E_previous is: $(time() - t)\n")
    flush(stdout)
    
    println("Now starting the iattDMRG algorithm\n")
    flush(stdout)

    println("Step: 0, E = $(E_previous)\n")
    flush(stdout)
    write(file, "mps_0", mps)

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

            E_current = real(inner(mps', H, mps))
            push!(energy_list, (step, E_current))
            push!(max_bond_list, (step, maxlinkdim(mps)))
            write(file, "mps_$(step)", mps)

            println("Step: $step, E = $E_current, Time = $(time()-t)")
            flush(stdout)
            t = time()
            
            e = abs(E_current-E_previous)/N 
            if e < tol
                println("The absolute value of the difference in energy at time step $step was found to be $e which is less than tol = $tol, hence the while loop breaks here.")
                flush(stdout)
                break
            end

            E_previous = E_current
            
        end

        if step == max_steps
            write(file, "mps_$(step)", mps)
        end

    end
    
    println("The absolute value of the difference in energy after $max_steps steps did not reach the desired tol = $tol, hence the function stops here.")
    flush(stdout)

    write(file, "energy_list", energy_list)
    write(file, "max_bond_list", max_bond_list)

end

function run_attDMRG()

    t = time()
    energy_list = []
    max_bond_list = []
    particle_number_list = []
    entanglement_entropy_list = []
    z_config_list = []
    q_config_list = []
    ef_config_list = []
    step_num = []
    
    H = get_Hamiltonian(sites, x, l_0, mg)
    PN = get_particle_number_MPO(sites)
    
    if h5_previous_path[end] != "5" 
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
        max_mps_key = "mps_$(max_mps_key_num)"
        mps = read(previous_file, max_mps_key, MPS)
        sites = read(previous_file, "sites", Vector{Index{Int64}})
    end
    println("The initial MPS needs to start from right canonical form: $(get_which_canonical_form(mps))\n")
    flush(stdout)

    Ho_mpo_list = get_exp_Ho(sites, -1im*tau/2, x) # odd/2
    Hz_mpo = get_exp_Hz(sites, -1im*tau/2, x, l_0, mg) # 1+aH_z/2
    He_mpo_list = get_exp_He(sites, -1im*tau, x, l_0) # even
    Ho_mpo_list_2 = get_exp_Ho(sites, -1im*tau, x) # odd

    push!(step_num, 0)
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
            push!(step_num, step)
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

        if step == max_steps
            push!(step_num, step)
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
    get_dmrg_results()
    run_iattDMRG()
else
    run_attDMRG() 
end
