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

    sites = siteinds("S=1/2", N) # , conserve_qns = true) 
    H = get_Hamiltonian(sites, x, l_0, mg)
    # state = [isodd(n) ? "0" : "1" for n = 1:N]
    # psi0 = randomMPS(sites, state, linkdims = 2)
    psi0 = randomMPS(sites, linkdims = 2)
    obs = DMRGObserver(;energy_tol = tol)
    nsweeps = max_steps

    dmrg_energy, dmrg_state = dmrg(H, psi0; nsweeps, cutoff, observer=obs, outputlevel=1)

    write(file, "dmrg_energy", dmrg_energy)
    write(file, "dmrg_state", dmrg_state)

end

function run_iatdDMRG()

    t = time()
    
    # Prepare initial rho
    if h5_previous_path == "None" 
        println("Initializing with the identity MPO\n")
        sites = siteinds("S=1/2", N, conserve_qns = true)
        rho = MPO(sites, "Id")
        rho = rho/tr(rho)
    else
        println("Initializing the MPO from h5_previous_path = $(h5_previous_path)\n")
        flush(stdout)
        previous_file = h5open(h5_previous_path, "r")
        keys_previous_file = keys(previous_file)
        rho_keys = filter(key -> occursin("rho_", key), keys_previous_file)
        max_rho_key_num = rho_keys[argmax(parse.(Int, [split(item, "_")[2] for item in rho_keys]))]
        max_rho_key = "$(max_rho_key_num)"
        rho = read(previous_file, max_rho_key, MPO)
        orthogonalize!(rho, 1) # put the MPO in right canonical form
        rho = rho/tr(rho)
        sites = dag(reduce(vcat, siteinds(rho; :plev => 0)))
    end

    # Get the Hamiltonian MPO
    H = get_Hamiltonian(sites, x, l_0, mg)

    # Get the exponential of the Hamiltonian terms for evolution
    odd_gates_4 = get_exp_Ho_list(sites, -tau/4, x) # odd/4
    exp_Hz_mpo = get_exp_Hz(sites, -tau/4, x, l_0, mg) # 1+aH_z/4
    even_gates_2 = get_exp_He_list(sites, -tau/2, x) # even/2
    odd_gates_2 = get_exp_Ho_list(sites, -tau/2, x) # odd/2

    energy_list = Float64[]
    max_bond_list = Int[]
    step_num_list = Int[]
    E_previous = inner(H, rho) 
    push!(step_num_list, 0)
    push!(energy_list, E_previous) # the 0 here is the step number
    push!(max_bond_list, maxlinkdim(rho))
    E_current = 0
    
    println("The time to get the Hamiltonian, initial rho, MPO lists and the first E_previous is: $(time() - t)\n")
    flush(stdout)
    
    println("Now starting the iattDMRG algorithm\n")
    flush(stdout)

    println("Step: 0, E = $(E_previous)\n")
    flush(stdout)
    write(file, "rho_0", rho)

    t = time()
    no_convergence = true # flag for some print statement after the for loop
    for step in 1:max_steps
    
        if step == 1

            apply_odd!(odd_gates_4, rho; cutoff = cutoff)
            rho = rho/tr(rho)
            rho = apply(apply(exp_Hz_mpo, rho; cutoff = cutoff, normalize = true), replaceprime(dag(exp_Hz_mpo'), 2 => 0); cutoff = cutoff, normalize = true)
            apply_even!(even_gates_2, rho; cutoff = cutoff)
            rho = rho/tr(rho)
            rho = apply(apply(exp_Hz_mpo, rho; cutoff = cutoff, normalize = true), replaceprime(dag(exp_Hz_mpo'), 2 => 0); cutoff = cutoff, normalize = true)
            
        elseif step == max_steps
        
            apply_odd!(odd_gates_2, rho; cutoff = cutoff)
            rho = rho/tr(rho)
            rho = apply(apply(exp_Hz_mpo, rho; cutoff = cutoff, normalize = true), replaceprime(dag(exp_Hz_mpo'), 2 => 0); cutoff = cutoff, normalize = true)
            apply_even!(even_gates_2, rho; cutoff = cutoff)
            rho = rho/tr(rho)
            rho = apply(apply(exp_Hz_mpo, rho; cutoff = cutoff, normalize = true), replaceprime(dag(exp_Hz_mpo'), 2 => 0); cutoff = cutoff, normalize = true)
            apply_odd!(odd_gates_4, rho; cutoff = cutoff)
            rho = rho/tr(rho)
                
        else

            apply_odd!(odd_gates_2, rho; cutoff = cutoff)
            rho = rho/tr(rho)
            rho = apply(apply(exp_Hz_mpo, rho; cutoff = cutoff, normalize = true), replaceprime(dag(exp_Hz_mpo'), 2 => 0); cutoff = cutoff, normalize = true)
            apply_even!(even_gates_2, rho; cutoff = cutoff)
            rho = rho/tr(rho)
            rho = apply(apply(exp_Hz_mpo, rho; cutoff = cutoff, normalize = true), replaceprime(dag(exp_Hz_mpo'), 2 => 0); cutoff = cutoff, normalize = true)

        end

        rho = (dag(swapprime(rho, 0, 1)) + rho)/2 # fix hermiticity
        
        if step % measure_every == 0

            rho = rho/tr(rho)
            E_current = inner(H, rho) 
            push!(step_num_list, step)
            push!(energy_list, E_current)
            push!(max_bond_list, maxlinkdim(rho))
            write(file, "rho_$(step)", rho)

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
                write(file, "rho_$(step)", rho)
            end
        end

    end
    
    if no_convergence
        println("The absolute value of the difference in energy after $max_steps steps did not reach the desired tol = $tol, hence the function stops here.\n")
        flush(stdout)
    end

end

if get_dmrg
    println("DMRG starting now\n")
    get_dmrg_results()
    println("DMRG is finished, starting iattDMRG now\n")
    run_iatdDMRG()
else
    run_iatdDMRG() 
end

close(file)
