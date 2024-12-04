# Starting the global clock
t_total_initial = time()
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

# Perform the time evolution
function evolve()

    # Input arguments for file and opening the results h5 file
    println("Now getting the input arguments for the file and opening the results h5 file ", now())
    flush(stdout)
    path_to_project_number = ARGS[1] 
    path_to_inputs_h5 = "$(path_to_project_number)/inputs.h5" # Where the h5 file containing the inputs is
    job_id = ARGS[2] # The job id to get the inputs within this h5 file for this specific job
    path_to_results = "$(path_to_project_number)/HDF5/$(job_id).h5" # The path to the h5 file to save all results of states
    path_to_backup = "$(path_to_project_number)/HDF5_b/$(job_id)_b.h5"
    h5open(path_to_results, "w") do results_file # cw is read/write create if needed mode
    rerun_mode = isfile(path_to_backup)
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
    which_initial_state = inputs["wis"]
    conserve_qns = parse(Bool, inputs["cqns"])
    N = inputs["N"]
    function get_initial_state()

        sites_initial_state = siteinds("S=1/2", N, conserve_qns = conserve_qns)

        if rerun_mode # case of rerunning the script from saved data

            h5open(path_to_backup, "r") do f
                return read(f, "saved_mps", MPS)
            end

        elseif which_initial_state == "dirac_vacuum"

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

        elseif which_initial_state == "first_naive"

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

            state[1] = "1"
            state[2] = "0"
            psi = MPS(sites_initial_state, state)
            Ms = [gs]
            first_energy, first = dmrg(H, [gs], psi, sweeps, weight = 10, ishermitian = true, observer = observer, outputlevel = 1)

            write(results_file, "first_energy", first_energy)
            write(results_file, "first", first)
            rho = outer(first', first) # Get the density matrix
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
    mps = get_initial_state()
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
    sites = siteinds(mps)
    l_0_1 = inputs["l_0_1"]
    side = "left"
    H = get_double_aH_Hamiltonian(sites, x, l_0_1, ma, lambda, side)
    H = MPO(H, sites)
    H_kin, H_el, H_m = get_double_aH_Hamiltonian_individual_terms(N, x, l_0_1, side)
    H_kin, H_el, H_m = MPO(H_kin, sites), MPO(H_el, sites), MPO(H_m, sites)
    for i in 2:2:length(sites) # This is done so that the odd, even gates and taylor MPO have physical legs matching the purified MPS and combining this with the swapprime done on the operators later the transpose is taken on the operators acting on the even sites which correspond to operators acting on the right of the density matrix
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
    already_done_t_over_a = 0 # starting the time variable
    l_0 = get_applied_field(which_applied_field, inputs, already_done_t_over_a)
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

    if rerun_mode # case when the script ran at least once and saved a state but did not finish the whole evolution

        println("Preparation for the case of re-running the script ", now())
        flush(stdout)

        h5open(path_to_backup, "r") do f

            taylor_mpo = read(f, "taylor_mpo", MPO)
            z_configs = read(f, "z_configs")
            link_dims = read(f, "link_dims")
            energy = read(f, "energy")
            kin_energy = read(f, "kin_energy")
            m_energy = read(f, "m_energy")
            el_energy = read(f, "el_energy")
            H = read(f, "H", MPO)
            H_kin = read(f, "H_kin", MPO)
            H_m = read(f, "H_m", MPO)
            H_el = read(f, "H_el", MPO)
            already_done_steps = read(f, "ads")
            x = inputs["x"]
            ma = inputs["ma"]
            lambda = inputs["lambda"]
            aT = inputs["aT"]
            aD = inputs["aD"]
            env_corr_type = inputs["env_corr_type"]
            sites = siteinds(mps)
            l_0_1 = inputs["l_0_1"]
            side = "left"
            tau = inputs["tau"]
            dissipator_sites = inputs["ds"]
            which_applied_field = inputs["waf"]
            already_done_t_over_a = 0 # starting the time variable
            l_0 = get_applied_field(which_applied_field, inputs, already_done_t_over_a)
            number_of_time_steps = inputs["nots"]
            taylor_expansion_order = inputs["teo"]
            taylor_expansion_cutoff_1 = inputs["tec_1"]
            taylor_expansion_cutoff_2 = inputs["tec_2"]
            time_varying_applied_field_flag = parse(Bool, inputs["tvaff"])
            if time_varying_applied_field_flag
                already_done_t_over_a = read(f, "adtoa")
            end
            
        end

        println("The taylor_mpo with taylor order $(taylor_expansion_order) and cutoffs $(taylor_expansion_cutoff_1), $(taylor_expansion_cutoff_2) has bond dimensions ", linkdims(taylor_mpo))
        println("Finished preparation for the case of re-running the script ", now())
        flush(stdout)

    else # case where the script runs for the first time

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
        z_configs = zeros(ComplexF64, number_of_time_steps+1, N)
        z_configs[1, :] = measure_z_config(mps)
        link_dims = zeros(Int64, number_of_time_steps+1, 2*N-1)
        link_dims[1, :] = linkdims(mps)
        energy = zeros(ComplexF64, number_of_time_steps+1)
        energy[1] = measure_mpo(mps, H)
        kin_energy = zeros(ComplexF64, number_of_time_steps+1)
        kin_energy[1] = measure_mpo(mps, H_kin; alg = "naive")
        m_energy = zeros(ComplexF64, number_of_time_steps+1)
        m_energy[1] = measure_mpo(mps, H_m; alg = "naive")
        el_energy = zeros(ComplexF64, number_of_time_steps+1)
        el_energy[1] = measure_mpo(mps, H_el; alg = "naive")
        already_done_steps = 0
        time_varying_applied_field_flag = parse(Bool, inputs["tvaff"])
        println("Finished getting the lists for the tracked observables ", now())
        flush(stdout)

    end

    cutoff = inputs["cutoff"]
    maxdim = inputs["md"]
    which_steps_to_save_state = inputs["wstss"]
    mem = inputs["mem"]
    day_to_make_backup = inputs["dtmb"]*86400
    back_up_flag = true
    steps_to_backup = inputs["stb"]

    if time_varying_applied_field_flag

        t_over_a = 0 + already_done_t_over_a
        l_0_list = [get_applied_field(which_applied_field, inputs, t_over_a)]

        for step in 1+already_done_steps:number_of_time_steps

            t = time() # Starting the time for the step

            if ((t - t_total_initial) > day_to_make_backup && back_up_flag) || (step in steps_to_backup) # day_to_make_backup was converted upstairs to seconds and t-time() returns seconds

                back_up_flag = false

                h5open(path_to_backup, "w") do f

                    write(f, "taylor_mpo", taylor_mpo)
                    write(f, "z_configs", z_configs)
                    write(f, "link_dims", link_dims)
                    write(f, "energy", energy)
                    write(f, "kin_energy", kin_energy)
                    write(f, "m_energy", m_energy)
                    write(f, "el_energy", el_energy)
                    write(f, "H", H)
                    write(f, "H_kin", H_kin)
                    write(f, "H_m", H_m)
                    write(f, "H_el", H_el)
                    write(f, "ads", step-1)
                    write(f, "adtoa", already_done_t_over_a)
                    write(f, "l_0_list", l_0_list)
                    write(f, "saved_mps", mps)

                end

            end
    
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

            # Fix trace
            mps /= trace_mps(mps)

            # Compute the tracked observables
            z_configs[step+1, :] = measure_z_config(mps)
            linkdims_of_step = linkdims(mps)
            link_dims[step+1, :] = linkdims_of_step

            # Save state to file
            if step in which_steps_to_save_state
                write(results_file, "$(step)", mps)
            end

            total_mem_tmp = (Base.gc_live_bytes()/2^20)/10^3
            # if total_mem_tmp >= mem
            #     println("GC Cleaning")
            #     flush(stdout)
            #     GC.gc(true)
            # end
            println("Step = $(step), Time = $(time() - t), Links = $(linkdims(mps)), Mem = $(total_mem_tmp)")
            flush(stdout)

        end

        # Write tracked observables to results h5 file
        println("Now writing the observables to results h5 file ", now())
        flush(stdout)
        write(results_file, "z_configs", z_configs)
        write(results_file, "link_dims", link_dims)
        write(results_file, "l_0_list", l_0_list)
        println("Finished writing the observables to results h5 file ", now())
        flush(stdout)

    else

        for step in 1+already_done_steps:number_of_time_steps

            t = time() # Starting the time for the step

            if ((t - t_total_initial) > day_to_make_backup && back_up_flag) || (step in steps_to_backup) # day_to_make_backup was converted upstairs to seconds and t-time() returns seconds

                back_up_flag = false

                println("Now making back up to file ", now())

                h5open(path_to_backup, "w") do f

                    write(f, "taylor_mpo", taylor_mpo)
                    write(f, "z_configs", z_configs)
                    write(f, "link_dims", link_dims)
                    write(f, "energy", energy)
                    write(f, "kin_energy", kin_energy)
                    write(f, "m_energy", m_energy)
                    write(f, "el_energy", el_energy)
                    write(f, "H", H)
                    write(f, "H_kin", H_kin)
                    write(f, "H_m", H_m)
                    write(f, "H_el", H_el)
                    write(f, "ads", step-1)
                    write(f, "saved_mps", mps)

                end

                println("Finished making back up to file ", now())


            end

            # One time step with ATDDMRG
            apply_odd!(odd, mps, cutoff, maxdim)
            mps = apply(taylor_mpo, mps; cutoff = cutoff, maxdim = maxdim)
            apply_even!(even, mps, cutoff, maxdim)
            mps = apply(taylor_mpo, mps; cutoff = cutoff, maxdim = maxdim)
            apply_odd!(odd, mps, cutoff, maxdim)

            # Fix trace
            mps /= trace_mps(mps)

            # Compute the tracked observables
            z_configs[step+1, :] = measure_z_config(mps)
            linkdims_of_step = linkdims(mps)
            link_dims[step+1, :] = linkdims_of_step
            energy[step+1] = measure_mpo(mps, H)
            kin_energy[step+1] = measure_mpo(mps, H_kin)
            m_energy[step+1] = measure_mpo(mps, H_m)
            el_energy[step+1] = measure_mpo(mps, H_el)

            # Save state to file
            if step in which_steps_to_save_state
                write(results_file, "$(step)", mps)
            end

            total_mem_tmp = (Base.gc_live_bytes()/2^20)/10^3
            # if total_mem_tmp >= 0.5*mem
            #     println("GC Cleaning")
            #     flush(stdout)
            #     GC.gc(true)
            # end
            println("Step = $(step), Time = $(time() - t), Links = $(linkdims(mps)), Mem = $(total_mem_tmp)")
            flush(stdout)

        end

        # Write tracked observables to results h5 file
        println("Now writing the observables to results h5 file ", now())
        flush(stdout)
        write(results_file, "z_configs", z_configs)
        write(results_file, "link_dims", link_dims)
        write(results_file, "energy", energy)
        write(results_file, "kin_energy", kin_energy)
        write(results_file, "m_energy", m_energy)
        write(results_file, "el_energy", el_energy)
        println("Finished writing the observables to results h5 file ", now())
        flush(stdout)

    end

    end # ending the do block for the results file

end

evolve()

println("Finished ", now())
flush(stdout)
