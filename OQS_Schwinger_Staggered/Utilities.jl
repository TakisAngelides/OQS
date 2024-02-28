function get_exp_L_taylor(sites, a, x, l_0, ma, aD_0, sigma_over_a, env_corr_type, aT)

    # a here is the coefficient namely exp(a * L_taylor)
    mpo = a*get_L_taylor(sites, x, l_0, ma, aD_0, sigma_over_a, env_corr_type, aT)

    # exp(a * L_taylor) -> I + a * L_taylor
    final_mpo = MPO(sites, "Id") + mpo 

    return final_mpo

end

function get_L_taylor(sites, x, l_0, ma, aD_0, sigma_over_a, env_corr_type, aT)

    # Gets the part of the Lindblad operator to be taylor expanded in the trotterization

    N = length(sites)

    opsum = OpSum()

    for n in 1:N-1
        
        for k in n+1:N

            opsum += (0.25/x)*(N-k),"Z",n,"Z",k

        end

        opsum += ((N/8 - 0.25*ceil((n-1)/2) + l_0*(N-n)/2)/x),"Z",n

        opsum += (0.5*ma*(-1)^(n-1)),"Z",n

    end

    opsum += (0.5*ma*(-1)^(N-1)),"Z",N

    opsum += ((l_0^2)*(N-1)/(2*x) + (l_0*N/(4*x)) + (N^2)/16),"Id",1

    mpo1 = MPO(opsum, sites)

    mpo2 = get_Lindblad_dissipative_part(aD_0, sigma_over_a, env_corr_type, aT, sites)

    return mpo1 + mpo2

end

function get_H_taylor(sites, x, l_0, ma)

    # Gets the part of the Lindblad operator to be taylor expanded in the trotterization

    N = length(sites)

    opsum = OpSum()

    for n in 1:N-1
        
        for m in n+1:N

            opsum += (0.25/x)*(N - m),"Z",n,"Z",m

        end

        opsum += ((N/8 - 0.25*ceil((n-1)/2) + l_0*(N-n)/2)/x),"Z",n

        opsum += (0.5*ma*(-1)^(n-1)),"Z",n

    end

    opsum += (0.5*ma*(-1)^(N-1)),"Z",N

    opsum += ((l_0^2)*(N-1)/(2*x) + (l_0*N/(4*x)) + (N^2)/(16*x)),"Id",1

    return MPO(opsum, sites)

end

function get_exp_Ho_list(sites, a)::Vector{ITensor}

    """

    a = prefactor of H_o eg: -i tau / 2 to give exp(-i * tau * Ho / 2)

    This list of operators incldues the odd terms of the kinetic term

    Note: (XX+YY)/2 = S+S- + S-S+

    """

    gates = []
    N = length(sites)

    for n=1:2:(N-1)

        hj = 0.5 * op("S-", sites[n]) * op("S+", sites[n+1])
        hj += 0.5 * op("S+", sites[n]) * op("S-", sites[n+1])
        Gj = exp(a * hj)
        push!(gates, Gj)

    end

    return gates

end

function get_exp_He_list(sites, a)::Vector{ITensor}

    """

    a = prefactor of H_e eg: -i tau to give exp(-i * tau * He)

    This list of operators incldues the even terms of the kinetic term and we also includes the identity operator of H

    Note: (XX+YY)/2 = S+S- + S-S+

    """

    gates = []
    N = length(sites)

    for n=2:2:(N-2)

        hj = 0.5 * op("S-", sites[n]) * op("S+", sites[n+1])
        hj += 0.5 * op("S+", sites[n]) * op("S-", sites[n+1])
        Gj = exp(a * hj)
        push!(gates, Gj)

    end

    return gates

end

function get_aH_Hamiltonian(sites, x, l_0, ma, lambda)

    """
    This gives aH Hamiltonian
    """

    N = length(sites)

    opsum = OpSum()

    for n in 1:N-1
        
        for m in n+1:N
            
            # Long range ZZ interaction term
            opsum += 0.25*(1/x)*(N-m+lambda),"Z",n,"Z",m

        end

        # Kinetic term
        opsum += 0.5,"S+",n,"S-",n+1
        opsum += 0.5,"S-",n,"S+",n+1

        opsum += (1/x)*(N/8 - 0.25*ceil((n-1)/2) + l_0*(N-n)/2),"Z",n
        
        opsum += (0.5*ma*(-1)^(n-1)),"Z",n

    end

    opsum += (0.5*ma*(-1)^(N-1)),"Z",N

    opsum += ((l_0^2)*(N-1)/(2*x) + (l_0*N)/(4*x) + (N^2)/(16*x) + (lambda*N/(8*x))),"Id",1

    return MPO(opsum, sites)

end

function get_which_canonical_form(mps)

    N = length(mps)
    canonical_form::Array{String} = []

    for site in 1:N

        mps_site = mps[site]
        
        a = mps_site
        adag = dag(mps_site)
        if site != 1
            adag_idx = commonind(a, mps[site-1])
            replaceind!(adag, adag_idx, prime(adag_idx))
        end
        res = a*adag
        inds_res = inds(res)
        res = ITensors.Array(res, inds_res...)
        s = size(res)
        if length(s) == 0
            l = 1
        else
            l = s[1]
        end
        is_right = isapprox(res, I(l))

        mps_site = mps[site]
        a = mps_site
        adag = dag(mps_site)
        if site != N
            adag_idx = commonind(a, mps[site+1])
            replaceind!(adag, adag_idx, prime(adag_idx))
        end
        res = a*adag
        inds_res = inds(res)
        res = ITensors.Array(res, inds_res...)
        s = size(res)
        if length(s) == 0
            l = 1
        else
            l = s[1]
        end
        is_left = isapprox(res, I(l))

        if is_left
            if is_right 
                push!(canonical_form, "L/R")
            else
                push!(canonical_form, "L")
            end
        elseif is_right
            push!(canonical_form, "R")
        else
            push!(canonical_form, "N")
        end

    end

    return canonical_form

end

function apply_odd!(Ho_mpo_list, mpo; cutoff = 1e-20)

    N = length(mpo)

    for (idx, gate) in enumerate(Ho_mpo_list)

        # idx will be 1, 2, 3, ... and we set it below to the left index of the gate to be applied 1, 3, 5, ... 
        idx = 2*idx-1

        # println("Odd gate left index is $(idx) and the MPO canonical form is ", get_MPO_canonical_form(mpo))

        tmp = replaceprime(prime(gate'*mpo[idx]*mpo[idx+1]; :tags => "Site")*hermitian_conjugate_mpo(gate), 3 => 1)

        U, S, V = ITensors.svd(tmp, commoninds(tmp, mpo[idx])..., lefttags = "Link,l=$(idx)", righttags = "Link,l=$(idx)", cutoff = cutoff)
        
        mpo[idx] = U

        # S = S/norm(S)
        
        mpo[idx+1] = S*V

        # Extra SVD for ATTDMRG compared to TEBD
        if idx != N-1
        
            idx += 1
        
            tmp = mpo[idx]*mpo[idx+1]
        
            U,S,V = ITensors.svd(tmp, commoninds(tmp, mpo[idx])..., lefttags = "Link,l=$(idx)", righttags = "Link,l=$(idx)", cutoff = cutoff)
        
            mpo[idx] = U

            # S = S/norm(S)
        
            mpo[idx+1] = S*V
        
        end

    end

end

function apply_even!(He_mpo_list, mpo; cutoff = 1e-20)

    """
    After we apply 1-tau*Hz/2 with the apply function we end up with right canonical form.
    To apply the first even gate which starts at the second site we need to have the first site in Left canonical form
    """

    # mps[1], mps[2] = ITensors.qr(mps[1]*mps[2], uniqueinds(mps[1], mps[2]); positive = true, tags = "Link,l=$(1)")
    mpo[1], S, V = ITensors.svd(mpo[1]*mpo[2], commoninds(mpo[1], mpo[1]*mpo[2])..., lefttags = "Link,l=$(1)", righttags = "Link,l=$(1)")
    # S = S/norm(S)
    mpo[2] = S*V

    for (idx, gate) in enumerate(He_mpo_list)

        idx = 2*idx

        # println("Even gate left index is $(idx) and the MPO canonical form is ", get_MPO_canonical_form(mpo))

        tmp = replaceprime(prime(gate'*mpo[idx]*mpo[idx+1]; :tags => "Site")*hermitian_conjugate_mpo(gate), 3 => 1)
        
        U, S, V = ITensors.svd(tmp, commoninds(tmp, mpo[idx])..., lefttags = "Link,l=$(idx)", righttags = "Link,l=$(idx)", cutoff = cutoff)
        
        mpo[idx] = U

        # S = S/norm(S)
        
        mpo[idx+1] = S*V

        # Extra SVD for ATTDMRG compared to TEBD
        
        idx += 1
    
        tmp = mpo[idx]*mpo[idx+1]
    
        U,S,V = ITensors.svd(tmp, commoninds(tmp, mpo[idx])..., lefttags = "Link,l=$(idx)", righttags = "Link,l=$(idx)", cutoff = cutoff)
    
        mpo[idx] = U

        # S = S/norm(S)
    
        mpo[idx+1] = S*V
    
    end

end

function get_entanglement_entropy(psi, site, tol = 1e-12)
    
    ITensors.orthogonalize!(psi, site)
    if site == 1
        U,S,V = svd(psi[site], siteind(psi, site); cutoff = 1e-20)
    else
        U,S,V = svd(psi[site], (linkind(psi, site-1), siteind(psi, site)); cutoff = 1e-20)
    end

    SvN = sum(-real(singular_value^2)*log(real(singular_value^2)) for singular_value in diag(S) if real(singular_value^2) >= tol)

    return SvN

end

function get_Z_site_operator(site)

    ampo::Sum{Scaled{ComplexF64, Prod{Op}}} = OpSum()

    ampo += "Z",site

    return ampo

end

function get_Z_configuration(psi)

    n = length(psi)

    sites = siteinds(psi)

    res = []

    for i in 1:n
    
        push!(res, inner(psi', get_MPO_from_operator_sum(get_Z_site_operator(i), sites), psi))

    end

    return res

end

function get_charge_and_electric_field_configurations(psi, l_0)

    Z_configuration = get_Z_configuration(psi)

    N = div(length(psi), 2)

    n_links = N - 1

    Q_configuration = []

    for n in 1:N

        Q_i = 0.5*(Z_configuration[n] + (-1)^(n-1))

        push!(Q_configuration, Q_i)

    end

    E_configuration = []

    for n in 1:n_links

        E_i = l_0 + sum(Q_configuration[1:n])

        push!(E_configuration, E_i)

    end

    return Q_configuration, E_configuration

end

function get_particle_number_MPO(sites)

    N = length(sites)

    opsum = OpSum()

    for n in 1:N
        
        opsum += 0.5*(-1)^(n-1),"Z",n

    end

    opsum += 0.5*N,"Id",1

    mpo = MPO(opsum, sites)

    return mpo

end

function get_initial_zero_charge_MPO(sites, state)

    """
    Prepare the mpo = |state><state| where |state> is a basis state given as a list of integers 1 and 2 e.g. state = [1,1,2,1] which would be the state |0010>
    """

    N = length(sites)
    mpo = MPO(sites)

    for i=1:N

        if i == 1 || i == N

            s, sp = inds(mpo[i]; :tags => "Site")
            l = inds(mpo[i]; :tags => "Link")[1]
            mpo[i][s => state[i], sp => state[i], l => 1] = 1.0

        else

            s, sp = inds(mpo[i]; :tags => "Site")
            l1, l2 = inds(mpo[i]; :tags => "Link")
            mpo[i][s => state[i], sp => state[i], l1 => 1, l2 => 1] = 1.0

        end

    end

    return mpo

end

function get_MPO_site_canonical_form(mpo_site, which_site, mpo_site_index)

    mpo_site_dag = mpo_site'
    noprime!(mpo_site_dag; :plev => 2)
    mpo_site_dag = dag(mpo_site_dag'; :tags => "Site")
    tmp = mpo_site_dag * mpo_site

    # Checking for LCF
    if which_site == "first"

        tmp_l = tmp * dag(delta(inds(tmp; :tags => "Site")))
        res = ITensors.Array(tmp_l, inds(tmp_l)...)
        s = size(res)
        if length(s) == 0
            l = 1
        else
            l = s[1]
        end
        is_left = isapprox(res, I(l))

    elseif which_site == "last"

        tmp_l = tmp * dag(delta(inds(tmp; :tags => "Site")))
        tmp_l = tmp_l * dag(delta(inds(tmp_l; :tags => "Link")))
        res = ITensors.Array(tmp_l, inds(tmp_l)...)
        s = size(res)
        if length(s) == 0
            l = 1
        else
            l = s[1]
        end
        is_left = isapprox(res, I(l))
        
    else

        tmp_l = tmp * dag(delta(inds(tmp; :tags => "Site")))
        tmp_l = tmp_l * dag(delta(inds(tmp_l; :tags => "Link,l=$(mpo_site_index-1)")))
        res = ITensors.Array(tmp_l, inds(tmp_l)...)
        s = size(res)
        if length(s) == 0
            l = 1
        else
            l = s[1]
        end
        is_left = isapprox(res, I(l))

    end

    # Checking for RCF
    if which_site == "first"

        tmp = tmp * dag(delta(inds(tmp; :tags => "Site")))
        tmp = tmp * dag(delta(inds(tmp; :tags => "Link")))
        res = ITensors.Array(tmp, inds(tmp)...)
        s = size(res)
        if length(s) == 0
            l = 1
        else
            l = s[1]
        end
        is_right = isapprox(res, I(l))

    elseif which_site == "last"

        tmp = tmp * dag(delta(inds(tmp; :tags => "Site")))
        res = ITensors.Array(tmp, inds(tmp)...)
        s = size(res)
        if length(s) == 0
            l = 1
        else
            l = s[1]
        end
        is_right = isapprox(res, I(l))
        
    else

        tmp = tmp * dag(delta(inds(tmp; :tags => "Site")))
        tmp = tmp * dag(delta(inds(tmp; :tags => "Link,l=$(mpo_site_index)")))
        res = ITensors.Array(tmp, inds(tmp)...)
        s = size(res)
        if length(s) == 0
            l = 1
        else
            l = s[1]
        end
        is_right = isapprox(res, I(l))

    end

    if is_left
        if is_right 
            return "L/R"
        else
            return "L"
        end
    elseif is_right
        return "R"
    else
        return "N"
    end

end

function get_MPO_canonical_form(mpo)
    
    res = String[]
    N = length(mpo)
    for (mpo_site_index, mpo_site) in enumerate(mpo)
        if mpo_site_index == 1
            push!(res, get_MPO_site_canonical_form(mpo_site, "first", mpo_site_index))
        elseif mpo_site_index == N
            push!(res, get_MPO_site_canonical_form(mpo_site, "last", mpo_site_index))
        else
            push!(res, get_MPO_site_canonical_form(mpo_site, "none", mpo_site_index))
        end
    end
    return res
end

function ishermitian(mpo; tol = 1e-14)

    return norm(mpo - dag(swapprime(mpo, 0, 1))) < tol

end

function ispositive(mpo; tol = 1e-14)

    """
    We check this by using DMRG to find the ground state energy of the "density matrix" MPO
    """

    n = length(mpo)
    dmrg_tol = tol
    cutoff = tol
    sites = dag(reduce(vcat, siteinds(mpo; :plev => 0)))
    state = [isodd(i) ? "0" : "1" for i = 1:n]
    psi0 = randomMPS(sites, state)
    obs = DMRGObserver(;energy_tol = dmrg_tol)
    nsweeps = 100
    dmrg_energy, _ = dmrg(mpo, psi0; nsweeps, cutoff, observer=obs, outputlevel=1)

    return abs(dmrg_energy) < tol

end

function mpo_to_matrix(mpo)

    n = length(mpo)
    a = contract(mpo)
    a = Array(a, inds(a; :plev => 1)..., inds(a; :plev => 0)...)
    a = reshape(a, 2^n, 2^n)

    return a

end

function get_aH_Hamiltonian_sparse_matrix(N, x, ma, l_0, lambda)

    """

    This gives aH as a sparse Hamiltonian and here we will also add the penalty term

    """

    eye(n::Int64) = sparse(I, n, n);

    H = spzeros(2^(N), 2^(N))

    # Kinetic term
    for n=1:N-1
        H += (1/4)*get_op(["X", "X"], [n, n+1], N)
        H += (1/4)*get_op(["Y", "Y"], [n, n+1], N)
    end

    # Long range ZZ interaction term
    for n = 1:N-1
        for m = n+1:N
            H += (0.25/x)*(N - m + lambda)*get_op(["Z", "Z"], [n, m], N)
        end
    end

    # # Mass term
    for n=1:N
        H += (0.5*ma)*((-1)^(n-1))*get_op(["Z"], [n], N)
    end

    # Electric single Z term
    for n=1:N-1
        H += ((N/8 - 0.25*ceil((n-1)/2) + l_0*(N-n)/2)/x)*get_op(["Z"], [n], N)
    end

    # # Constant term
    H += ((l_0^2)*(N-1)/(2*x) + (l_0*N/(4*x)) + (N^2)/(16*x) + lambda*N/(8*x))*eye(2^N)

    return H
end

function environment_correlator(type, n, m, aD_0, sigma_over_a)

    if type == "constant"
        return aD_0
    elseif type == "delta"
        if n == m
            return aD_0
        else
            return 0.0
        end
    else # gaussian case
        return exp(-0.5*(1/sigma_over_a)^2*(n-m)^2)
    end

end

function get_kinetic_part_aH_Hamiltonian_sparse_matrix(N)
    
    eye(n::Int64) = sparse(I, n, n);

    Hk = spzeros(2^(N), 2^(N))

    # Kinetic term
    for n=1:N-1
        Hk += (1/4)*get_op(["X", "X"], [n, n+1], N)
        Hk += (1/4)*get_op(["Y", "Y"], [n, n+1], N)
    end

    return Hk

end

function my_kron(A, B)
    
    m, n = size(A)
    p, q = size(B)

    C = zeros(ComplexF64, m * p, n * q)

    for i in 1:p
        for j in 1:q
            C[(i-1)*m+1 : i*m, (j-1)*n+1 : j*n] = A * B[i, j]
        end
    end

    return C
end

function get_op(ops, positions, N; reverse_flag = true)

    op_dict = Dict("X" => sparse([0 1; 1 0]), "Y" => sparse([0 -1im; 1im 0]), "Z" => sparse([1 0; 0 -1]))
    zipped = TupleTools.sort(Tuple(zip(1:length(ops), positions, ops)); by = x -> x[2])
    old_positions = [element[2] for element in zipped] 
    old_ops = [element[3] for element in zipped]

    positions = []
    ops = []

    if length(Set(old_positions)) != length(old_positions) # case where we have duplicate positions
        
        flag = false

        for (idx, pos) in enumerate(old_positions)

            if flag

                flag = false
                continue

            end

            if idx != length(old_positions)

                if pos != old_positions[idx+1]

                    push!(positions, pos)
                    push!(ops, op_dict[old_ops[idx]])

                else

                    push!(positions, pos)
                    push!(ops, op_dict[old_ops[idx]]*op_dict[old_ops[idx+1]])
                    flag = true

                end

            else

                push!(positions, pos)
                push!(ops, op_dict[old_ops[idx]])

            end

        end

    else

        for (idx, pos) in enumerate(old_positions)

            push!(positions, pos)
            push!(ops, op_dict[old_ops[idx]])
        
        end

    end

    eye(n) = sparse(I, n, n)

    res = eye(1)

    for (i, pos) in enumerate(positions)

        if i == 1
            how_many_I_before = pos-1
        else
            how_many_I_before = pos - positions[i-1] - 1
        end

        pos = positions[i]
        op = ops[i]
    
        if reverse_flag
            res = my_kron(res, eye(2^how_many_I_before))
            res = my_kron(res, op)
        else
            res = kron(res, eye(2^how_many_I_before))
            res = kron(res, op)
        end

    end

    if reverse_flag
        res = my_kron(res, eye(2^(N - positions[end])))
    else
        res = kron(res, eye(2^(N - positions[end])))
    end

    return res

end

function get_LdagL_sparse_matrix(N, n, m, aT)

    # This is returning a^2 * L^\dagger_n L_m
    
    eye(n::Int64) = sparse(I, n, n);

    res = spzeros(2^(N), 2^(N))

    res += 0.25*(-1)^(n+m) * get_op(["Z", "Z"], [n, m], N) 
    res += 0.5*(-1)^(n+m) * get_op(["Z"], [n], N) 
    res += 0.25*(-1)^(n+m) * eye(2^N)

    if (n != 1)
        res += (-1im*(-1)^(n + m)/(32*aT)) * get_op(["X", "Y", "Z"], [n-1, n, m], N) 
        res += (1im*(-1)^(n + m)/(32*aT)) * get_op(["Y", "X", "Z"], [n-1, n, m], N) 
    end
    
    if (n != N)
        res += (1im*(-1)^(n + m)/(32*aT)) * get_op(["X", "Y", "Z"], [n, n+1, m], N) 
        res += (-1im*(-1)^(n + m)/(32*aT)) * get_op(["Y", "X", "Z"], [n, n+1, m], N) 
    end

    if (m != 1)
        res += (1im*(-1)^(n + m)/(32*aT)) * get_op(["Z", "X", "Y"], [n, m-1, m], N) 
        res += (-1im*(-1)^(n + m)/(32*aT)) * get_op(["Z", "Y", "X"], [n, m-1, m], N) 
    end
    
    if (m != N)
        res += (-1im*(-1)^(n + m)/(32*aT)) * get_op(["Z", "X", "Y"], [n, m, m+1], N) 
        res += (1im*(-1)^(n + m)/(32*aT)) * get_op(["Z", "Y", "X"], [n, m, m+1], N) 
    end

    if (n != 1) && (m != 1)
        res += (-(-1)^(n + m)/(256*aT^2)) * get_op(["Y", "X", "X", "Y"], [n-1, n, m-1, m], N) 
        res += ((-1)^(n + m)/(256*aT^2)) * get_op(["Y", "X", "Y", "X"], [n-1, n, m-1, m], N)
        res += ((-1)^(n + m)/(256*aT^2)) * get_op(["X", "Y", "X", "Y"], [n-1, n, m-1, m], N) 
        res += (-(-1)^(n + m)/(256*aT^2)) * get_op(["X", "Y", "Y", "X"], [n-1, n, m-1, m], N)  
    end

    if (n != 1) && (m != N)
        res += (-(-1)^(n + m)/(256*aT^2)) * get_op(["Y", "X", "Y", "X"], [n-1, n, m, m+1], N) 
        res += ((-1)^(n + m)/(256*aT^2)) * get_op(["Y", "X", "X", "Y"], [n-1, n, m, m+1], N) 
        res += (-(-1)^(n + m)/(256*aT^2)) * get_op(["X", "Y", "X", "Y"], [n-1, n, m, m+1], N) 
        res += ((-1)^(n + m)/(256*aT^2)) * get_op(["X", "Y", "Y", "X"], [n-1, n, m, m+1], N) 
    end

    if (n != N) && (m != 1)
        res += (-(-1)^(n + m)/(256*aT^2)) * get_op(["X", "Y", "X", "Y"], [n, n+1, m-1, m], N) 
        res += ((-1)^(n + m)/(256*aT^2)) * get_op(["X", "Y", "Y", "X"], [n, n+1, m-1, m], N) 
        res += (-(-1)^(n + m)/(256*aT^2)) * get_op(["Y", "X", "Y", "X"], [n, n+1, m-1, m], N) 
        res += ((-1)^(n + m)/(256*aT^2)) * get_op(["Y", "X", "X", "Y"], [n, n+1, m-1, m], N) 
    end
    
    if (n != N) && (m != N)
        res += (-(-1)^(n + m)/(256*aT^2)) * get_op(["X", "Y", "Y", "X"], [n, n+1, m, m+1], N) 
        res += ((-1)^(n + m)/(256*aT^2)) * get_op(["X", "Y", "X", "Y"], [n, n+1, m, m+1], N) 
        res += (-(-1)^(n + m)/(256*aT^2)) * get_op(["Y", "X", "X", "Y"], [n, n+1, m, m+1], N) 
        res += ((-1)^(n + m)/(256*aT^2)) * get_op(["Y", "X", "Y", "X"], [n, n+1, m, m+1], N) 
    end

    return res

end

function get_Lindblad_jump_operator_old(m, aT, sites)

    N = length(sites)

    opsum = OpSum()

    opsum += 0.5*(-1)^m,"Z",m

    opsum += 0.5*(-1)^m,"Id",1

    if m != 1
        opsum += ( 1im*(-1)^m/(16*aT)),"X",m-1,"Y",m
        opsum += (-1im*(-1)^m/(16*aT)),"Y",m-1,"X",m
    end
    
    if m != N
        opsum += (-1im*(-1)^m/(16*aT)),"X",m,"Y",m+1 
        opsum += ( 1im*(-1)^m/(16*aT)),"Y",m,"X",m+1
    end

    return MPO(opsum, sites)

end

function get_Lindblad_jump_operator(m, aT, sites)

    N = length(sites)

    opsum = OpSum()

    opsum += 0.5*(-1)^m,"Z",m

    opsum += 0.5*(-1)^m,"Id",1

    c = (-1)^m/(8*aT)

    if m != 1
        opsum +=  c,"S-",m-1,"S+",m
        opsum += -c,"S+",m-1,"S-",m
    end
    
    if m != N
        opsum += -c,"S-",m,"S+",m+1 
        opsum +=  c,"S+",m,"S-",m+1
    end

    return MPO(opsum, sites)

end

function get_LdagL_old(n, m, aT, sites)

    N = length(sites)
    
    opsum = OpSum()

    opsum += 0.25*(-1)^(n+m),"Z",n,"Z",m

    opsum += 0.5*(-1)^(n+m),"Z",n

    opsum += 0.25*(-1)^(n+m),"Id",1

    if (n != 1)
        opsum += (-1im*(-1)^(n + m)/(32*aT)),"X",n-1,"Y",n,"Z",m
        opsum += (1im*(-1)^(n + m)/(32*aT)),"Y",n-1,"X",n,"Z",m
    end
    
    if (n != N)
        opsum += (1im*(-1)^(n + m)/(32*aT)),"X",n,"Y",n+1,"Z",m
        opsum += (-1im*(-1)^(n + m)/(32*aT)),"Y",n,"X",n+1,"Z",m
    end

    if (m != 1)
        opsum += (1im*(-1)^(n + m)/(32*aT)),"Z",n,"X",m-1,"Y",m
        opsum += (-1im*(-1)^(n + m)/(32*aT)),"Z",n,"Y",m-1,"X",m
    end
    
    if (m != N)
        opsum += (-1im*(-1)^(n + m)/(32*aT)),"Z",n,"X",m,"Y",m+1
        opsum += (1im*(-1)^(n + m)/(32*aT)),"Z",n,"Y",m,"X",m+1
    end

    if (n != 1) && (m != 1)
        opsum += (-(-1)^(n + m)/(256*aT^2)),"Y",n-1,"X",n,"X",m-1,"Y",m
        opsum += ((-1)^(n + m)/(256*aT^2)),"Y",n-1,"X",n,"Y",m-1,"X",m
        opsum += ((-1)^(n + m)/(256*aT^2)),"X",n-1,"Y",n,"X",m-1,"Y",m
        opsum += (-(-1)^(n + m)/(256*aT^2)),"X",n-1,"Y",n,"Y",m-1,"X",m
    end

    if (n != 1) && (m != N)
        opsum += (-(-1)^(n + m)/(256*aT^2)),"Y",n-1,"X",n,"Y",m,"X",m+1
        opsum += ((-1)^(n + m)/(256*aT^2)),"Y",n-1,"X",n,"X",m,"Y",m+1
        opsum += (-(-1)^(n + m)/(256*aT^2)),"X",n-1,"Y",n,"X",m,"Y",m+1
        opsum += ((-1)^(n + m)/(256*aT^2)),"X",n-1,"Y",n,"Y",m,"X",m+1
    end

    if (n != N) && (m != 1)
        opsum += (-(-1)^(n + m)/(256*aT^2)),"X",n,"Y",n+1,"X",m-1,"Y",m
        opsum += ((-1)^(n + m)/(256*aT^2)),"X",n,"Y",n+1,"Y",m-1,"X",m
        opsum += (-(-1)^(n + m)/(256*aT^2)),"Y",n,"X",n+1,"Y",m-1,"X",m
        opsum += ((-1)^(n + m)/(256*aT^2)),"Y",n,"X",n+1,"X",m-1,"Y",m
    end
    
    if (n != N) && (m != N)
        opsum += (-(-1)^(n + m)/(256*aT^2)),"X",n,"Y",n+1,"Y",m,"X",m+1
        opsum += ((-1)^(n + m)/(256*aT^2)),"X",n,"Y",n+1,"X",m,"Y",m+1
        opsum += (-(-1)^(n + m)/(256*aT^2)),"Y",n,"X",n+1,"X",m,"Y",m+1
        opsum += ((-1)^(n + m)/(256*aT^2)),"Y",n,"X",n+1,"Y",m,"X",m+1
    end

    return MPO(opsum, sites)

end

function get_LdagL(n, m, aT, sites)

    N = length(sites)
    
    opsum = OpSum()

    opsum += 0.25*(-1)^(n+m),"Z",n,"Z",m

    opsum += 0.5*(-1)^(n+m),"Z",n

    opsum += 0.25*(-1)^(n+m),"Id",1

    if (n != 1)
        opsum += (-(-1)^(n + m)/(16*aT)),"S-",n-1,"S+",n,"Z",m
        opsum += ((-1)^(n + m)/(16*aT)),"S+",n-1,"S-",n,"Z",m
    end
    
    if (n != N)
        opsum += ((-1)^(n + m)/(16*aT)),"S-",n,"S+",n+1,"Z",m
        opsum += (-(-1)^(n + m)/(16*aT)),"S+",n,"S-",n+1,"Z",m
    end

    if (m != 1)
        opsum += ((-1)^(n + m)/(16*aT)),"Z",n,"S-",m-1,"S+",m
        opsum += (-(-1)^(n + m)/(16*aT)),"Z",n,"S+",m-1,"S-",m
    end
    
    if (m != N)
        opsum += (-(-1)^(n + m)/(16*aT)),"Z",n,"S-",m,"S+",m+1
        opsum += ((-1)^(n + m)/(16*aT)),"Z",n,"S+",m,"S-",m+1
    end

    if (n != 1) && (m != 1)
        opsum += (-(-1)^(n + m)/(64*aT^2)),"S+",n-1,"S-",n,"S+",m-1,"S-",m
        opsum += ((-1)^(n + m)/(64*aT^2)),"S+",n-1,"S-",n,"S-",m-1,"S+",m
        opsum += ((-1)^(n + m)/(64*aT^2)),"S-",n-1,"S+",n,"S+",m-1,"S-",m
        opsum += (-(-1)^(n + m)/(64*aT^2)),"S-",n-1,"S+",n,"S-",m-1,"S+",m
    end

    if (n != 1) && (m != N)
        opsum += (-(-1)^(n + m)/(64*aT^2)),"S+",n-1,"S-",n,"S-",m,"S+",m+1
        opsum += ((-1)^(n + m)/(64*aT^2)),"S+",n-1,"S-",n,"S+",m,"S-",m+1
        opsum += (-(-1)^(n + m)/(64*aT^2)),"S-",n-1,"S+",n,"S+",m,"S-",m+1
        opsum += ((-1)^(n + m)/(64*aT^2)),"S-",n-1,"S+",n,"S-",m,"S+",m+1
    end

    if (n != N) && (m != 1)
        opsum += (-(-1)^(n + m)/(64*aT^2)),"S+",n,"S-",n+1,"S-",m-1,"S+",m
        opsum += ((-1)^(n + m)/(64*aT^2)),"S+",n,"S-",n+1,"S+",m-1,"S-",m
        opsum += (-(-1)^(n + m)/(64*aT^2)),"S-",n,"S+",n+1,"S+",m-1,"S-",m
        opsum += ((-1)^(n + m)/(64*aT^2)),"S-",n,"S+",n+1,"S-",m-1,"S+",m
    end
    
    if (n != N) && (m != N)
        opsum += (-(-1)^(n + m)/(64*aT^2)),"S+",n,"S-",n+1,"S+",m,"S-",m+1
        opsum += ((-1)^(n + m)/(64*aT^2)),"S+",n,"S-",n+1,"S-",m,"S+",m+1
        opsum += (-(-1)^(n + m)/(64*aT^2)),"S-",n,"S+",n+1,"S-",m,"S+",m+1
        opsum += ((-1)^(n + m)/(64*aT^2)),"S-",n,"S+",n+1,"S+",m,"S-",m+1
    end

    return MPO(opsum, sites)

end

function hermitian_conjugate_mpo(mpo)

    return dag(swapprime(mpo, 0, 1))

end

function transpose_mpo(mpo)

    return swapprime(dag(conj(mpo)), 0 => 1)

end

function get_Lindblad_dissipative_part(aD_0, sigma_over_a, env_corr_type, aT, sites)

    # This returns sum over n, m aD_0 * env_corr * a^2 * L^\dagger_n * L_m

    N = length(sites)

    mpo = aD_0 * environment_correlator(env_corr_type, 1, 1, aD_0, sigma_over_a) * get_LdagL(1, 1, aT, sites)

    for n=1:N
        for m=1:N

            if !((n == 1) && (m == 1))
                mpo += aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * get_LdagL(n, m, aT, sites)
            end

        end
    end

    return mpo

end

function get_double_size_Lindblad_operator(N, sites, x, ma, l_0, lambda, aD_0, sigma_over_a, aT, env_corr_type)

    opsum = OpSum()

    # The -iH_left which goes from 1 to N
    for n in 1:N-1
        for m in n+1:N
            # Long range ZZ interaction term
            opsum += -1im*0.25*(1/x)*(N-m),"Z",n,"Z",m
        end
        # Kinetic term
        opsum += -1im*0.25,"X",n,"X",n+1
        opsum += -1im*0.25,"Y",n,"Y",n+1
        # Single Z terms
        opsum += -1im*(1/x)*(N/8 - 0.25*ceil((n-1)/2) + l_0*(N-n)/2),"Z",n
        opsum += -1im*(0.5*ma*(-1)^(n-1)),"Z",n
    end
    opsum += -1im*(0.5*ma*(-1)^(N-1)),"Z",N
    opsum += -1im*((l_0^2)*(N-1)/(2*x) + (l_0*N)/(4*x) + (N^2)/(16*x)),"Id",1

    # # The +iH_right which goes from N+1 to 2N
    for n in 1:N-1
        for m in n+1:N
            # Long range ZZ interaction term
            opsum += 1im*0.25*(1/x)*(N-m),"Z",n+N,"Z",m+N
        end
        # Kinetic term
        opsum += 1im*0.25,"X",n+N,"X",n+1+N
        opsum += 1im*0.25,"Y",n+N,"Y",n+1+N
        # Single Z terms
        opsum += 1im*(1/x)*(N/8 - 0.25*ceil((n-1)/2) + l_0*(N-n)/2),"Z",n+N
        opsum += 1im*(0.5*ma*(-1)^(n-1)),"Z",n+N
    end
    opsum += 1im*(0.5*ma*(-1)^(N-1)),"Z",N+N
    opsum += 1im*((l_0^2)*(N-1)/(2*x) + (l_0*N)/(4*x) + (N^2)/(16*x)),"Id",1+N

    h1 = MPO(opsum, sites; cutoff = 0)

    if aD_0 != 0

        for n=1:N
            for m=1:N

                # The L_m_left L_n_right
                opsum = OpSum()
                opsum += 0.5*(-1)^m,"Z",m
                opsum += 0.5*(-1)^m,"Id",1
                if m != 1
                    opsum += ( 1im*(-1)^m/(16*aT)),"X",m-1,"Y",m
                    opsum += (-1im*(-1)^m/(16*aT)),"Y",m-1,"X",m
                end
                if m != N
                    opsum += (-1im*(-1)^m/(16*aT)),"X",m,"Y",m+1 
                    opsum += ( 1im*(-1)^m/(16*aT)),"Y",m,"X",m+1
                end 
                h2 = MPO(aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a)*opsum, sites; cutoff = 0)
                
                opsum = OpSum()
                opsum += 0.5*(-1)^n,"Z",n+N
                opsum += 0.5*(-1)^n,"Id",1
                if n != 1
                    opsum += (-1im*(-1)^n/(16*aT)),"X",n-1,"Y",n+N
                    opsum += ( 1im*(-1)^n/(16*aT)),"Y",n-1,"X",n+N
                end
                if n != N
                    opsum += ( 1im*(-1)^n/(16*aT)),"X",n,"Y",n+N+1 
                    opsum += (-1im*(-1)^n/(16*aT)),"Y",n,"X",n+N+1
                end 
                h3 = hermitian_conjugate_mpo(MPO(aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a)*opsum, sites; cutoff = 0))
                h1 = add(h1, apply(h2, h3, cutoff = 0); cutoff = 0)
                
                # The Ldag_n_L_m_left
                opsum = OpSum()
                opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * 0.25*(-1)^(n+m),"Z",n,"Z",m
                opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * 0.5*(-1)^(n+m),"Z",n
                opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * 0.25*(-1)^(n+m),"Id",1
                if (n != 1)
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (-1im*(-1)^(n + m)/(32*aT)),"X",n-1,"Y",n,"Z",m
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (1im*(-1)^(n + m)/(32*aT)),"Y",n-1,"X",n,"Z",m
                end            
                if (n != N)
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (1im*(-1)^(n + m)/(32*aT)),"X",n,"Y",n+1,"Z",m
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (-1im*(-1)^(n + m)/(32*aT)),"Y",n,"X",n+1,"Z",m
                end
                if (m != 1)
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (1im*(-1)^(n + m)/(32*aT)),"Z",n,"X",m-1,"Y",m
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (-1im*(-1)^(n + m)/(32*aT)),"Z",n,"Y",m-1,"X",m
                end            
                if (m != N)
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (-1im*(-1)^(n + m)/(32*aT)),"Z",n,"X",m,"Y",m+1
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (1im*(-1)^(n + m)/(32*aT)),"Z",n,"Y",m,"X",m+1
                end
                if (n != 1) && (m != 1)
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (-(-1)^(n + m)/(256*aT^2)),"Y",n-1,"X",n,"X",m-1,"Y",m
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * ((-1)^(n + m)/(256*aT^2)),"Y",n-1,"X",n,"Y",m-1,"X",m
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * ((-1)^(n + m)/(256*aT^2)),"X",n-1,"Y",n,"X",m-1,"Y",m
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (-(-1)^(n + m)/(256*aT^2)),"X",n-1,"Y",n,"Y",m-1,"X",m
                end
                if (n != 1) && (m != N)
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (-(-1)^(n + m)/(256*aT^2)),"Y",n-1,"X",n,"Y",m,"X",m+1
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * ((-1)^(n + m)/(256*aT^2)),"Y",n-1,"X",n,"X",m,"Y",m+1
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (-(-1)^(n + m)/(256*aT^2)),"X",n-1,"Y",n,"X",m,"Y",m+1
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * ((-1)^(n + m)/(256*aT^2)),"X",n-1,"Y",n,"Y",m,"X",m+1
                end
                if (n != N) && (m != 1)
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (-(-1)^(n + m)/(256*aT^2)),"X",n,"Y",n+1,"X",m-1,"Y",m
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * ((-1)^(n + m)/(256*aT^2)),"X",n,"Y",n+1,"Y",m-1,"X",m
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (-(-1)^(n + m)/(256*aT^2)),"Y",n,"X",n+1,"Y",m-1,"X",m
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * ((-1)^(n + m)/(256*aT^2)),"Y",n,"X",n+1,"X",m-1,"Y",m
                end            
                if (n != N) && (m != N)
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (-(-1)^(n + m)/(256*aT^2)),"X",n,"Y",n+1,"Y",m,"X",m+1
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * ((-1)^(n + m)/(256*aT^2)),"X",n,"Y",n+1,"X",m,"Y",m+1
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (-(-1)^(n + m)/(256*aT^2)),"Y",n,"X",n+1,"X",m,"Y",m+1
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * ((-1)^(n + m)/(256*aT^2)),"Y",n,"X",n+1,"Y",m,"X",m+1
                end
                
                # The Ldag_n_L_m_right
                opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * 0.25*(-1)^(n+m),"Z",n+N,"Z",m+N
                opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * 0.5*(-1)^(n+m),"Z",n+N
                opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * 0.25*(-1)^(n+m),"Id",1+N
                if (n != 1)
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (-1im*(-1)^(n + m)/(32*aT)),"X",n-1+N,"Y",n+N,"Z",m+N
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (1im*(-1)^(n + m)/(32*aT)),"Y",n-1+N,"X",n+N,"Z",m+N
                end            
                if (n != N)
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (1im*(-1)^(n + m)/(32*aT)),"X",n+N,"Y",n+1+N,"Z",m+N
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (-1im*(-1)^(n + m)/(32*aT)),"Y",n+N,"X",n+1+N,"Z",m+N
                end
                if (m != 1)
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (1im*(-1)^(n + m)/(32*aT)),"Z",n+N,"X",m-1+N,"Y",m+N
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (-1im*(-1)^(n + m)/(32*aT)),"Z",n+N,"Y",m-1+N,"X",m+N
                end            
                if (m != N)
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (-1im*(-1)^(n + m)/(32*aT)),"Z",n+N,"X",m+N,"Y",m+1+N
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (1im*(-1)^(n + m)/(32*aT)),"Z",n+N,"Y",m+N,"X",m+1+N
                end
                if (n != 1) && (m != 1)
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (-(-1)^(n + m)/(256*aT^2)),"Y",n-1+N,"X",n+N,"X",m-1+N,"Y",m+N
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * ((-1)^(n + m)/(256*aT^2)),"Y",n-1+N,"X",n+N,"Y",m-1+N,"X",m+N
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * ((-1)^(n + m)/(256*aT^2)),"X",n-1+N,"Y",n+N,"X",m-1+N,"Y",m+N
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (-(-1)^(n + m)/(256*aT^2)),"X",n-1+N,"Y",n+N,"Y",m-1+N,"X",m+N
                end
                if (n != 1) && (m != N)
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (-(-1)^(n + m)/(256*aT^2)),"Y",n-1+N,"X",n+N,"Y",m+N,"X",m+1+N
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * ((-1)^(n + m)/(256*aT^2)),"Y",n-1+N,"X",n+N,"X",m+N,"Y",m+1+N
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (-(-1)^(n + m)/(256*aT^2)),"X",n-1+N,"Y",n+N,"X",m+N,"Y",m+1+N
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * ((-1)^(n + m)/(256*aT^2)),"X",n-1+N,"Y",n+N,"Y",m+N,"X",m+1+N
                end
                if (n != N) && (m != 1)
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (-(-1)^(n + m)/(256*aT^2)),"X",n+N,"Y",n+1+N,"X",m-1+N,"Y",m+N
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * ((-1)^(n + m)/(256*aT^2)),"X",n+N,"Y",n+1+N,"Y",m-1+N,"X",m+N
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (-(-1)^(n + m)/(256*aT^2)),"Y",n+N,"X",n+1+N,"Y",m-1+N,"X",m+N
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * ((-1)^(n + m)/(256*aT^2)),"Y",n+N,"X",n+1+N,"X",m-1+N,"Y",m+N
                end            
                if (n != N) && (m != N)
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (-(-1)^(n + m)/(256*aT^2)),"X",n+N,"Y",n+1+N,"Y",m+N,"X",m+1+N
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * ((-1)^(n + m)/(256*aT^2)),"X",n+N,"Y",n+1+N,"X",m+N,"Y",m+1+N
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (-(-1)^(n + m)/(256*aT^2)),"Y",n+N,"X",n+1+N,"X",m+N,"Y",m+1+N
                    opsum += -0.5 * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * ((-1)^(n + m)/(256*aT^2)),"Y",n+N,"X",n+1+N,"Y",m+N,"X",m+1+N
                end
                h1 = add(h1, MPO(opsum, sites; cutoff = 0); cutoff = 0)

            end
        end

    end

    return h1 
    
end

function get_Lindblad_jump_operator_sparse_matrix(N, m, aT)

    """
    This is getting a*L_m where L_m is the jump operator at site m
    """

    eye(n::Int64) = sparse(I, n, n);
    res = spzeros(2^(N), 2^(N))        
    
    res += 0.5*((-1)^m)*get_op(["Z"], [m], N)
    res += 0.5*((-1)^m)*eye(2^N)
    
    if m != 1
        res += (1im*(-1)^m/(16*aT))*get_op(["X", "Y"], [m-1, m], N)
        res += (-1im*(-1)^m/(16*aT))*get_op(["Y", "X"], [m-1, m], N)
    end
    
    if m != N
        res += (-1im*(-1)^m/(16*aT))*get_op(["X", "Y"], [m, m+1], N)
        res += (1im*(-1)^m/(16*aT))*get_op(["Y", "X"], [m, m+1], N)
    end

    return res

end

function get_Lindblad_sparse_matrix(N, x, ma, l_0, lambda, aD_0, sigma_over_a, aT, env_corr_type)

    """

    This gets the Lindblad operator as a sparse matrix in the purified version
    see eg eq 15 16 in Numerical evaluation of two-time correlation functions
    in open quantum systems with matrix product state methods:
    a comparison - kollath et al

    """

    L = spzeros(2^(2*N), 2^(2*N))

    eye(n::Int64) = sparse(I, n, n);

    H = get_aH_Hamiltonian_sparse_matrix(N, x, ma, l_0, lambda)

    # Unitary part of Lindbladian
    L += -1im * my_kron(H, eye(2^N)) + 1im * my_kron(eye(2^N), transpose(H)) 

    for n in 1:N
        for m in 1:N

            tmp1 = get_Lindblad_jump_operator_sparse_matrix(N, n, aT)
            tmp2 = get_Lindblad_jump_operator_sparse_matrix(N, m, aT)

            tmp3 = tmp1' * tmp2 # the dash is the dagger
            
            L += aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * ((my_kron(tmp2, transpose(tmp1'))) - 0.5*(my_kron(tmp3, eye(2^N))) -0.5*(my_kron(eye(2^N), transpose(tmp3))))

        end
    end

    return L

end

function get_entanglement_entropy_vector(state, trace_indices, dof_list, tol=1e-12)
    
    # Inputs:
    # state = numpy array statevector
    # trace_indices = list of indices of sites to be traced out
    # dof_list = list of number of degrees of freedom per site for all sites
    # tol = any eigenvalue of the reduced density matrix that is smaller than this tolerance will be neglected
    # Outputs:
    # ee, rho_reduced = entanglement entropy and reduced density matrix
    
    # Make sure input is in the right type form and state is normalized to 1
    state = state / norm(state)
    
    # Just a simple list containing the indices from 1 to N where N is the total number of sites
    site_indices = 1:length(dof_list)
    
    # The dimension of each index to be traced
    trace_dof = [dof_list[i] for i in trace_indices]

    # List containing the indices of the sites not to be traced
    untraced_indices = setdiff(site_indices, trace_indices)

    # The dimension of each index in the list of untraced indices
    untraced_dof = [dof_list[i] for i in untraced_indices]

    # Reshape statevector into tensor of rank N with each index having some degrees of freedom specified by the dof_list
    # for example if it is a spin-1/2 chain then each site has 2 degrees of freedom and the dof_list should be [2]*N = [2, 2, 2, 2, ..., 2]
    state = reshape(state, dof_list)

    # Revert indices of the reshaped tensor to meet the convention of qiskit eg for 4 qubits q3 q2 q1 q0 is the ordering
    state = permutedims(state, site_indices[end:-1:1]) # TODO: check whether this is necessary for julia

    # Permute the indices of the rank N tensor so the untraced indices are placed on the left and the ones to be traced on the right
    state = permutedims(state, vcat(untraced_indices, trace_indices))

    # Reshape the rank N tensor into a matrix where you merge the untraced indices into 1 index and you merge the traced indices into 1 index
    # if the former index is called I and the latter J then we have state_{I, J}
    state = reshape(state, (prod(untraced_dof), prod(trace_dof)))

    # The reduced density matrix is given by state_{I, J}*state_complex_conjugated_{J, K}, so we see from here that the indices to be
    # traced out ie the ones contained in the merged big index J are summed over in the matrix multiplication
    rho_reduced = state * adjoint(state)

    evals = eigen(rho_reduced).values

    ee = sum(-real(eval)*log(real(eval)) for eval in evals if real(eval) >= tol)

    return ee, rho_reduced # return both the entanglement entropy and the reduced density matrix
end

function apply_taylor_part(rho, cutoff, tau, sites, x, l_0, ma, aD_0, sigma_over_a, env_corr_type, aT)

    H_T = get_H_taylor(sites, x, l_0, ma)

    tmp = -0.5*1im*tau*H_T

    # rho + rho * idt/2 H_T * rho - idt/2 H_T * rho
    rho_final = rho + apply(rho, hermitian_conjugate_mpo(tmp); cutoff = cutoff) + apply(tmp, rho; cutoff = cutoff)

    # second order term in the taylor expansion only for the Hamiltonian part
    # rho_final += (-tau^2/8)*apply(H_T, apply(H_T, rho; cutoff = cutoff); cutoff = cutoff) + (-tau^2/8)*apply(rho, apply(H_T, H_T; cutoff = cutoff); cutoff = cutoff) + (tau^2/4)*apply(H_T, apply(rho, H_T; cutoff = cutoff); cutoff = cutoff)

    if aD_0 != 0

        # - 0.5 * dt * 0.5 * L_n^\dagger L_m
        LdagnLm = -0.5 * tau * 0.5 * get_Lindblad_dissipative_part(aD_0, sigma_over_a, env_corr_type, aT, sites)

        # sum over n and m: - 0.5 * dt/2 * L_n^\dagger L_m * rho
        rho_final += apply(LdagnLm, rho; cutoff = cutoff) 
        
        # sum over n and m: - 0.5 * dt/2 * rho * L_n^\dagger L_m 
        rho_final += apply(rho, LdagnLm; cutoff = cutoff) 

        # sum over n and m: aD_0 * f(a(n-m)) * L_m * rho * L_n^\dagger
        for n=1:N
            for m=1:N

                if env_corr_type == "delta" && (n != m)

                    continue
                
                else

                    # aD_0 * f(a(n-m)) * L_m
                    Lm = 0.5 * tau * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * get_Lindblad_jump_operator(m, aT, sites)

                    # L_n^\dagger 
                    Ldagn = hermitian_conjugate_mpo(get_Lindblad_jump_operator(n, aT, sites))

                    # 0.5 * dt * aD_0 * f(a(n-m)) * L_m * rho * L_n^\dagger
                    Lmrho = apply(Lm, rho; cutoff = cutoff)
                    rho_final += apply(Lmrho, Ldagn; cutoff = cutoff)

                end
        
            end
        end

    end

    return rho_final

end

function get_entanglement_entropy_mpo(rho, trace_indices, sites; tol = 1e-12)

    N = length(sites) - length(trace_indices)

    tmp = []
    for (idx, element) in enumerate(rho)
        if idx in trace_indices
            push!(tmp, element * delta(dag(sites[idx]'), sites[idx]))
        else
            push!(tmp, element)
        end
    end

    a = contract(tmp)
    a = Array(a, inds(a; :plev => 1)..., inds(a; :plev => 0)...)
    a = reshape(a, 2^N, 2^N)

    evals, _ = eigen(a)

    ee = sum(-real(eval)*log(real(eval)) for eval in evals if real(eval) >= tol)

    return ee

end

function get_entanglement_entropy_matrix(N, rho_m, keep_indices; tol = 1e-12)

    a = partial_trace(rho_m, keep_indices)

    a = reshape(a, 2^(div(N, 2)), 2^(div(N, 2)))

    evals, _ = eigen(a)

    ee2 = sum(-real(eval)*log(real(eval)) for eval in evals if real(eval) >= tol)

    return ee2

end

function check_zeroq(n, N)
    return sum((digits(n, base=2, pad = N).*2).-1) == 0
end

function project_zeroq(M)

    nrow, ncol = size(M)
    n = Int(log2(nrow))
    new_nrow = binomial(n, div(n, 2))
    res = zeros(ComplexF64, new_nrow, new_nrow)

    row_count = 0
    for row in 1:nrow

        if !(check_zeroq(row-1, n))
            continue
        else
            row_count += 1
        end
        # println(row_count, " ", digits(row-1, base = 2, pad = 4))
        col_count = 0
        for col in 1:ncol

            if !(check_zeroq(col-1, n))
                continue
            else
                col_count += 1
                # println("row_count, col_count: ", (row_count, col_count), ", row, col: ", (digits(row-1, base = 2, pad = n), digits(col-1, base = 2, pad = n)))
                res[row_count, col_count] = M[row, col]
            end

        end

    end

    return res

end

function get_Lindblad_reduced_sparse_matrix(N, x, ma, l_0, lambda, aD_0, sigma_over_a, aT, env_corr_type)

    ldim = binomial(N, div(N, 2))^2

    L = spzeros(ldim, ldim)

    eye(n::Int64) = sparse(I, n, n);

    H = get_aH_Hamiltonian_sparse_matrix(N, x, ma, l_0, lambda)
    H_r = project_zeroq(H)
    idnt = eye(2^N)
    idnt_r = project_zeroq(idnt)

    # Unitary part of Lindbladian
    L += -1im * my_kron(H_r, idnt_r) + 1im * my_kron(idnt_r, transpose(H_r)) 

    if aD_0 != 0

        for n in 1:N
            for m in 1:N

                tmp1 = project_zeroq(get_Lindblad_jump_operator_sparse_matrix(N, n, aT))
                tmp2 = project_zeroq(get_Lindblad_jump_operator_sparse_matrix(N, m, aT))

                tmp3 = tmp1' * tmp2 # the dash is the dagger
                
                L += aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * ((my_kron(tmp2, conj(tmp1))) - 0.5*(my_kron(tmp3, idnt_r)) -0.5*(my_kron(idnt_r, transpose(tmp3))))

            end
        end

    end

    return L

end

function get_entanglement_entropy_reduced_matrix(N, rho_m; tol = 1e-12)

    # If you have 4 qubits 1234 it computes the von Neumann entropy by partitioning 
    # the system in half 12 34 so it always assumes you have even number of sites

    dimres = 2^(div(N, 2))
    
    res = zeros(ComplexF64, dimres, dimres)

    zero_q_list = [join(digits(i-1, base = 2, pad = N)) for i in 1:2^N if check_zeroq(i-1, N)] 

    for row in 1:dimres
        for col in 1:dimres

            for trace_idx in 1:dimres 

                bigrow = join(vcat(digits(row-1, base=2, pad = div(N,2)), digits(trace_idx-1, base=2, pad = div(N,2)))) # this is bin(row))bin(trace_idx)
                bigcol = join(vcat(digits(col-1, base=2, pad = div(N,2)), digits(trace_idx-1, base=2, pad = div(N,2)))) # this is bin(col))bin(trace_idx)

                if !(bigrow in zero_q_list) || !(bigcol in zero_q_list)
                    continue
                else
                    bigrow_idx = findfirst(x -> x == bigrow, zero_q_list)
                    bigcol_idx = findfirst(x -> x == bigcol, zero_q_list)
                    # println("r_a = ", row, ", c_a = ", col, " ", bigrow, " ", bigcol, " r_a r_b = ", bigrow_idx, ", c_a c_b = ", bigcol_idx, " value = ", real(rho_m[bigrow_idx, bigcol_idx]))
                    # println(bigrow_idx, " ", bigcol_idx)
                    res[row, col] += rho_m[bigrow_idx, bigcol_idx]
                end

            end
                
        end

    end

    evals, _ = eigen(res)

    ee2 = sum(-real(eval)*log(real(eval)) for eval in evals if real(eval) >= tol)

    return ee2

end

function swap(i, j, N)

    res = sparse(I, 2^N, 2^N)

    idx1 = min(i, j)
    idx2 = max(i, j)

    local_swap = sparse([[1 0 0 0]; [0 0 1 0]; [0 1 0 0]; [0 0 0 1]])
    full_swap(idx) = kron(sparse(I, 2^(idx-1), 2^(idx-1)), kron(local_swap, sparse(I, 2^(N-idx-1), 2^(N-idx-1))))

    for k in idx1:idx2-1

        res *= full_swap(k)

    end

    if idx2-idx1 > 1

        for k in reverse(idx1:idx2-2)

            res *= full_swap(k)

        end

    end

    return res

end

function get_CP_operator_sparse(N)

    x(idx) = get_op(["X"], [idx], N; reverse_flag = false)

    res = sparse(I, 2^N, 2^N)

    for j in 1:div(N, 2)

        res *= x(j)*x(N+1-j)*swap(j, N+1-j, N)

    end

    return res

end

function decimal_to_padded_binary_list(decimal, bit_length)
    binary_list = Int[]

    while decimal > 0 || length(binary_list) < bit_length
        pushfirst!(binary_list, (decimal % 2) + 1)
        decimal = div(decimal, 2)
    end

    # Pad with leading zeros if needed
    while length(binary_list) < bit_length
        pushfirst!(binary_list, 0)
    end

    return binary_list 
end

function mps_to_list(mps)
    
    res = []
    tmp = contract(mps)
    N = length(mps)
    for i in 1:2^N
        # Convert a decimal into a list of integers that represent the number in binary 
        # (instead of 0 and 1 we have 1 and 2 in this list to fit with Julia)
        binary_list = decimal_to_padded_binary_list(i-1, N) 
        push!(res, tmp[binary_list...])
    end
    return res

end

function get_charge_config_sparse(s)

    config = []

    q_n(n) = (get_op(["Z"], [n], N) + (-1)^(n-1)*sparse(I, 2^N, 2^N))*0.5

    for i in 1:N

        push!(config, s'*q_n(i)*s)

    end

    return config

end

function get_product_mps(state, sites)

    N = length(state)
    mps = MPS(sites)
    links = [Index(QN() => 1; tags = "Link,l=$(n)") for n in 1:N-1]

    # Index(QN() => 1, QN("Sz", -2) => 1, QN("Sz", 2) => 1, QN("Sz", 0) => 1, QN("Sz", 0) => 1; tags = join(["Link,l=", string(n)]))

    for n in 1:N

        if n == 1

            s, lr = sites[n], links[n]
            
            mps[n] = ITensor(ComplexF64, s, lr)

            if state[n] == "0"
                mps[n][s => 1, lr => 1] = 1
                # mps[n][s => 2, lr => 1] = 0
            else
                # mps[n][s => 1, lr => 1] = 0
                mps[n][s => 2, lr => 1] = 1
            end
            
        elseif n == N

            s, ll = sites[n], dag(links[n-1])
            
            mps[n] = ITensor(ComplexF64, s, ll)

            if state[n] == "0"
                mps[n][s => 1, ll => 1] = 1
                # mps[n][s => 2, ll => 1] = 0
            else
                # mps[n][s => 1, ll => 1] = 0
                mps[n][s => 2, ll => 1] = 1
            end

        else

            s, ll, lr = sites[n], dag(links[n-1]), links[n]

            mps[n] = ITensor(ComplexF64, s, ll, lr)

            if state[n] == "0"
                mps[n][s => 1, ll => 1, lr => 1] = 1
                # mps[n][s => 2, ll => 1, lr => 1] = 0
            else
                # mps[n][s => 1, ll => 1, lr => 1] = 0
                mps[n][s => 2, ll => 1, lr => 1] = 1
            end

        end
    end

    return mps

end

function get_particle_number_operator_sparse(N)

    op = 0.5*N*sparse(I, 2^N, 2^N)
    for n in 1:N
        op += 0.5*(-1)^(n-1)*get_op(["Z"], [n], N)
    end

    return op

end

function get_particle_number_operator_sparse_old(N, x)

    op = zeros(2^N, 2^N)
    for n in 1:N
        op += (sqrt(x)/(4*N))*get_op(["Z"], [n], N)
    end

    return op

end

function get_charge_config_from_zeroq_density_matrix(N, rho)

    res = []
    for i in 1:N
        op = project_zeroq(0.5*(-1)^(i-1)*sparse(I, 2^N, 2^N) + 0.5*get_op(["Z"], [i], N))
        push!(res, real(tr(rho*op)))
    end

    return res

end

function get_electric_field_from_zeroq_density_matrix(N, rho, l_0)

    res = []
    charge_config = get_charge_config_from_zeroq_density_matrix(N, rho)
    for i in 1:N
        push!(res, l_0 + sum(charge_config[1:i]))
    end

    return res

end

function get_entanglement_entropy_reduced_from_environment(rho; tol = 1e-12)

    evals = eigen(Matrix(rho)).values

    return sum(-real(eval)*log(real(eval)) for eval in evals if real(eval) >= tol)


end

function get_CP_operator(sites)

    N = length(sites)

    # MPO for X gate
    function x_gate(idx)
        opsum = OpSum()
        opsum += "X",idx
        x_mpo = MPO(opsum, sites)
        return x_mpo
    end

    # MPO for swap gate for site, site + 1
    function swap_nearest_neighbour_gate_mpo(idx)
    
        opsum = OpSum()
        opsum += 0.5,"I",1
        opsum += 0.5,"X",idx,"X",idx+1
        opsum += 0.5,"Y",idx,"Y",idx+1
        opsum += 0.5,"Z",idx,"Z",idx+1 
        swap_mpo = MPO(opsum, sites)

        return swap_mpo
        
    end

    # MPO for swap operator between sites i and j
    function swap(i, j)

        idx1 = min(i, j)
        idx2 = max(i, j)

        res = MPO(sites, "Id")

        for k in idx1:idx2-1

            res = apply(swap_nearest_neighbour_gate_mpo(k), res)
    
        end
    
        if idx2-idx1 > 1
    
            for k in reverse(idx1:idx2-2)
    
                res = apply(swap_nearest_neighbour_gate_mpo(k), res)
    
            end
    
        end

        return res

    end

    final_res = MPO(sites, "Id")

    for j in 1:div(N, 2)

        final_res = apply(x_gate(N+1-j), final_res)
        final_res = apply(x_gate(j), final_res)
        final_res = apply(swap(j, N+1-j), final_res)

    end

    return final_res

end

function get_dirac_vacuum_zeroq_density_matrix_sparse(N)

    state = join([isodd(n) ? "0" : "1" for n = 1:N])
    decimal_number = parse(Int, state, base=2) + 1
    rho = zeros(2^N, 2^N)
    rho[decimal_number, decimal_number] = 1 
    return project_zeroq(rho)

end
