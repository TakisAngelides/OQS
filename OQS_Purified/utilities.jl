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

    if aD_0 != 0.0
        mpo2 = get_Lindblad_dissipative_part(aD_0, sigma_over_a, env_corr_type, aT, sites)
        res = mpo1 + mpo2
    else
        res = mpo1
    end

    return res

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

function apply_odd!(odd, mps, cutoff, maxdim)

    l = length(mps)

    for (n_idx, n) in enumerate(1:4:(l-3)) # n is the left most site the gate acts on, the gates are always 4 site in span

        gate = odd[n_idx]

        t = noprime(gate*prod(mps[n:n+3]))

        # Separate the t tensor into individual MPS site tensors
        for idx in n:n+2
            t_indices = inds(t)
            U_indices = filter(i -> hastags(i, "Site,n=$(idx)") || hastags(i, "Link,l=$(idx-1)"), t_indices)
            U, S, V = ITensors.svd(t, U_indices; cutoff = cutoff, maxdim = maxdim, lefttags = "Link,l=$(idx)", righttags = "Link,l=$(idx)")
            t = S*V
            mps[idx] = U
        end
        mps[n+3] = t

        if n+4 <= l
            # Extra SVD as required by ATD DMRG
            t = mps[n+3]*mps[n+4]
            mps[n+3], S, V = ITensors.svd(t, uniqueinds(t, mps[n+4]); cutoff = cutoff, maxdim = maxdim, lefttags = "Link,l=$(n+3)", righttags = "Link,l=$(n+3)")
            mps[n+4] = S*V
        end

    end

end

function apply_even!(even, mps, cutoff, maxdim)

    l = length(mps)
    
    t = mps[1]*mps[2]
    mps[1], mps[2] = ITensors.qr(t, uniqueinds(t, mps[2]); tags = "Link,l=$(1)")

    t = mps[2]*mps[3]
    mps[2], mps[3] = ITensors.qr(t, uniqueinds(t, mps[3]); tags = "Link,l=$(2)")

    for (n_idx, n) in enumerate(3:4:(l-3)) # n is the left most site the gate acts on, the gates are always 4 site in span

        gate = even[n_idx]

        t = noprime(gate*prod(mps[n:n+3]))

        # Separate the t tensor into individual MPS site tensors
        for idx in n:n+2
            t_indices = inds(t)
            U_indices = filter(i -> hastags(i, "Site,n=$(idx)") || hastags(i, "Link,l=$(idx-1)"), t_indices)
            U, S, V = ITensors.svd(t, U_indices; cutoff = cutoff, maxdim = maxdim, lefttags = "Link,l=$(idx)", righttags = "Link,l=$(idx)")
            t = S*V
            mps[idx] = U
        end
        mps[n+3] = t

        if n_idx != length(3:4:(l-3))
            # Extra SVD as required by ATD DMRG
            t = mps[n+3]*mps[n+4]
            mps[n+3], S, V = ITensors.svd(t, uniqueinds(t, mps[n+4]); cutoff = cutoff, maxdim = maxdim, lefttags = "Link,l=$(n+3)", righttags = "Link,l=$(n+3)")
            mps[n+4] = S*V
        end

    end

end

function get_entanglement_entropy(psi, site, tol = 1e-12)
    
    ITensors.orthogonalize!(psi, site)
    if site == 1
        U,S,V = svd(psi[site], siteind(psi, site); cutoff = 0)
    else
        U,S,V = svd(psi[site], (linkind(psi, site-1), siteind(psi, site)); cutoff = 0)
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

    N = length(psi)

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

function lowest_eval_mpo(mpo; tol = 1e-9)

    n = length(mpo)
    dmrg_tol = tol
    cutoff = tol
    sites = dag(reduce(vcat, siteinds(mpo; :plev => 0)))
    state = [isodd(i) ? "0" : "1" for i = 1:n]
    psi0 = randomMPS(sites, state)
    nsweeps = 100
    dmrg_energy, _ = dmrg(mpo, psi0; nsweeps, cutoff, outputlevel=1, ishermitian = false)

    return dmrg_energy

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
        H += 0.5*get_op(["S-", "S+"], [n, n+1], N)
        H += 0.5*get_op(["S+", "S-"], [n, n+1], N)
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

function environment_correlator(type, n, m, aD, inputs)

    if type == "constant"
        return aD
    elseif type == "delta"
        if n == m
            return aD
        else
            return 0.0
        end
    else # gaussian case
        sigma_over_a = inputs["sigma_over_a"]
        return aD*exp(-0.5*(1/sigma_over_a)^2*(n-m)^2)
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

    op_dict = Dict("X" => sparse([0 1; 1 0]), "Y" => sparse([0 -1im; 1im 0]), "Z" => sparse([1 0; 0 -1]), "S+" => sparse([0 1; 0 0]), "S-" => sparse([0 0; 1 0]))
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

    """
    Gives the operator L_n_dagger * L_m
    """

    N = length(sites)
    
    opsum = OpSum()

    opsum += 0.25*(-1)^(n+m),"Z",n,"Z",m

    opsum += 0.25*(-1)^(n+m),"Z",n

    opsum += 0.25*(-1)^(n+m),"Z",m

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

    if (n + m) % 2 != 0 # this is the case of the first part of the dissipator with aL_2m-1 * aL_2n_dagger where the two act on left and right (top and bottom) respectively

        if (n != 1)
            opsum += (-(-1)^(n + m)/(16*aT)),"S-",n-1,"S+",n
            opsum += ((-1)^(n + m)/(16*aT)),"S+",n-1,"S-",n
        end

        if (n != N)
            opsum += ((-1)^(n + m)/(16*aT)),"S-",n,"S+",n+1
            opsum += (-(-1)^(n + m)/(16*aT)),"S+",n,"S-",n+1
        end
    
        if (m != 1)
            opsum += ((-1)^(n + m)/(16*aT)),"S-",m-1,"S+",m
            opsum += (-(-1)^(n + m)/(16*aT)),"S+",m-1,"S-",m
        end
        
        if (m != N)
            opsum += (-(-1)^(n + m)/(16*aT)),"S-",m,"S+",m+1
            opsum += ((-1)^(n + m)/(16*aT)),"S+",m,"S-",m+1
        end

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

    mpo = environment_correlator(env_corr_type, 1, 1, aD_0, sigma_over_a) * get_LdagL(1, 1, aT, sites)

    for n=1:N
        for m=1:N

            if !((n == 1) && (m == 1))
                mpo += environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * get_LdagL(n, m, aT, sites)
            end

        end
    end

    return mpo

end

function get_double_size_Lindblad_operator(N, sites, x, ma, l_0, lambda, aD_0, sigma_over_a, aT, env_corr_type)

    # TODO: as it is now this function does not work when sites have conserved_qns, the XY need to change to S-S+ etc

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

    h1 = MPO(opsum, sites)

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
                h1 = add(h1, apply(h2, h3; cutoff = 0); alg = "directsum")
                
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
                h1 = add(h1, MPO(opsum, sites; cutoff = 0); alg = "directsum")

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
        res += (-(-1)^m/(8*aT))*get_op(["S+", "S-"], [m-1, m], N)
        res += ((-1)^m/(8*aT))*get_op(["S-", "S+"], [m-1, m], N)
    end
    
    if m != N
        res += ((-1)^m/(8*aT))*get_op(["S+", "S-"], [m, m+1], N)
        res += (-(-1)^m/(8*aT))*get_op(["S-", "S+"], [m, m+1], N)
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
    L += -1im * my_kron(H, eye(2^N)) # + 1im * my_kron(eye(2^N), transpose(H)) 

    for n in 1:N
        for m in 1:N

            tmp1 = get_Lindblad_jump_operator_sparse_matrix(N, n, aT)
            tmp2 = get_Lindblad_jump_operator_sparse_matrix(N, m, aT)

            tmp3 = tmp1' * tmp2 # the dash is the dagger
            
            L += environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * ((my_kron(tmp2, transpose(tmp1'))) - 0.5*(my_kron(tmp3, eye(2^N))) -0.5*(my_kron(eye(2^N), transpose(tmp3))))

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

function apply_taylor_part(rho, tau, sites, x, l_0, ma, aD_0, sigma_over_a, env_corr_type, aT; cutoff = 0, max_rho_D = 500)

    H_T = get_H_taylor(sites, x, l_0, ma)

    tmp = -0.5*1im*tau*H_T

    # rho + rho * idt/2 H_T * rho - idt/2 H_T * rho
    rho_final = rho + apply(rho, hermitian_conjugate_mpo(tmp); cutoff = cutoff, maxdim = max_rho_D) + apply(tmp, rho; cutoff = cutoff, maxdim = max_rho_D)

    # second order term in the taylor expansion only for the Hamiltonian part
    # rho_final += (-tau^2/8)*apply(H_T, apply(H_T, rho; cutoff = cutoff, maxdim = max_rho_D); cutoff = cutoff, maxdim = max_rho_D) + (-tau^2/8)*apply(rho, apply(H_T, H_T; cutoff = cutoff, maxdim = max_rho_D); cutoff = cutoff, maxdim = max_rho_D) + (tau^2/4)*apply(H_T, apply(rho, H_T; cutoff = cutoff, maxdim = max_rho_D); cutoff = cutoff, maxdim = max_rho_D)

    if aD_0 != 0

        # - 0.5 * dt * 0.5 * L_n^\dagger L_m
        LdagnLm = -0.5 * tau * 0.5 * get_Lindblad_dissipative_part(aD_0, sigma_over_a, env_corr_type, aT, sites)

        # sum over n and m: - 0.5 * dt/2 * L_n^\dagger L_m * rho
        rho_final += apply(LdagnLm, rho; cutoff = cutoff, maxdim = max_rho_D) 
        
        # sum over n and m: - 0.5 * dt/2 * rho * L_n^\dagger L_m 
        rho_final += apply(rho, LdagnLm; cutoff = cutoff, maxdim = max_rho_D) 

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
                    Lmrho = apply(Lm, rho; cutoff = cutoff, maxdim = max_rho_D)
                    rho_final += apply(Lmrho, Ldagn; cutoff = cutoff, maxdim = max_rho_D)

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

function get_Lindblad_reduced_sparse_matrix(N, x, ma, l_0, lambda, aD_0, sigma_over_a, aT, env_corr_type; edges_only = false)

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

        if edges_only == false

            for n in 1:N
                for m in 1:N

                    tmp1 = project_zeroq(get_Lindblad_jump_operator_sparse_matrix(N, n, aT))
                    tmp2 = project_zeroq(get_Lindblad_jump_operator_sparse_matrix(N, m, aT))

                    tmp3 = tmp1' * tmp2 # the dash is the dagger
                    
                    L += environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (my_kron(tmp2, conj(tmp1)) - 0.5*my_kron(tmp3, idnt_r) -0.5*my_kron(idnt_r, transpose(tmp3)))

                end
            end

        else

            for n in [1, N]
                for m in [1, N]

                    tmp1 = project_zeroq(get_Lindblad_jump_operator_sparse_matrix(N, n, aT))
                    tmp2 = project_zeroq(get_Lindblad_jump_operator_sparse_matrix(N, m, aT))

                    tmp3 = tmp1' * tmp2 # the dash is the dagger
                    
                    L += environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (my_kron(tmp2, conj(tmp1)) - 0.5*my_kron(tmp3, idnt_r) -0.5*my_kron(idnt_r, transpose(tmp3)))

                end
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

function get_z_config_from_zeroq_density_matrix(N, rho)

    res = []
    for i in 1:N
        op = project_zeroq(get_op(["Z"], [i], N))
        push!(res, real(tr(rho*op)))
    end

    return res

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
    for i in 1:N-1
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

function get_dirac_vacuum_with_string_zeroq_density_matrix_sparse(N, q_sites)

    state = [isodd(n) ? "0" : "1" for n = 1:N]
    for i in 1:N
        if i in q_sites
            if isodd(i)
                state[i] = "1"
            else
                state[i] = "0"
            end
        end
    end
    state = join(state)
    decimal_number = parse(Int, state, base=2) + 1
    rho = zeros(2^N, 2^N)
    rho[decimal_number, decimal_number] = 1 
    return project_zeroq(rho)

end

function get_dirac_vacuum_density_matrix(sites)

    N = length(sites)
    state = [isodd(n) ? "1" : "0" for n = 1:N] 
    mps = MPS(sites, state)
   
    return outer(mps', mps)

end

function get_dirac_vacuum_mps(sites; flip_sites = [])

    N = length(sites)
    state = [isodd(n) ? "1" : "0" for n = 1:N]
    state = []
    for n in 1:N
        if isodd(n)
            if n in flip_sites
                push!(state, "0")
            else
                push!(state, "1")
            end
        else
            if n in flip_sites
                push!(state, "1")
            else
                push!(state, "0") 
            end
        end
    end
    mps = MPS(sites, state)
   
    return mps

end

function get_total_charge_reduced_operator_sparse(N)

    res = sparse(I, binomial(N, div(N, 2)), binomial(N, div(N, 2)))
    for i in 1:N
        res += project_zeroq(0.5*(-1)^(i-1)*sparse(I, 2^N, 2^N) + 0.5*get_op(["Z"], [i], N))
    end

    return res

end

function get_delta_tensor(idx)

    s1, s2 = dag(idx), dag(idx')
    d = ITensor(diagm(0 => ones(Int, dim(s1)))[:, end:-1:1], s1, s2)

    return d

end

function reverse_dir(t, idx)
    return noprime(t*get_delta_tensor(idx); :tags => tags(idx))
end

function rho_vec_to_mps(rho_vec)

    N = length(rho_vec)

    mps = MPS(2*N)

    for i in 1:N

        # this object has 2 physical legs and we want to svd between them to seperate them into two mps tensors
        M = rho_vec[i]

        # picks the left link and the dashed site indices for the U of SVD
        left_inds_M = (inds(M; :tags => "Link,l=$(i-1)"), inds(M; :plev => 1))

        U, S, V = ITensors.svd(M, left_inds_M; leftdir = ITensors.In, rightdir = ITensors.In)
        V = S*V

        # Fix U indices
        if i != 1
            replacetags!(U, "Link,l=$(i-1)", "Link,l=$(2*i-2)") # change tag of left link to mps convention
        end
        replacetags!(U, "Link,u", "Link,l=$(2*i-1)") # change tag of svd link to mps convention
        replacetags!(U, "S=1/2,Site,n=$(i)", "S=1/2,Site,n=$(2*i-1)") # change the label of the second physical index from i to 2*i-1
        mps[2*i-1] = U

        # Fix V indices
        if i != N # if i == N there is no right link
            replacetags!(V, "Link,l=$(i)", "Link,l=$(2*i)") # same as the line above but for the right link
        end
        replacetags!(V, "Link,u", "Link,l=$(2*i-1)") # change tag of svd link to mps convention (here it is Link,u because of V = S*V above)
        replacetags!(V, "S=1/2,Site,n=$(i)", "S=1/2,Site,n=$(2*i)") # change the label of the second physical index from i to 2*i
        # V = reverse_dir(V, inds(V; :tags => "Site")[1]) # reverse the physical leg direction from in to out - this is mandatory for constructing MPO with autoMPO
        mps[2*i] = V

    end

    return noprime(mps)

end

function compare(z1::Complex, z2::Complex)
    real(z1) > real(z2)
end

function get_double_aH_Hamiltonian(sites, x, l_0, ma, lambda, side)

    """
    This gives aH Hamiltonian

    side specifies "left" or "right" to imply H tensor product I or vice versa

    """

    N = div(length(sites), 2)

    opsum = OpSum()

    for n in 1:N-1

        if side == "left"
            n_idx = 2*n-1
        else
            n_idx = 2*n 
        end
        
        for m in n+1:N

            if side == "left"
                m_idx = 2*m-1
            else
                m_idx = 2*m
            end
            
            # Long range ZZ interaction term
            opsum += 0.25*(1/x)*(N-m+lambda),"Z",n_idx,"Z",m_idx

        end

        # Kinetic term
        opsum += 0.5,"S+",n_idx,"S-",n_idx+2
        opsum += 0.5,"S-",n_idx,"S+",n_idx+2

        opsum += (1/x)*(N/8 - 0.25*ceil((n-1)/2) + l_0*(N-n)/2),"Z",n_idx
        
        opsum += (0.5*ma*(-1)^(n-1)),"Z",n_idx

    end

    if side == "left"
        opsum += (0.5*ma*(-1)^(N-1)),"Z",2*N-1
        opsum += ((l_0^2)*(N-1)/(2*x) + (l_0*N)/(4*x) + (N^2)/(16*x) + (lambda*N/(8*x))),"Id",1
    else
        opsum += (0.5*ma*(-1)^(N-1)),"Z",2*N
        opsum += ((l_0^2)*(N-1)/(2*x) + (l_0*N)/(4*x) + (N^2)/(16*x) + (lambda*N/(8*x))),"Id",2
    end

    return opsum

end

function get_double_aH_Hamiltonian_taylor_part(sites, x, l_0, ma, lambda, side)

    """
    This gives aH Hamiltonian

    side specifies "left" or "right" to imply H tensor product I or vice versa

    """

    N = div(length(sites), 2)

    opsum = OpSum()

    for n in 1:N-1

        if side == "left"
            n_idx = 2*n-1
        else
            n_idx = 2*n 
        end
        
        for m in n+1:N

            if side == "left"
                m_idx = 2*m-1
            else
                m_idx = 2*m
            end
            
            # Long range ZZ interaction term
            opsum += 0.25*(1/x)*(N-m+lambda),"Z",n_idx,"Z",m_idx

        end

        opsum += (1/x)*(N/8 - 0.25*ceil((n-1)/2) + l_0*(N-n)/2),"Z",n_idx
        
        opsum += (0.5*ma*(-1)^(n-1)),"Z",n_idx

    end

    if side == "left"
        opsum += (0.5*ma*(-1)^(N-1)),"Z",2*N-1
        opsum += ((l_0^2)*(N-1)/(2*x) + (l_0*N)/(4*x) + (N^2)/(16*x) + (lambda*N/(8*x))),"Id",1
    else
        opsum += (0.5*ma*(-1)^(N-1)),"Z",2*N
        opsum += ((l_0^2)*(N-1)/(2*x) + (l_0*N)/(4*x) + (N^2)/(16*x) + (lambda*N/(8*x))),"Id",2
    end

    return opsum

end

function get_Lindblad_dissipator(sites, x, l_0, ma, lambda, aT, aD, env_corr_type, sigma_over_a)

    N = div(length(sites), 2)
    n, m = 1, 1
    res = environment_correlator(env_corr_type, n, m, aD, sigma_over_a) * (-0.5*get_aLndag_aLm(2*n-1, 2*m-1, aT, sites, "left") - 0.5*get_aLndag_aLm(2*n, 2*m, aT, sites, "right"))

    for n in 1:N
        for m in 1:N
            if (n == 1) && (m == 1)
                continue
            end
            res += environment_correlator(env_corr_type, n, m, aD, sigma_over_a) * ( get_aLm_aLndag(2*n, 2*m-1, aT, sites) - 0.5 * get_aLndag_aLm(2*n-1, 2*m-1, aT, sites, "left") - 0.5 * get_aLndag_aLm(2*n, 2*m, aT, sites, "right") )
        end
    end

    return MPO(res, sites)

end

function get_Lindblad(sites, x, l_0, ma, lambda, aT, aD, env_corr_type, sigma_over_a)

    N = div(length(sites), 2)
    res = -1im*get_double_aH_Hamiltonian(sites, x, l_0, ma, lambda, "left")
    res += 1im*get_double_aH_Hamiltonian(sites, x, l_0, ma, lambda, "right")
    if aD != 0
        for n in 1:N
            for m in 1:N
                res += environment_correlator(env_corr_type, n, m, aD, sigma_over_a) * ( get_aLm_aLndag(2*n, 2*m-1, aT, sites) - 0.5 * get_aLndag_aLm(2*n-1, 2*m-1, aT, sites, "left") - 0.5 * get_aLndag_aLm(2*n, 2*m, aT, sites, "right") )
            end
        end
    end

    return MPO(res, sites)

end

function trace_mps(mps)

    sites = siteinds(mps)
    l = length(mps)

    i = 1
    left, right = 2*i-1, 2*i
    res = mps[left]*mps[right]*dag(delta(sites[left], sites[right]))

    for i in 2:div(l,2)
        left, right = 2*i-1, 2*i
        res *= mps[left]*mps[right]*dag(delta(sites[left], sites[right]))
    end

    return res[1]

end

function measure_z(mps, site)

    sites = siteinds(mps)
    l = length(mps)

    z = op("Z", sites[site])
    mps = apply(z, mps)

    i = 1
    left, right = 2*i-1, 2*i
    res = mps[left]*mps[right]*dag(delta(sites[left], sites[right]))

    for i in 2:div(l,2)
        left, right = 2*i-1, 2*i
        res *= mps[left]*mps[right]*dag(delta(sites[left], sites[right]))
    end

    return res[1]

end

function measure_z_config(mps; left = true)

    n = length(mps)
    z_config = []
    if left
        for site in 1:2:n
            push!(z_config, measure_z(mps, site))
        end
    else
        for site in 2:2:n
            push!(z_config, measure_z(mps, site))
        end
    end

    return z_config

end

function measure_particle_number(mps; normalize = false)

    n = length(mps)
    N = div(n, 2)
    z_config = []
    for site in 1:2:n-1
        push!(z_config, measure_z(mps, site))
    end

    res = 0
    for (idx, element) in enumerate(z_config)

        res += 0.5*(-1)^(idx-1)*element

    end

    if normalize
        return 0.5*N + res/trace_mps(mps)
    else
        return 0.5*N + res
    end

end

function get_aLm_aLndag(n, m, aT, sites)

    """
    Gives the operator aLm tensor product aLndagger so aLm acts on the left side system and aLndagger on the right
    """

    N = div(length(sites), 2)
    
    n_phys, m_phys = div(n, 2), div(m + 1, 2)

    opsum = OpSum()

    opsum += 0.25*(-1)^(n_phys+m_phys),"Z",m,"Z",n

    opsum += 0.25*(-1)^(n_phys+m_phys),"Z",n

    opsum += 0.25*(-1)^(n_phys+m_phys),"Z",m

    opsum += 0.25*(-1)^(n_phys+m_phys),"Id",1

    if (n_phys != 1)
        opsum += (-(-1)^(n_phys + m_phys)/(16*aT)),"S-",n-2,"S+",n
        opsum += ((-1)^(n_phys + m_phys)/(16*aT)),"S+",n-2,"S-",n
    end

    if (n_phys != N)
        opsum += ((-1)^(n_phys + m_phys)/(16*aT)),"S-",n,"S+",n+2
        opsum += (-(-1)^(n_phys + m_phys)/(16*aT)),"S+",n,"S-",n+2
    end

    if (m_phys != 1)
        opsum += ((-1)^(n_phys + m_phys)/(16*aT)),"S-",m-2,"S+",m
        opsum += (-(-1)^(n_phys + m_phys)/(16*aT)),"S+",m-2,"S-",m
    end
    
    if (m_phys != N)
        opsum += (-(-1)^(n_phys + m_phys)/(16*aT)),"S-",m,"S+",m+2
        opsum += ((-1)^(n_phys + m_phys)/(16*aT)),"S+",m,"S-",m+2
    end

    if (m_phys != 1)
        opsum += ((-1)^(n_phys + m_phys)/(16*aT)),"S-",m-2,"S+",m,"Z",n
        opsum += (-(-1)^(n_phys + m_phys)/(16*aT)),"S+",m-2,"S-",m,"Z",n
    end
    
    if (m_phys != N)
        opsum += (-(-1)^(n_phys + m_phys)/(16*aT)),"S-",m,"S+",m+2,"Z",n
        opsum += ((-1)^(n_phys + m_phys)/(16*aT)),"S+",m,"S-",m+2,"Z",n
    end

    if (n_phys != 1)
        opsum += (-(-1)^(n_phys + m_phys)/(16*aT)),"Z",m,"S-",n-2,"S+",n
        opsum += ((-1)^(n_phys + m_phys)/(16*aT)),"Z",m,"S+",n-2,"S-",n
    end
    
    if (n_phys != N)
        opsum += ((-1)^(n_phys + m_phys)/(16*aT)),"Z",m,"S-",n,"S+",n+2
        opsum += (-(-1)^(n_phys + m_phys)/(16*aT)),"Z",m,"S+",n,"S-",n+2
    end

    if (n_phys != 1) && (m_phys != 1)
        opsum += (-(-1)^(n_phys + m_phys)/(64*aT^2)),"S+",m-2,"S-",m,"S+",n-2,"S-",n
        opsum += ((-1)^(n_phys + m_phys)/(64*aT^2)),"S+",m-2,"S-",m,"S-",n-2,"S+",n
        opsum += ((-1)^(n_phys + m_phys)/(64*aT^2)),"S-",m-2,"S+",m,"S+",n-2,"S-",n
        opsum += (-(-1)^(n_phys + m_phys)/(64*aT^2)),"S-",m-2,"S+",m,"S-",n-2,"S+",n
    end

    if (m_phys != 1) && (n_phys != N)
        opsum += (-(-1)^(n_phys + m_phys)/(64*aT^2)),"S+",m-2,"S-",m,"S-",n,"S+",n+2
        opsum += ((-1)^(n_phys + m_phys)/(64*aT^2)),"S+",m-2,"S-",m,"S+",n,"S-",n+2
        opsum += (-(-1)^(n_phys + m_phys)/(64*aT^2)),"S-",m-2,"S+",m,"S+",n,"S-",n+2
        opsum += ((-1)^(n_phys + m_phys)/(64*aT^2)),"S-",m-2,"S+",m,"S-",n,"S+",n+2
    end

    if (m_phys != N) && (n_phys != 1)
        opsum += (-(-1)^(n_phys + m_phys)/(64*aT^2)),"S+",m,"S-",m+2,"S-",n-2,"S+",n
        opsum += ((-1)^(n_phys + m_phys)/(64*aT^2)),"S+",m,"S-",m+2,"S+",n-2,"S-",n
        opsum += (-(-1)^(n_phys + m_phys)/(64*aT^2)),"S-",m,"S+",m+2,"S+",n-2,"S-",n
        opsum += ((-1)^(n_phys + m_phys)/(64*aT^2)),"S-",m,"S+",m+2,"S-",n-2,"S+",n
    end
    
    if (n_phys != N) && (m_phys != N)
        opsum += (-(-1)^(n_phys + m_phys)/(64*aT^2)),"S+",m,"S-",m+2,"S+",n,"S-",n+2
        opsum += ((-1)^(n_phys + m_phys)/(64*aT^2)),"S+",m,"S-",m+2,"S-",n,"S+",n+2
        opsum += (-(-1)^(n_phys + m_phys)/(64*aT^2)),"S-",m,"S+",m+2,"S-",n,"S+",n+2
        opsum += ((-1)^(n_phys + m_phys)/(64*aT^2)),"S-",m,"S+",m+2,"S+",n,"S-",n+2
    end

    return opsum
    
end

function get_aLndag_aLm(n, m, aT, sites, side)

    """
    Gives the operator aL_n_dagger tensor product aL_m can be either acting all on the left or all on the right
    """

    N = div(length(sites), 2)

    if side == "left"
        n_phys, m_phys = div(n + 1, 2), div(m + 1, 2)
    else
        n_phys, m_phys = div(n, 2), div(m, 2)
    end
    
    opsum = OpSum()

    opsum += 0.25*(-1)^(n_phys + m_phys),"Z",n,"Z",m

    opsum += 0.25*(-1)^(n_phys + m_phys),"Z",n

    opsum += 0.25*(-1)^(n_phys + m_phys),"Z",m

    if side == "left"
        opsum += 0.25*(-1)^(n_phys + m_phys),"Id",1
    else
        opsum += 0.25*(-1)^(n_phys + m_phys),"Id",2
    end

    if (n_phys != 1)
        opsum += (-(-1)^(n_phys + m_phys)/(16*aT)),"S-",n-2,"S+",n
        opsum += ((-1)^(n_phys + m_phys)/(16*aT)),"S+",n-2,"S-",n
    end

    if (n_phys != N)
        opsum += ((-1)^(n_phys + m_phys)/(16*aT)),"S-",n,"S+",n+2
        opsum += (-(-1)^(n_phys + m_phys)/(16*aT)),"S+",n,"S-",n+2
    end

    if (m_phys != 1)
        opsum += ((-1)^(n_phys + m_phys)/(16*aT)),"S-",m-2,"S+",m
        opsum += (-(-1)^(n_phys + m_phys)/(16*aT)),"S+",m-2,"S-",m
    end
    
    if (m_phys != N)
        opsum += (-(-1)^(n_phys + m_phys)/(16*aT)),"S-",m,"S+",m+2
        opsum += ((-1)^(n_phys + m_phys)/(16*aT)),"S+",m,"S-",m+2
    end

    if (n_phys != 1)
        opsum += (-(-1)^(n_phys + m_phys)/(16*aT)),"S-",n-2,"S+",n,"Z",m
        opsum += ((-1)^(n_phys + m_phys)/(16*aT)),"S+",n-2,"S-",n,"Z",m
    end
    
    if (n_phys != N)
        opsum += ((-1)^(n_phys + m_phys)/(16*aT)),"S-",n,"S+",n+2,"Z",m
        opsum += (-(-1)^(n_phys + m_phys)/(16*aT)),"S+",n,"S-",n+2,"Z",m
    end

    if (m_phys != 1)
        opsum += ((-1)^(n_phys + m_phys)/(16*aT)),"Z",n,"S-",m-2,"S+",m
        opsum += (-(-1)^(n_phys + m_phys)/(16*aT)),"Z",n,"S+",m-2,"S-",m
    end
    
    if (m_phys != N)
        opsum += (-(-1)^(n_phys + m_phys)/(16*aT)),"Z",n,"S-",m,"S+",m+2
        opsum += ((-1)^(n_phys + m_phys)/(16*aT)),"Z",n,"S+",m,"S-",m+2
    end

    if (n_phys != 1) && (m_phys != 1)
        opsum += (-(-1)^(n_phys + m_phys)/(64*aT^2)),"S+",n-2,"S-",n,"S+",m-2,"S-",m
        opsum += ((-1)^(n_phys + m_phys)/(64*aT^2)),"S+",n-2,"S-",n,"S-",m-2,"S+",m
        opsum += ((-1)^(n_phys + m_phys)/(64*aT^2)),"S-",n-2,"S+",n,"S+",m-2,"S-",m
        opsum += (-(-1)^(n_phys + m_phys)/(64*aT^2)),"S-",n-2,"S+",n,"S-",m-2,"S+",m
    end

    if (n_phys != 1) && (m_phys != N)
        opsum += (-(-1)^(n_phys + m_phys)/(64*aT^2)),"S+",n-2,"S-",n,"S-",m,"S+",m+2
        opsum += ((-1)^(n_phys + m_phys)/(64*aT^2)),"S+",n-2,"S-",n,"S+",m,"S-",m+2
        opsum += (-(-1)^(n_phys + m_phys)/(64*aT^2)),"S-",n-2,"S+",n,"S+",m,"S-",m+2
        opsum += ((-1)^(n_phys + m_phys)/(64*aT^2)),"S-",n-2,"S+",n,"S-",m,"S+",m+2
    end

    if (n_phys != N) && (m_phys != 1)
        opsum += (-(-1)^(n_phys + m_phys)/(64*aT^2)),"S+",n,"S-",n+2,"S-",m-2,"S+",m
        opsum += ((-1)^(n_phys + m_phys)/(64*aT^2)),"S+",n,"S-",n+2,"S+",m-2,"S-",m
        opsum += (-(-1)^(n_phys + m_phys)/(64*aT^2)),"S-",n,"S+",n+2,"S+",m-2,"S-",m
        opsum += ((-1)^(n_phys + m_phys)/(64*aT^2)),"S-",n,"S+",n+2,"S-",m-2,"S+",m
    end
    
    if (n_phys != N) && (m_phys != N)
        opsum += (-(-1)^(n_phys + m_phys)/(64*aT^2)),"S+",n,"S-",n+2,"S+",m,"S-",m+2
        opsum += ((-1)^(n_phys + m_phys)/(64*aT^2)),"S+",n,"S-",n+2,"S-",m,"S+",m+2
        opsum += (-(-1)^(n_phys + m_phys)/(64*aT^2)),"S-",n,"S+",n+2,"S-",m,"S+",m+2
        opsum += ((-1)^(n_phys + m_phys)/(64*aT^2)),"S-",n,"S+",n+2,"S+",m,"S-",m+2
    end

    return opsum

end

function get_Lindblad_Hamiltonian_part(sites, x, l_0, ma, lambda, aT, aD, env_corr_type, sigma_over_a)

    res = -1im*get_double_aH_Hamiltonian(sites, x, l_0, ma, lambda, "left")
    res += 1im*get_double_aH_Hamiltonian(sites, x, l_0, ma, lambda, "right")
    
    return res

end

function get_mpo_taylor_expansion(mpo, order, cutoff, sites)

    """
    Returns the taylor expansion to the given input order of the mpo as 1 + mpo + mpo * mpo / factorial(2) + mpo * mpo * mpo / factorial(3) + ... etc
    """

    l = [mpo]
    for i in 2:order
        push!(l, apply(l[end], mpo/i; cutoff = cutoff))
    end

    tmp1 = sum(l)
    tmp2 = MPO(sites, "Id")
    
    for i in 2:2:length(tmp2)
        tmp2[i] = swapprime(tmp2[i], 0, 1; :tags => "Site")
    end

    return add(tmp1, tmp2; cutoff = 0)

end

function get_Lindblad_taylor_part(sites, x, l_0, ma, lambda, aT, aD, env_corr_type, sigma_over_a)

    N = div(length(sites), 2)
    res = -1im*get_double_aH_Hamiltonian_taylor_part(sites, x, l_0, ma, lambda, "left")
    res += 1im*get_double_aH_Hamiltonian_taylor_part(sites, x, l_0, ma, lambda, "right")
    if aD != 0
        for n in 1:N
            for m in 1:N
                res += environment_correlator(env_corr_type, n, m, aD, sigma_over_a) * ( get_aLm_aLndag(2*n, 2*m-1, aT, sites) - 0.5 * get_aLndag_aLm(2*n-1, 2*m-1, aT, sites, "left") - 0.5 * get_aLndag_aLm(2*n, 2*m, aT, sites, "right") )
            end
        end
    end

    return MPO(res, sites)

end

function get_odd(sites, a, gate_lists)

    """
    sites are the site indices of the MPS to which the odd gates will act on

    a is usually the time step or the time step divided by 2

    gate_lists is a list of lists specifically nn_odd from the function get_odd_even_taylor_groups in which nn_odd is explained
    """

    gates = []

    l = length(sites)

    for (n_idx, n) in enumerate(1:4:(l-3))

        group_opsum = gate_lists[n_idx] # an opsum with terms lying in the range n:n+3, the actual operators in this opsum always have site numbers from 1 to 4
        group_opsum_indices = sites[n:n+3] # these are the actual indices that will be given to the ITensor below to apply the gate to the appropriate sites

        gate = ITensor(group_opsum, group_opsum_indices)
        Gj = exp(a * gate)

        # Take the transpose on the even sites which correspond to legs which would have been bottom legs on an MPO
        swapprime!(Gj, 0, 1; :tags => "n=$(n+1)")
        swapprime!(Gj, 0, 1; :tags => "n=$(n+3)")

        push!(gates, Gj)

    end

    return gates

end

function get_even(sites, a, gate_lists)

    """
    sites are the site indices of the MPS to which the even gates will act on

    a is usually the time step or the time step divided by 2

    gate_lists is a list of lists specifically nn_even from the function get_odd_even_taylor_groups in which nn_even is explained
    """

    gates = []

    l = length(sites)

    for (n_idx, n) in enumerate(3:4:(l-3))

        group_opsum = gate_lists[n_idx] # an opsum with terms lying in the range n:n+3, the actual operators in this opsum always have site numbers from 1 to 4
        group_opsum_indices = sites[n:n+3] # these are the actual indices that will be given to the ITensor below to apply the gate to the appropriate sites

        gate = ITensor(group_opsum, group_opsum_indices)
        Gj = exp(a * gate)

        # Take the transpose on the even sites which correspond to legs which would have been bottom legs on an MPO
        swapprime!(Gj, 0, 1; :tags => "n=$(n+1)")
        swapprime!(Gj, 0, 1; :tags => "n=$(n+3)")

        push!(gates, Gj)

    end

    return gates
end

function get_Lindblad_opsum(sites, x, l_0, ma, lambda, aT, aD, env_corr_type, sigma_over_a)

    N = div(length(sites), 2)
    res = -1im*get_double_aH_Hamiltonian(sites, x, l_0, ma, lambda, "left")
    res += 1im*get_double_aH_Hamiltonian(sites, x, l_0, ma, lambda, "right")
    if aD != 0
        for n in 1:N
            for m in 1:N
                res += environment_correlator(env_corr_type, n, m, aD, sigma_over_a) * ( get_aLm_aLndag(2*n, 2*m-1, aT, sites) - 0.5 * get_aLndag_aLm(2*n-1, 2*m-1, aT, sites, "left") - 0.5 * get_aLndag_aLm(2*n, 2*m, aT, sites, "right") )
            end
        end
    end

    return res

end

function get_odd_even_taylor_groups(opsum, sites)

    """
    This function takes as input an OpSum and the sites of the MPS and returns nn_odd (nearest neighbour odd), nn_even and taylor.
    
    The gates considered odd are 4 site gates having a left most index 1:4:l-3 where l is the length of sites. 
    
    The gates considered even are 4 site gates having a left most index 3:4:l-3. 
    
    Everything else is considered to be in the taylor group which will be taylor expanded in the
    trotterization scheme: odd/2 taylor/2 even taylor/2 odd/2 (second order trotterization).

    The nn_odd is a list of lists where each list holds all the gates lying within the group of sites "left most index of group up to left most index of group + 3"
    where the left most index of each group is defined by 1:4:l-3. Similarly for nn_even.
    """

    function pad_op(element, group_left)

        """
        This function takes an element from an OpSum e.g. a Scaled{ComplexF64, Prod{Op}} and a group_left index which
        stands for the beginning of a group of 4 sites.

        The purpose of this function is most easily explained with an example:

        element = 0.5,"Z",2,"S-",4
        group_left = 2

        Here group_left is 2 so our sites span 2, 3, 4, 5 and the operators in the element are going to be shifted
        to have indices 1, 2, 3, 4 and then accordingly we will place the operators on the sites putting identities where we had no
        operators resulting in:

        0.5,"I",1,"I",2,"I",3,"I",4,"Z",1,"S-",3 where the first part is padding with identities and then the operators are appropriately placed within 1, 2, 3, 4

        The reason we need this function is to do build the exponential of gates for time evolution where the ITensor function will get the opsum of many of this
        elements but the indices will not be 1, 2, 3, 4 but rather sites[group_left:group_left+3] where sites are the actual site indices of the MPS
        """

        element_coefficient = element.args[1]      
        element_list = [[ITensors.which_op(element[i]), ITensors.site(element[i]) - group_left + 1] for i in 1:length(element)] # result e.g. [["Z", 1], ["S-", 3]] using the example above
        el = [["I", i] for i in 1:4] # put identities on all sites to cover the ones which don't have a non-trivial operator acting on them in the element (padding)
        for e in element_list # put in the non-trivial operators as gathered from the element input
            push!(el, [e[1], e[2]]) 
        end
        # Put both the identities and the non-trivial operators from the element input into one opsum
        complete_opsum = OpSum()
        complete_opsum += element_coefficient,el[1][1],el[1][2]
        for i in 2:length(el)
            complete_opsum *= el[i][1],el[i][2]
        end

        return complete_opsum

    end

    l = length(sites)
    taylor = OpSum() # defined as anything that does not fall into the following two categories below
    nn_odd = [OpSum() for _ in 1:4:l-3] # defined as anything that lies between 1 and 4 inclusive, 5 and 8, 9 and 12 etc
    nn_even = [OpSum() for _ in 3:4:l-3] # defined as anything that lies between 3 and 6 inclusive, 7 and 10
    
    odd_index_groups_list = [[i, i+1, i+2, i+3] for i in 1:4:l-3]
    odd_index_groups = reshape(reduce(vcat, odd_index_groups_list'), (length(odd_index_groups_list), 4)) # make the list of lists into a matrix for convenience
    even_index_groups_list = [odd_index_groups_list[i] .+ 2 for i in 1:length(odd_index_groups_list)-1]
    even_index_groups = reshape(reduce(vcat, even_index_groups_list'), (length(even_index_groups_list), 4)) # make the list of lists into a matrix for convenience
    
    for element in opsum 

        idxs = ITensors.sites(element) # where the element acts non-trivially
        left, right = minimum(idxs), maximum(idxs) # left and right most sites the gate is acting on
        
        gate_span = right - left + 1
        if gate_span <= 4 # if it can fit within an even or odd gate

            odd_group_of_left_index = findfirst(x -> x == left, odd_index_groups) # the index of the group to which the element can belong to (this can also return nothing when the element cannot belong to an odd group)
            if odd_group_of_left_index !== nothing
                odd_exists = true
                odd_group_of_left_index = odd_group_of_left_index[1] # just take the index of the group and discard the second element which is the position within the group
            else
                odd_exists = false
            end
            if odd_exists
                # the left and right most sites of the odd group the element can belong to
                min_of_odd_group, max_of_odd_group = minimum(odd_index_groups[odd_group_of_left_index, :]), maximum(odd_index_groups[odd_group_of_left_index, :])
                if left >= min_of_odd_group && right <= max_of_odd_group
                    odd_possible = true
                else
                    odd_possible = false
                end
            else
                odd_possible = false
            end
    
            # same as above but for even
            even_group_of_left_index = findfirst(x -> x == left, even_index_groups)
            if even_group_of_left_index !== nothing
                even_exists = true
                even_group_of_left_index = even_group_of_left_index[1]
            else
                even_exists = false
            end
            if even_exists
                min_of_even_group, max_of_even_group = minimum(even_index_groups[even_group_of_left_index, :]), maximum(even_index_groups[even_group_of_left_index, :])
                if left >= min_of_even_group && right <= max_of_even_group
                    even_possible = true
                else
                    even_possible = false
                end
            else
                even_possible = false
            end
    
            if odd_possible && even_possible
    
                if length(nn_odd[odd_group_of_left_index]) < length(nn_even[even_group_of_left_index]) # if the gate can go in both an odd and an even choose the one with less number of gates already
    
                    nn_odd[odd_group_of_left_index] += pad_op(element, min_of_odd_group)
                    
                else
    
                    nn_even[even_group_of_left_index] += pad_op(element, min_of_even_group)

                end
            
            elseif odd_possible
    
                nn_odd[odd_group_of_left_index] += pad_op(element, min_of_odd_group)

    
            elseif even_possible
    
                nn_even[even_group_of_left_index] += pad_op(element, min_of_even_group)
    
            else
    
                taylor += element
    
            end
            
        else
    
            taylor += element
    
        end
    end

    return nn_odd, nn_even, taylor

end

function apply_swap_to_pairs(mps, swap_tensor)

    sites = siteinds(mps)

    res = MPS(sites)

    for i in 1:2:length(mps)-1

        gate = ITensor(swap_tensor, dag(sites[i]), dag(sites[i+1]), dag(sites[i]'), dag(sites[i+1]'))

        tmp = noprime(mps[i]*mps[i+1]*gate)
        
        res[i], res[i+1] = ITensors.qr(tmp, commoninds(mps[i], tmp); tags = "Link,l=$(i)")
    
    end

    return res

end

function hermitian_conjugate_purified_density_matrix_mps(mps, swap_tensor)

    return dag(apply_swap_to_pairs(mps, swap_tensor))

end

function purified_density_matrix_from_mps_to_matrix(mps)

    n = div(length(mps), 2)
    a = contract(mps)
    top_inds = []
    bottom_inds = []
    i = 0
    for ind in inds(a)
        if i % 2 == 0
            push!(top_inds, ind)
        else
            push!(bottom_inds, ind)
        end
        i += 1
    end
    a = Array(a, top_inds..., bottom_inds...)
    a = reshape(a, 2^n, 2^n)
    
    return a

end

function get_double_aH_Hamiltonian_without_l0_terms(sites, x, ma, lambda, side)

    """
    This gives aH Hamiltonian

    side specifies "left" or "right" to imply H tensor product I or vice versa

    """

    N = div(length(sites), 2)

    opsum = OpSum()

    for n in 1:N-1

        if side == "left"
            n_idx = 2*n-1
        else
            n_idx = 2*n 
        end
        
        for m in n+1:N

            if side == "left"
                m_idx = 2*m-1
            else
                m_idx = 2*m
            end
            
            # Long range ZZ interaction term
            opsum += 0.25*(1/x)*(N-m+lambda),"Z",n_idx,"Z",m_idx

        end

        # Kinetic term
        opsum += 0.5,"S+",n_idx,"S-",n_idx+2
        opsum += 0.5,"S-",n_idx,"S+",n_idx+2
        
        opsum += (0.5*ma*(-1)^(n-1)),"Z",n_idx

    end

    if side == "left"
        opsum += (0.5*ma*(-1)^(N-1)),"Z",2*N-1
    else
        opsum += (0.5*ma*(-1)^(N-1)),"Z",2*N
    end

    return opsum

end

function get_Lindblad_opsum_without_l0_terms(sites, x, ma, lambda, aT, aD, env_corr_type, inputs, dissipator_sites)

    N = div(length(sites), 2)
    res = -1im*get_double_aH_Hamiltonian_without_l0_terms(sites, x, ma, lambda, "left")
    res += 1im*get_double_aH_Hamiltonian_without_l0_terms(sites, x, ma, lambda, "right")
    
    if aD != 0
        for n in dissipator_sites
            for m in dissipator_sites
                res += environment_correlator(env_corr_type, n, m, aD, inputs) * ( get_aLm_aLndag(2*n, 2*m-1, aT, sites) - 0.5 * get_aLndag_aLm(2*n-1, 2*m-1, aT, sites, "left") - 0.5 * get_aLndag_aLm(2*n, 2*m, aT, sites, "right") )
            end
        end
    end

    return res

end

function get_double_aH_Hamiltonian_just_l0_terms(sites, x, l_0, lambda, side)

    """
    This gives aH Hamiltonian

    side specifies "left" or "right" to imply H tensor product I or vice versa

    """

    N = div(length(sites), 2)

    opsum = OpSum()

    for n in 1:N-1

        if side == "left"
            n_idx = 2*n-1
        else
            n_idx = 2*n 
        end
        
        opsum += (1/x)*(N/8 - 0.25*ceil((n-1)/2) + l_0*(N-n)/2),"Z",n_idx
        
    end

    if side == "left"
        opsum += ((l_0^2)*(N-1)/(2*x) + (l_0*N)/(4*x) + (N^2)/(16*x) + (lambda*N/(8*x))),"Id",1
    else
        opsum += ((l_0^2)*(N-1)/(2*x) + (l_0*N)/(4*x) + (N^2)/(16*x) + (lambda*N/(8*x))),"Id",2
    end

    return opsum

end

function get_Lindblad_opsum_just_l0_terms(sites, x, l_0, lambda)

    res = -1im*get_double_aH_Hamiltonian_just_l0_terms(sites, x, l_0, lambda, "left")
    res += 1im*get_double_aH_Hamiltonian_just_l0_terms(sites, x, l_0, lambda, "right")
    
    return res

end

function get_applied_field(which_applied_field, inputs, t_over_a)

    l_0_1 = inputs["l_0_1"]
    
    if which_applied_field == "constant"
        return l_0_1
    else
        l_0_2 = inputs["l_0_2"]
        a_omega = inputs["a_omega"]
        if which_applied_field == "sauter"
            return l_0_1 + l_0_2/cosh(a_omega*t_over_a)^2
        elseif which_applied_field == "gaussian"
            return l_0_1 + l_0_2*exp(-(a_omega*t_over_a)^2)
        else # oscillatory case
            return l_0_1 + l_0_2*cos(a_omega*t_over_a)
        end
    end

end

function get_rho_dagger_rho_purified(mps)

    """
    This function takes a purified mps as input which represents rho the density matrix and performs 
    the operation rho -> rho^dagger * rho to return another purified mps representing the result
    """

    mps_left = hermitian_conjugate_purified_density_matrix_mps(mps)

    mps_right_contracted = MPS(div(length(mps), 2))
    mps_left_contracted = MPS(div(length(mps_left), 2))
    res = MPS(div(length(mps_left), 2))

    for i in 1:length(mps_right_contracted)

        mps_right_contracted[i] = mps[2*i-1]*mps[2*i] 
        mps_left_contracted[i] = prime(mps_left[2*i-1]*mps_left[2*i])

    end

    for i in 1:length(mps_right_contracted)

        idxs = inds(mps_left_contracted[i]; :tags => "Site,n=$(2*i)"), inds(mps_right_contracted[i]; :tags => "Site,n=$(2*i-1)")
        D = dag(delta(idxs))
        res[i] = mps_right_contracted[i] * mps_left_contracted[i] * D

    end

    for i in 2:length(res)-1

        C_left = combiner(inds(res[i]; :tags => "Link,l=$(2*i-2)"), tags = "Link,l=$(2*i-2)", dir = ITensors.In)
        C_right = combiner(inds(res[i]; :tags => "Link,l=$(2*i)"), tags = "Link,l=$(2*i)", dir = ITensors.Out)

        res[i-1] *= dag(C_left)
        res[i] = res[i] * C_left * C_right
        res[i+1] *= dag(C_right)

    end

    for i in 1:length(res)
        if i == 1
            swaptags!(res[i], "Link,l=$(2*i)", "Link,l=$(i)")
        elseif i == length(res)
            swaptags!(res[i], "Link,l=$(2*i-2)", "Link,l=$(i-1)")
        else
            swaptags!(res[i], "Link,l=$(2*i-2)", "Link,l=$(i-1)")
            swaptags!(res[i], "Link,l=$(2*i)", "Link,l=$(i)")
        end
    end

    for i in 1:length(res)

        swaptags!(res[i], "S=1/2,Site,n=$(2*i)", "S=1/2,Site,n=$(i)"; :plev => 0)
        swaptags!(res[i], "S=1/2,Site,n=$(2*i-1)", "S=1/2,Site,n=$(i)"; :plev => 1)

    end

    return rho_vec_to_mps(res)
    
end

function mpo_from_purified_mps(mps)

    mps = deepcopy(mps)

    N = length(mps)
    for i in 1:2:N
        prime!(mps[i]; :tags => "Site")
    end
    half_mps = MPS(div(N, 2))
    for i in 1:2:N
        half_mps[div(i, 2) + 1] = mps[i]*mps[i+1]
    end
    for i in 1:length(half_mps)
        if i == 1
            swaptags!(half_mps[i], "Link,l=$(2*i)", "Link,l=$(i)")
        elseif i == length(half_mps)
            swaptags!(half_mps[i], "Link,l=$(2*i-2)", "Link,l=$(i-1)")
        else
            swaptags!(half_mps[i], "Link,l=$(2*i-2)", "Link,l=$(i-1)")
            swaptags!(half_mps[i], "Link,l=$(2*i)", "Link,l=$(i)")
        end
    end
    for i in 1:length(half_mps)
        swaptags!(half_mps[i], "S=1/2,Site,n=$(2*i)", "S=1/2,Site,n=$(i)"; :plev => 0)
        swaptags!(half_mps[i], "S=1/2,Site,n=$(2*i-1)", "S=1/2,Site,n=$(i)"; :plev => 1)
    end
    mpo = convert(MPO, half_mps)
    return mpo

end
