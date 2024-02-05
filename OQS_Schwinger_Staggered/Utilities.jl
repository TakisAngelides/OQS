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

function apply_odd!(Ho_mpo_list, mpo; cutoff = 1e-9)

    N = length(mpo)

    for (idx, gate) in enumerate(Ho_mpo_list)

        idx = 2*idx-1

        tmp = replaceprime(prime(gate'*mpo[idx]*mpo[idx+1]; :tags => "Site")*conj(gate), 3 => 1)

        U, S, V = ITensors.svd(tmp, commoninds(tmp, mpo[idx])..., lefttags = "Link,l=$(idx)", righttags = "Link,l=$(idx)", cutoff = cutoff)
        
        mpo[idx] = U

        S = S/norm(S)
        
        mpo[idx+1] = S*V

        # Extra SVD for ATTDMRG compared to TEBD
        if idx != N-1
        
            idx += 1
        
            tmp = mpo[idx]*mpo[idx+1]
        
            U,S,V = ITensors.svd(tmp, commoninds(tmp, mpo[idx])..., lefttags = "Link,l=$(idx)", righttags = "Link,l=$(idx)", cutoff = cutoff)
        
            mpo[idx] = U

            S = S/norm(S)
        
            mpo[idx+1] = S*V
        
        end

    end

end

function apply_even!(He_mpo_list, mpo; cutoff = 1e-9)

    """
    After we apply 1-tau*Hz/2 with the apply function we end up with right canonical form.
    To apply the first even gate which starts at the second site we need to have the first site in Left canonical form
    """

    # mps[1], mps[2] = ITensors.qr(mps[1]*mps[2], uniqueinds(mps[1], mps[2]); positive = true, tags = "Link,l=$(1)")
    mpo[1], S, V = ITensors.svd(mpo[1]*mpo[2], commoninds(mpo[1], mpo[1]*mpo[2])..., lefttags = "Link,l=$(1)", righttags = "Link,l=$(1)")
    S = S/norm(S)
    mpo[2] = S*V

    for (idx, gate) in enumerate(He_mpo_list)

        idx = 2*idx

        tmp = replaceprime(prime(gate'*mpo[idx]*mpo[idx+1]; :tags => "Site")*conj(gate), 3 => 1)

        U, S, V = ITensors.svd(tmp, commoninds(tmp, mpo[idx])..., lefttags = "Link,l=$(idx)", righttags = "Link,l=$(idx)", cutoff = cutoff)
        
        mpo[idx] = U

        S = S/norm(S)
        
        mpo[idx+1] = S*V

        # Extra SVD for ATTDMRG compared to TEBD
        
        idx += 1
    
        tmp = mpo[idx]*mpo[idx+1]
    
        U,S,V = ITensors.svd(tmp, commoninds(tmp, mpo[idx])..., lefttags = "Link,l=$(idx)", righttags = "Link,l=$(idx)", cutoff = cutoff)
    
        mpo[idx] = U

        S = S/norm(S)
    
        mpo[idx+1] = S*V
    
    end

end

function get_entanglement_entropy(psi, site)
    
    ITensors.orthogonalize!(psi, site)
    if site == 1
        U,S,V = svd(psi[site], siteind(psi, site))
    else
        U,S,V = svd(psi[site], (linkind(psi, site-1), siteind(psi, site)))
    end
    SvN = 0.0
    for n=1:dim(S, 1)
        p = S[n,n]^2
        SvN -= p * log(p)
    end
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
        
        opsum += (-1)^(n-1),"Z",n

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
    X = sparse(Float64[0 1; 1 0])
    Y = sparse(ComplexF64[0 -1im; 1im 0])
    Z = sparse(Float64[1 0; 0 -1])

    # Kinetic term
    for n=1:N-1
        Hk += (1/4)*kron(eye(2^(n-1)), kron(X, kron(X, eye(2^(N-(n+1))))))
        Hk += (1/4)*kron(eye(2^(n-1)), kron(Y, kron(Y, eye(2^(N-(n+1))))))
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

function get_Lindblad_jump_operator(m, aT, sites)

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

function get_LdagL(n, m, aT, sites)

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
        opsum += -1im*0.5,"S+",n,"S-",n+1
        opsum += -1im*0.5,"S-",n,"S+",n+1
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
        opsum += 1im*0.5,"S+",n+N,"S-",n+1+N
        opsum += 1im*0.5,"S-",n+N,"S+",n+1+N
        opsum += 1im*(1/x)*(N/8 - 0.25*ceil((n-1)/2) + l_0*(N-n)/2),"Z",n+N
        opsum += 1im*(0.5*ma*(-1)^(n-1)),"Z",n+N
    end
    opsum += 1im*(0.5*ma*(-1)^(N-1)),"Z",N+N
    opsum += 1im*((l_0^2)*(N-1)/(2*x) + (l_0*N)/(4*x) + (N^2)/(16*x)),"Id",1+N

    h1 = MPO(opsum, sites)

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
            h2 = MPO(aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a)*opsum, sites)
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
            h3 = MPO(aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a)*opsum, sites)
            h1 += apply(h2, h3, cutoff = 1e-20)
            
            # # # The Ldag_n_L_m_left
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
            
            # # The Ldag_n_L_m_right
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
            h1 += MPO(opsum, sites)

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
    L += -1im * kron(H, eye(2^N)) + 1im * kron(eye(2^N), transpose(H)) 

    for n in 1:N
        for m in 1:N

            tmp1 = get_Lindblad_jump_operator_sparse_matrix(N, n, aT)
            tmp2 = get_Lindblad_jump_operator_sparse_matrix(N, m, aT)

            tmp3 = tmp1' * tmp2 # the dash is the dagger
            
            L += aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * ((kron(tmp2, conj(tmp1))) - 0.5*(kron(tmp3, eye(2^N))) -0.5*(kron(eye(2^N), transpose(tmp3))))

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

    # tmp1 = MPO(sites, "Id") + 0.5*1im*tau*H_T
    # tmp2 = MPO(sites, "Id") - 0.5*1im*tau*H_T
    # rho = apply(apply(tmp2, rho; cutoff = cutoff), tmp1; cutoff = cutoff)

    tmp = 0.5*1im*tau*H_T

    # rho + rho * idt/2 H_T * rho - idt/2 H_T * rho = ()
    rho_final = rho + apply(rho, tmp; cutoff = cutoff) - apply(tmp, rho; cutoff = cutoff)

    # second order term in the taylor expansion only for the Hamiltonian part
    # rho_final += (-tau^2/8)*apply(H_T, apply(H_T, rho; cutoff = cutoff); cutoff = cutoff) + (-tau^2/8)*apply(rho, apply(H_T, H_T; cutoff = cutoff); cutoff = cutoff) + (tau^2/4)*apply(H_T, apply(rho, H_T; cutoff = cutoff); cutoff = cutoff)

    # - 0.5 * dt * 0.5 * L_n^\dagger L_m
    LdagnLm = -0.5 * tau * 0.5 * get_Lindblad_dissipative_part(aD_0, sigma_over_a, env_corr_type, aT, sites)

    # sum over n and m: - 0.5 * dt/2 * L_n^\dagger L_m * rho
    rho_final += apply(LdagnLm, rho; cutoff = cutoff) 
    
    # sum over n and m: - 0.5 * dt/2 * rho * L_n^\dagger L_m 
    rho_final += apply(rho, LdagnLm; cutoff = cutoff) 

    # sum over n and m: aD_0 * f(a(n-m)) * L_m * rho * L_n^\dagger
    for n=1:N
        for m=1:N

            # aD_0 * f(a(n-m)) * L_m
            Lm = 0.5 * tau * aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * get_Lindblad_jump_operator(m, aT, sites)

            # L_n^\dagger (we only conjugate and not transpose cause we will apply it on the right of rho anyway)
            Ldagn = conj(get_Lindblad_jump_operator(n, aT, sites))

            # 0.5 * dt * aD_0 * f(a(n-m)) * L_m * rho * L_n^\dagger
            rho_final += apply(apply(Lm, rho; cutoff = cutoff), Ldagn; cutoff = cutoff)
    
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
