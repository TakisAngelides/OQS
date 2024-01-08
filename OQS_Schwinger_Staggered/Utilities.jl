function get_exp_Hz(sites, a, x, l_0, mg)

    """
    a = prefactor of H_z eg: -i tau / 2 to give exp(-i * tau * Hz / 2)

    This operator includes the ZZ long range interactions, the mass term and the single Z term from the electric field term

    """

    N = length(sites)

    opsum = OpSum()

    for n in 1:N-1
        
        for k in n+1:N

            # Z_n Z_k long range interaction
            opsum += a*4*0.5*(N-k),"Sz",n,"Sz",k

        end

        # Single Z terms coming from the electric field term
        opsum += a*2*(N/4 - 0.5*ceil((n-1)/2) + l_0*(N-n)),"Sz",n

        # Single Z terms from mass term
        opsum += a*2*(mg*sqrt(x)*(-1)^(n-1)),"Sz",n

    end

    # For loop above goes from n = 1 to N-1 but mass term also has n = N
    opsum += a*2*(mg*sqrt(x)*(-1)^(N-1)),"Sz",N

    # Constant term
    opsum += a*((l_0^2)*(N-1) + (l_0*N/2) + (N^2)/8),"Id",1

    mpo = MPO(opsum, sites)

    # exp(aHz/2) -> I + aHz/2
    final_mpo = MPO(sites, "Id") + mpo 

    return final_mpo

end

function get_Hz(sites, x, l_0, mg)

    """
    a = prefactor of H_z eg: -i tau / 2 to give exp(-i * tau * Hz / 2)

    This operator includes the ZZ long range interactions, the mass term and the single Z term from the electric field term

    """

    N = length(sites)

    opsum = OpSum()

    for n in 1:N-1
        
        for k in n+1:N

            opsum += 4*0.5*(N-k),"Sz",n,"Sz",k

        end

        opsum += 2*(N/4 - 0.5*ceil((n-1)/2) + l_0*(N-n)),"Sz",n

        opsum += 2*(mg*sqrt(x)*(-1)^(n-1)),"Sz",n

    end

    opsum += 2*(mg*sqrt(x)*(-1)^(N-1)),"Sz",N

    opsum += ((l_0^2)*(N-1) + (l_0*N/2) + (N^2)/8),"Id",1

    mpo = MPO(opsum, sites)

    return mpo

end

function get_exp_Ho_list(sites, a, x)::Vector{ITensor}

    """

    a = prefactor of H_o eg: -i tau / 2 to give exp(-i * tau * Ho / 2)

    This list of operators incldues the odd terms of the kinetic term

    Note: (XX+YY)/2 = S+S- + S-S+

    """

    gates = []
    N = length(sites)

    for n=1:2:(N-1)

        hj = x * op("S-", sites[n]) * op("S+", sites[n+1])
        hj += x * op("S+", sites[n]) * op("S-", sites[n+1])
        Gj = exp(a * hj)
        push!(gates, Gj)

    end

    return gates

end

function get_exp_He_list(sites, a, x)::Vector{ITensor}

    """

    a = prefactor of H_e eg: -i tau to give exp(-i * tau * He)

    This list of operators incldues the even terms of the kinetic term and we also includes the identity operator of H

    Note: (XX+YY)/2 = S+S- + S-S+

    """

    gates = []
    N = length(sites)

    for n=2:2:(N-2)

        hj = x * op("S-", sites[n]) * op("S+", sites[n+1])
        hj += x * op("S+", sites[n]) * op("S-", sites[n+1])
        Gj = exp(a * hj)
        push!(gates, Gj)

    end

    return gates

end

function get_W_Hamiltonian(sites, x, l_0, mg)

    """
    This gives W = 2*H/(ag^2)
    """

    N = length(sites)

    opsum = OpSum()

    for n in 1:N-1
        
        for k in n+1:N
            
            opsum += 0.5*(N-k),"Z",n,"Z",k

        end

        opsum += x,"S+",n,"S-",n+1
        opsum += x,"S-",n,"S+",n+1

        opsum += (N/4 - 0.5*ceil((n-1)/2) + l_0*(N-n)),"Z",n
        
        opsum += (mg*sqrt(x)*(-1)^(n-1)),"Z",n

    end

    opsum += (mg*sqrt(x)*(-1)^(N-1)),"Z",N

    opsum += ((l_0^2)*(N-1) + (l_0*N/2) + (N^2)/8),"Id",1

    return MPO(opsum, sites)

end

function get_aH_Hamiltonian(sites, x, l_0, ma)

    """
    This gives aH Hamiltonian
    """

    N = length(sites)

    opsum = OpSum()

    for n in 1:N-1
        
        for m in n+1:N
            
            # Long range ZZ interaction term
            opsum += 0.25*(1/x)*(N-m),"Z",n,"Z",m

        end

        # Kinetic term
        opsum += 0.5,"S+",n,"S-",n+1
        opsum += 0.5,"S-",n,"S+",n+1

        opsum += (1/x)*(N/8 - 0.25*ceil((n-1)/2) + l_0*(N-n)/2),"Z",n
        
        opsum += (0.5*ma*(-1)^(n-1)),"Z",n

    end

    opsum += (0.5*ma*(-1)^(N-1)),"Z",N

    opsum += ((l_0^2)*(N-1)/(2*x) + (l_0*N)/(4*x) + (N^2)/(16*x)),"Id",1

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

        tmp = replaceprime(prime(gate'*mpo[idx]*mpo[idx+1]; :tags => "Site")*gate, 3 => 1)

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

        tmp = replaceprime(prime(gate'*mpo[idx]*mpo[idx+1]; :tags => "Site")*gate, 3 => 1)

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
    
    orthogonalize!(psi, site)
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

    ampo += 2,"Sz",site

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
        
        opsum += 0.5*(-1)^(n-1),"Sz",n

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
    X = sparse(Float64[0 1; 1 0])
    Y = sparse(ComplexF64[0 -1im; 1im 0])
    Z = sparse(Float64[1 0; 0 -1])

    # Kinetic term
    for n=1:N-1
        H += (1/4)*kron(eye(2^(n-1)), kron(X, kron(X, eye(2^(N-(n+1))))))
        H += (1/4)*kron(eye(2^(n-1)), kron(Y, kron(Y, eye(2^(N-(n+1))))))
    end

    # Long range ZZ interaction term
    for n = 1:N-1
        for m = n+1:N
            H += (1/x)*(1/4)*(N - m + lambda)*kron(eye(2^(n-1)), kron(Z, kron(eye(2^(m-n-1)), kron(Z, eye(2^(N-m))))))
        end
    end

    # Mass term
    for n=1:N
        H += (ma/2)*((-1)^(n-1))*kron(eye(2^(n-1)), kron(Z, eye(2^(N-n))))
    end

    # Electric single Z term
    for n=1:N-1
        H += ((1/x)*(N/8 - (1/4)*ceil((n-1)/2) + l_0*(N-n)/2))*kron(eye(2^(n-1)), kron(Z, eye(2^(N-n))))
    end

    # Constant term
    H += (l_0^2*(N-1)/(2*x) + l_0*N/(4*x) + (N^2)/(16*x) + lambda*N/(8*x))*eye(2^N)

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

function get_Lindblad_jump_operator_sparse_matrix(N, m, aT)

    eye(n::Int64) = sparse(I, n, n);
    X = sparse(Float64[0 1; 1 0])
    Y = sparse(ComplexF64[0 -1im; 1im 0])
    Z = sparse(Float64[1 0; 0 -1])

    res = spzeros(2^(N), 2^(N))        
            
    res += 0.5*((-1)^m)*kron(eye(2^(m-1)), kron(Z, eye(2^(N-m))))
    res += 0.5*((-1)^m)*eye(2^N)
    
    if m != 1
        res += (1im*(-1)^m/(16*aT))*kron(eye(2^(m-2)), kron(X, kron(Y, eye(2^(N-m)))))
        res += (-1im*(-1)^m/(16*aT))*kron(eye(2^(m-2)), kron(Y, kron(X, eye(2^(N-m)))))
    end
    
    if m != N
        res += (-1im*(-1)^m/(16*aT))*kron(eye(2^(m-1)), kron(X, kron(Y, eye(2^(N-m-1)))))
        res += (1im*(-1)^m/(16*aT))*kron(eye(2^(m-1)), kron(Y, kron(X, eye(2^(N-m-1)))))
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
    X = sparse(Float64[0 1; 1 0])
    Y = sparse(ComplexF64[0 -1im; 1im 0])
    Z = sparse(Float64[1 0; 0 -1])

    H = get_aH_Hamiltonian_sparse_matrix(N, x, ma, l_0, lambda)

    # Unitary part of Lindbladian
    L += -1im * kron(H, eye(2^N)) + 1im * kron(eye(2^N), H) 

    for n in 1:N
        for m in 1:N

            tmp1 = get_Lindblad_jump_operator_sparse_matrix(N, n, aT)
            tmp2 = get_Lindblad_jump_operator_sparse_matrix(N, m, aT)

            tmp3 = conj.(tmp1) * tmp2

            L += aD_0 * environment_correlator(env_corr_type, n, m, aD_0, sigma_over_a) * (kron((tmp1), (tmp2')) - 0.5*kron((tmp3), eye(2^N)) -0.5*kron(eye(2^N), (tmp3')))

        end
    end

    return L

end
