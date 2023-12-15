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

function get_Hamiltonian(sites, x, l_0, mg)

    N = length(sites)

    opsum = OpSum()

    for n in 1:N-1
        
        for k in n+1:N
            
            opsum += 4*0.5*(N-k),"Sz",n,"Sz",k

        end

        opsum += x,"S+",n,"S-",n+1
        opsum += x,"S-",n,"S+",n+1

        opsum += 2*(N/4 - 0.5*ceil((n-1)/2) + l_0*(N-n)),"Sz",n
        
        opsum += 2*(mg*sqrt(x)*(-1)^(n-1)),"Sz",n

    end

    opsum += 2*(mg*sqrt(x)*(-1)^(N-1)),"Sz",N

    opsum += ((l_0^2)*(N-1) + (l_0*N/2) + (N^2)/8),"Id",1

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
