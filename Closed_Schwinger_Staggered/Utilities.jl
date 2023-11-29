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

function get_exp_Ho(sites, a, x)::Vector{ITensor}

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

function get_exp_He(sites, a, x, l_0)::Vector{ITensor}

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
        res = ITensors.Array(res, inds_res)
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
        res = ITensors.Array(res, inds_res)
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

function apply_Ho_mpo_list!(Ho_mpo_list, mps; cutoff = 1e-9)

    for (idx_num, idx) in enumerate(1:2:N-1)

        gate = Ho_mpo_list[idx_num]
        
        tmp = noprime!(gate*mps[idx]*mps[idx+1])

        # println("In odd, for idx = $idx, ", get_which_canonical_form(mps))
        
        U, S, V = ITensors.svd(tmp, uniqueinds(mps[idx], mps[idx+1]), lefttags = "Link,l=$(idx)", righttags = "Link,l=$(idx)", cutoff = cutoff)
        
        mps[idx] = U

        S = S/norm(S)
        
        mps[idx+1] = S*V

        # Extra SVD for ATTDMRG compared to TEBD
        if idx != N-1
        
            idx += 1
        
            tmp = mps[idx]*mps[idx+1]
        
            U,S,V = ITensors.svd(tmp, uniqueinds(mps[idx], mps[idx+1]), lefttags = "Link,l=$(idx)", righttags = "Link,l=$(idx)", cutoff = cutoff)
        
            mps[idx] = U

            S = S/norm(S)
        
            mps[idx+1] = S*V
        
        end
    end

end

function apply_He_mpo_list!(He_mpo_list, mps; cutoff = 1e-9)

    """
    After we apply 1-tau*Hz/2 with the apply function we end up with right canonical form.
    To apply the first even gate which starts at the second site we need to have the first site in Left canonical form
    """

    # mps[1], mps[2] = ITensors.qr(mps[1]*mps[2], uniqueinds(mps[1], mps[2]); positive = true, tags = "Link,l=$(1)")
    mps[1], S, V = ITensors.svd(mps[1]*mps[2], uniqueinds(mps[1], mps[2]), lefttags = "Link,l=$(1)", righttags = "Link,l=$(1)")
    mps[2] = S*V

    for (idx_num, idx) in enumerate(2:2:N-2)

        gate = He_mpo_list[idx_num]
        
        tmp = noprime!(gate*mps[idx]*mps[idx+1])

        # println("In even, for idx = $idx, ", get_which_canonical_form(mps))
        
        U, S, V = ITensors.svd(tmp, uniqueinds(mps[idx], mps[idx+1]), lefttags = "Link,l=$(idx)", righttags = "Link,l=$(idx)", cutoff = cutoff)
        
        mps[idx] = U

        S = S/norm(S)
        
        mps[idx+1] = S*V
        
        # Extra SVD for ATTDMRG compared to TEBD

        idx += 1
    
        tmp = noprime!(mps[idx]*mps[idx+1])
    
        U,S,V = ITensors.svd(tmp, uniqueinds(mps[idx], mps[idx+1]), lefttags = "Link,l=$(idx)", righttags = "Link,l=$(idx)", cutoff = cutoff)
    
        mps[idx] = U

        S = S/norm(S)
        
        mps[idx+1] = S*V
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
