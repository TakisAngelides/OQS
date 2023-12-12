function get_exp_Hz_list(sites, a, x, l_0, mg)

    """
    a = prefactor of H_z eg: -i tau / 2 to give exp(-i * tau * Hz / 2)

    This operator includes the ZZ long range interactions, the mass term and the single Z term from the electric field term

    Its a list because we have first the no prime sites and then the primed

    """

    n = length(sites)

    N = div(n, 2)

    opsum_noprime = OpSum()
    opsum_primed = OpSum()

    for n in 1:N-1

        n_noprime = 2*n-1
        n_primed = 2*n
        
        for k in n+1:N

            k_noprime = 2*k-1
            k_primed = 2*k

            # Z_n Z_k long range interaction
            opsum_noprime += a*4*0.5*(N-k),"Sz",n_noprime,"Sz",k_noprime
            opsum_primed += a*4*0.5*(N-k),"Sz",n_primed,"Sz",k_primed

        end

        # Single Z terms coming from the electric field term
        opsum_noprime += a*2*(N/4 - 0.5*ceil((n-1)/2) + l_0*(N-n)),"Sz",n_noprime
        opsum_primed += a*2*(N/4 - 0.5*ceil((n-1)/2) + l_0*(N-n)),"Sz",n_primed

        # Single Z terms from mass term
        opsum_noprime += a*2*(mg*sqrt(x)*(-1)^(n-1)),"Sz",n_noprime
        opsum_primed += a*2*(mg*sqrt(x)*(-1)^(n-1)),"Sz",n_primed

    end

    # For loop above goes from n = 1 to N-1 but mass term also has n = N
    opsum_noprime += a*2*(mg*sqrt(x)*(-1)^(N-1)),"Sz",(2*N-1)
    opsum_primed += a*2*(mg*sqrt(x)*(-1)^(N-1)),"Sz",(2*N)

    # Constant term
    opsum_noprime += a*((l_0^2)*(N-1) + (l_0*N/2) + (N^2)/8),"Id",1
    opsum_primed += a*((l_0^2)*(N-1) + (l_0*N/2) + (N^2)/8),"Id",2

    mpo_noprime = MPO(opsum_noprime, sites)
    mpo_primed = MPO(opsum_primed, sites)


    # exp(aHz/2) -> I + aHz/2
    final_mpo_noprime = MPO(sites, "Id") + mpo_noprime 
    final_mpo_primed = MPO(sites, "Id") + mpo_primed 

    return [final_mpo_noprime, final_mpo_primed]

end

function get_Hz_noprime(sites, x, l_0, mg)

    """
    a = prefactor of H_z eg: -i tau / 2 to give exp(-i * tau * Hz / 2)

    This operator includes the ZZ long range interactions, the mass term and the single Z term from the electric field term

    """

    N = length(sites)

    opsum = OpSum()

    for n in 1:N-1

        n_noprime = 2*n-1
        
        for k in n+1:N

            k_noprime = 2*k-1

            opsum += 4*0.5*(N-k),"Sz",n_noprime,"Sz",k_noprime

        end

        opsum += 2*(N/4 - 0.5*ceil((n-1)/2) + l_0*(N-n)),"Sz",n_noprime

        opsum += 2*(mg*sqrt(x)*(-1)^(n-1)),"Sz",n_noprime

    end

    opsum += 2*(mg*sqrt(x)*(-1)^(N-1)),"Sz",(2*N-1)

    opsum += ((l_0^2)*(N-1) + (l_0*N/2) + (N^2)/8),"Id",1

    mpo = MPO(opsum, sites)

    return mpo

end

function get_exp_Ho_list(sites, a, x)

    """

    a = prefactor of H_o eg: -i tau / 2 to give exp(-i * tau * Ho / 2)

    This list of operators incldues the odd terms of the kinetic term

    Note: (XX+YY)/2 = S+S- + S-S+

    Since the MPS is purified in the order 1, 1', 2, 2', ... where the ' is the column index of the density matrix
    the odd gates of the unprimed start at i and end at i+2 for i in 1:4:n-3 and the gates of the primed are just slided down one 
    so at i+1 to i+3

    """

    gates = []
    n = length(sites)

    for idx=1:4:(n-3)

        # no primes gates
        hj = x * op("S-", sites[idx]) * op("I", sites[idx+1]) * op("S+", sites[idx+2])
        hj += x * op("S+", sites[idx]) * op("I", sites[idx+1]) * op("S-", sites[idx+2])
        Gj = exp(a * hj)
        push!(gates, (Gj, idx))

        # primes gates
        hj = x * op("S-", sites[idx+1]) * op("I", sites[idx+2]) * op("S+", sites[idx+3])
        hj += x * op("S+", sites[idx+1]) * op("I", sites[idx+2]) * op("S-", sites[idx+3])
        Gj = exp(a * hj)
        push!(gates, (Gj, idx+1))

    end

    return gates

end

function get_exp_He_list(sites, a, x)

    """

    a = prefactor of H_e eg: -i tau to give exp(-i * tau * He)

    This list of operators incldues the even terms of the kinetic term and we also includes the identity operator of H

    Note: (XX+YY)/2 = S+S- + S-S+

    """

    gates = []
    n = length(sites)

    for idx=3:4:(n-5)

        # no primes gates
        hj = x * op("S-", sites[idx]) * op("I", sites[idx+1]) * op("S+", sites[idx+2])
        hj += x * op("S+", sites[idx]) * op("I", sites[idx+1]) * op("S-", sites[idx+2])
        Gj = exp(a * hj)
        push!(gates, (Gj, idx))

        # primes gates
        hj = x * op("S-", sites[idx+1]) * op("I", sites[idx+2]) * op("S+", sites[idx+3])
        hj += x * op("S+", sites[idx+1]) * op("I", sites[idx+2]) * op("S-", sites[idx+3])
        Gj = exp(a * hj)
        push!(gates, (Gj, idx+1))

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

function get_Hamiltonian_noprime(sites, x, l_0, mg)

    n = length(sites)

    N = div(n, 2)

    opsum = OpSum()

    for n in 1:N-1

        n_noprime = 2*n-1
        
        for k in n+1:N

            k_noprime = 2*k-1
            
            opsum += 4*0.5*(N-k),"Sz",n_noprime,"Sz",k_noprime

        end

        opsum += x,"S+",n_noprime,"S-",n_noprime+1
        opsum += x,"S-",n,"S+",n_noprime+1

        opsum += 2*(N/4 - 0.5*ceil((n-1)/2) + l_0*(N-n)),"Sz",n_noprime
        
        opsum += 2*(mg*sqrt(x)*(-1)^(n-1)),"Sz",n_noprime

    end

    opsum += 2*(mg*sqrt(x)*(-1)^(N-1)),"Sz",(2*N-1)

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

    l = length(Ho_mpo_list)

    for (element_num, element) in enumerate(Ho_mpo_list)

        if (element_num % 2) != 0

            gate, left_idx = element

            tmp = noprime(gate*mps[left_idx]*mps[left_idx+1]*mps[left_idx+2])
            
            U1, S, V = ITensors.svd(tmp, commoninds(tmp, mps[left_idx]), cutoff = cutoff, lefttags = "Link,l=$(left_idx)", righttags = "Link,l=$(left_idx)") 

            R = S*V/norm(S)

            U2, R = ITensors.qr(R, [commoninds(R, mps[left_idx + 1]; :tags => "Site")..., commoninds(R, U1)...]) 

            settags!(U2, "Link,l=$(left_idx+1)"; :tags => "qr")

            settags!(R, "Link,l=$(left_idx+1)"; :tags => "qr")

            mps[left_idx], mps[left_idx+1], mps[left_idx+2] = U1, U2, R

        else 

            gate, left_idx = element

            tmp = noprime(gate*mps[left_idx]*mps[left_idx+1]*mps[left_idx+2])
            
            U1, S1, V1 = ITensors.svd(tmp, commoninds(tmp, mps[left_idx]), cutoff = cutoff, lefttags = "Link,l=$(left_idx)", righttags = "Link,l=$(left_idx)") 

            V1 = S1*V1/norm(S1)

            U2, S2, V2 = ITensors.svd(V1, [commoninds(V1, mps[left_idx + 1]; :tags => "Site")..., commoninds(V1, U1)...], lefttags = "Link,l=$(left_idx+1)", righttags = "Link,l=$(left_idx+1)") 

            V2 = S2*V2/norm(S2)

            if element_num == l

                mps[left_idx], mps[left_idx+1], mps[left_idx+2] = U1, U2, V2

            else

                U3, S3, V3 = ITensors.svd(V2, [commoninds(V2, mps[left_idx + 2]; :tags => "Site")..., commoninds(V2, U2)...], lefttags = "Link,l=$(left_idx+2)", righttags = "Link,l=$(left_idx+2)") 

                V3 = S3*V3/norm(S3)

                mps[left_idx], mps[left_idx+1], mps[left_idx+2] = U1, U2, U3
                mps[left_idx + 3] *= V3

            end

        end
    
    end

end

function apply_He_mpo_list!(He_mpo_list, mps; cutoff = 1e-9)

    """
    After we apply 1-tau*Hz/2 with the apply function we end up with right canonical form.
    To apply the first even gate which starts at the third site we need to have the first two site in Left canonical form
    """

    # mps[1], S, V = ITensors.svd(mps[1]*mps[2], uniqueinds(mps[1], mps[2]), lefttags = "Link,l=$(1)", righttags = "Link,l=$(1)")
    # mps[2] = S*V
    # mps[2], S, V = ITensors.svd(mps[2]*mps[3], uniqueinds(mps[2], mps[3]), lefttags = "Link,l=$(2)", righttags = "Link,l=$(2)")
    # mps[3] = S*V
    mps[1], mps[2] = ITensors.qr(mps[1]*mps[2], uniqueinds(mps[1], mps[2]); positive = true, tags = "Link,l=$(1)")
    mps[2], mps[3] = ITensors.qr(mps[2]*mps[3], uniqueinds(mps[2], mps[3]); positive = true, tags = "Link,l=$(2)")

    l = length(Ho_mpo_list)

    for (element_num, element) in enumerate(Ho_mpo_list)

        if (element_num % 2) != 0

            gate, left_idx = element

            tmp = noprime(gate*mps[left_idx]*mps[left_idx+1]*mps[left_idx+2])
            
            U1, S, V = ITensors.svd(tmp, commoninds(tmp, mps[left_idx]), cutoff = cutoff, lefttags = "Link,l=$(left_idx)", righttags = "Link,l=$(left_idx)") 

            R = S*V/norm(S)

            U2, R = ITensors.qr(R, [commoninds(R, mps[left_idx + 1]; :tags => "Site")..., commoninds(R, U1)...]) 

            settags!(U2, "Link,l=$(left_idx+1)"; :tags => "qr")

            settags!(R, "Link,l=$(left_idx+1)"; :tags => "qr")

            mps[left_idx], mps[left_idx+1], mps[left_idx+2] = U1, U2, R

        else 

            gate, left_idx = element

            tmp = noprime(gate*mps[left_idx]*mps[left_idx+1]*mps[left_idx+2])
            
            U1, S1, V1 = ITensors.svd(tmp, commoninds(tmp, mps[left_idx]), cutoff = cutoff, lefttags = "Link,l=$(left_idx)", righttags = "Link,l=$(left_idx)") 

            V1 = S1*V1/norm(S1)

            U2, S2, V2 = ITensors.svd(V1, [commoninds(V1, mps[left_idx + 1]; :tags => "Site")..., commoninds(V1, U1)...], lefttags = "Link,l=$(left_idx+1)", righttags = "Link,l=$(left_idx+1)") 

            V2 = S2*V2/norm(S2)

            if element_num == l

                mps[left_idx], mps[left_idx+1], mps[left_idx+2] = U1, U2, V2

            else

                U3, S3, V3 = ITensors.svd(V2, [commoninds(V2, mps[left_idx + 2]; :tags => "Site")..., commoninds(V2, U2)...], lefttags = "Link,l=$(left_idx+2)", righttags = "Link,l=$(left_idx+2)") 

                V3 = S3*V3/norm(S3)

                mps[left_idx], mps[left_idx+1], mps[left_idx+2] = U1, U2, U3
                mps[left_idx + 3] *= V3

            end

        end
    
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
