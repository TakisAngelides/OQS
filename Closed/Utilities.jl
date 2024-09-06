function get_exp_Hz(sites, a, x)

    """
    a = prefactor of H_z eg: -i tau / 2 to give exp(-i * tau * Hz / 2)

    This operator includes the ZZ long range interactions, the mass term and the single Z term from the electric field term

    """
    N = length(sites)

    opsum = OpSum()

    for n in 1:N-1
        
        for k in n+1:N

            opsum += (0.25/x)*(N-k),"Z",n,"Z",k

        end

    end

    mpo = a*MPO(opsum, sites)
 
    return mpo

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

    return add(tmp1, tmp2; cutoff = 0)

end

function get_Hz(sites, x, l_0, mg)

    # READ: This will not be the same as the get_exp_Hz function it was changed to match OQS Hamiltonian

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

function get_exp_Ho(sites, a, x, ma, l_0)::Vector{ITensor}

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
        hj += ((N/8 - 0.25*ceil((n-1)/2) + l_0*(N-n)/2)/x) * op("Z", sites[n]) * op("I", sites[n+1])
        hj += (0.5*ma*(-1)^(n-1)) * op("Z", sites[n]) * op("I", sites[n+1])
        if n == N-1
            hj += (0.5*ma*(-1)^(N-1)) * op("I", sites[n]) * op("Z", sites[n+1]) # last mass term on site N
            hj += ((l_0^2)*(N-1)/(2*x) + (l_0*N/(4*x)) + (N^2)/16) * op("I", sites[n]) * op("I", sites[n+1]) # once the constant term
        end
        Gj = exp(a * hj)
        push!(gates, Gj)

    end

    return gates

end

function get_exp_He(sites, a, x, ma, l_0)::Vector{ITensor}

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
        hj += ((N/8 - 0.25*ceil((n-1)/2) + l_0*(N-n)/2)/x) * op("Z", sites[n]) * op("I", sites[n+1])
        hj += (0.5*ma*(-1)^(n-1)) * op("Z", sites[n]) * op("I", sites[n+1])
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

function mpo_to_matrix(mpo)

    n = length(mpo)
    a = contract(mpo)
    a = Array(a, inds(a; :plev => 1)..., inds(a; :plev => 0)...)
    a = reshape(a, 2^n, 2^n)

    return a

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

function apply_Ho_mpo_list!(Ho_mpo_list, mps; cutoff = 0, maxdim = 1000)

    N = length(mps)

    for (idx_num, idx) in enumerate(1:2:N-1)

        gate = Ho_mpo_list[idx_num]
        
        tmp = noprime!(gate*mps[idx]*mps[idx+1])

        # println("In odd, for idx = $idx, ", get_which_canonical_form(mps))
        
        U, S, V = ITensors.svd(tmp, uniqueinds(mps[idx], mps[idx+1]), lefttags = "Link,l=$(idx)", righttags = "Link,l=$(idx)", cutoff = cutoff, maxdim = maxdim)
        
        mps[idx] = U

        S = S/norm(S)
        
        mps[idx+1] = S*V

        # Extra SVD for ATTDMRG compared to TEBD
        if idx != N-1
        
            idx += 1
        
            tmp = mps[idx]*mps[idx+1]
        
            U,S,V = ITensors.svd(tmp, uniqueinds(mps[idx], mps[idx+1]), lefttags = "Link,l=$(idx)", righttags = "Link,l=$(idx)", cutoff = cutoff, maxdim = maxdim)
        
            mps[idx] = U

            S = S/norm(S)
        
            mps[idx+1] = S*V
        
        end
    end

end

function apply_He_mpo_list!(He_mpo_list, mps; cutoff = 0, maxdim = 1000)

    """
    After we apply 1-tau*Hz/2 with the apply function we end up with right canonical form.
    To apply the first even gate which starts at the second site we need to have the first site in Left canonical form
    """

    N = length(mps)

    # mps[1], mps[2] = ITensors.qr(mps[1]*mps[2], uniqueinds(mps[1], mps[2]); positive = true, tags = "Link,l=$(1)")
    mps[1], S, V = ITensors.svd(mps[1]*mps[2], uniqueinds(mps[1], mps[2]), lefttags = "Link,l=$(1)", righttags = "Link,l=$(1)"; cutoff = cutoff, maxdim = maxdim)
    mps[2] = S*V

    for (idx_num, idx) in enumerate(2:2:N-2)

        gate = He_mpo_list[idx_num]
        
        tmp = noprime!(gate*mps[idx]*mps[idx+1])

        # println("In even, for idx = $idx, ", get_which_canonical_form(mps))
        
        U, S, V = ITensors.svd(tmp, uniqueinds(mps[idx], mps[idx+1]), lefttags = "Link,l=$(idx)", righttags = "Link,l=$(idx)", cutoff = cutoff, maxdim = maxdim)
        
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
    
        push!(res, inner(psi', MPO(get_Z_site_operator(i), sites), psi))

    end

    return res

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

    op_dict = Dict("I" => sparse([1 0; 0 1]), "X" => sparse([0 1; 1 0]), "Y" => sparse([0 -1im; 1im 0]), "Z" => sparse([1 0; 0 -1]), "S-" => sparse([0 0; 1 0]), "S+" => sparse([0 1; 0 0]))
    
    zipped = [(i, pos, op) for ((i, op), pos) in zip(enumerate(ops), positions)]
    zipped = sort(zipped, by = x -> (x[2], x[1]))
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

function get_pseudomomentum_operator(N, x)

    op = spzeros(2^(N), 2^(N))

    for n in 1:N-2

        op += get_op(["S-", "Z", "S+"], [n, n+1, n+2], N)
        op += get_op(["S+", "Z", "S-"], [n, n+1, n+2], N)

    end

    op *= -1im*x

    return op

end
