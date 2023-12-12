using ITensors
using LinearAlgebra
include("Utilities.jl")

# --

# sites = siteinds("S=1/2", 4)
# x = 1.0
# a = 1.0
# n = 1

# # hj = x * op("S-", sites[n]) * op("S+", sites[n+1]) * op("I", sites[n+2]) * op("I", sites[n+3])
# # hj += x * op("S+", sites[n]) * op("S-", sites[n+1]) * op("I", sites[n+2]) * op("I", sites[n+3])
# # hj += x * op("I", sites[n]) * op("I", sites[n+1]) * op("S-", sites[n+2]) * op("S+", sites[n+3])
# # hj += x * op("I", sites[n]) * op("I", sites[n+1]) * op("S+", sites[n+2]) * op("S-", sites[n+3])

# hj = x * op("S-", sites[n]) * op("S+", sites[n+1])
# hj += x * op("S+", sites[n]) * op("S-", sites[n+1])

# Gj = exp(a * hj)
# mpo = MPO(Gj, sites)
# println(linkinds(mpo))
# # A = ITensors.Array(Gj, inds(Gj; :plev => 1)..., inds(Gj; :plev => 0)...)
# # A = reshape(A, 4, 4)
# # @show A

# # u, s, v = ITensors.svd(Gj, inds(Gj; :tags => 1)..., inds(Gj; :tags => 2))

# # @show diag(s)

# --

N = 8
n = 2*N
sites = siteinds("S=1/2", n, conserve_qns = true)
for (site_idx, site) in enumerate(sites)
    if site_idx % 2 == 0
        sites[site_idx] = dag(site)
    end
end 
state = fill("1", n)
for state_idx in 3:4:n-1
    state[state_idx] = "0"
    state[state_idx+1] = "0"
end
mps = randomMPS(sites, state; linkdims = 2)
orthogonalize!(mps, 1)

a = 1.0
tau = 1.0
x = 1.0
l_0 = 1.0
mg = 1.0
cut_off = 1e-16

# gates = []
# idx = 1
# hj = x * op("S-", sites[idx]) * op("I", sites[idx+1]) * op("S+", sites[idx+2])
# hj += x * op("S+", sites[idx]) * op("I", sites[idx+1]) * op("S-", sites[idx+2])
# Gj = exp(a * hj)
# left_idx, right_idx = idx, idx+2
# push!(gates, (Gj, left_idx))

# gate, left_idx = gates[1]

# # mpo = MPO(gate, [sites[left_idx:right_idx]])
# # println(linkdims(mpo))

# println(inds(mps))
# println(norm(mps))
# println(get_which_canonical_form(mps), " ", linkdims(mps))
# tmp = noprime(gate*mps[left_idx]*mps[left_idx+1]*mps[left_idx+2])
# U1, S, V = ITensors.svd(tmp, commoninds(tmp, mps[left_idx]), cutoff = 1e-9, lefttags = "Link,l=$(idx)", righttags = "Link,l=$(idx)") 
# R = S*V/norm(S)
# println("Norm of S: ", norm(S), " ", diag(S))
# U2, R = ITensors.qr(R, [commoninds(R, mps[left_idx + 1]), commoninds(R, U1)]) 
# settags!(U2, "Link,l=$(idx+1)"; :tags => "qr")
# settags!(R, "Link,l=$(idx+1)"; :tags => "qr")
# mps[left_idx], mps[left_idx+1], mps[left_idx+2] = U1, U2, R
# println(inds(mps))
# println(norm(mps))
# println(get_which_canonical_form(mps), " ", linkdims(mps))

# println(linkdims(mps))
# println(get_which_canonical_form(mps))
Ho_mpo_list = get_exp_Ho_list(sites, -tau*1im/2, x)
apply_Ho_mpo_list!(Ho_mpo_list, mps; cutoff = 1e-9)

# println(linkdims(mps))
# println(get_which_canonical_form(mps))
# for element in mps
#     println(inds(element; :tags => "Link"))
# end

Hz_mpo_noprime, Hz_mpo_primed = get_exp_Hz_list(sites, a, x, l_0, mg)

# println(get_which_canonical_form(mps))

mps = apply(Hz_mpo_noprime, mps; cutoff = cut_off, normalize = true)

# println(linkdims(mps))
# println(get_which_canonical_form(mps))
# println(get_which_canonical_form(mps))

mps = apply(Hz_mpo_primed, mps; cutoff = cut_off, normalize = true)

# println(linkdims(mps))
# println(get_which_canonical_form(mps))
# println(get_which_canonical_form(mps))

# println()
# for element in mps
#     println(inds(element; :tags => "Link"))
# end

He_mpo_list = get_exp_He_list(sites, -tau*1im/2, x)
apply_He_mpo_list!(He_mpo_list, mps; cutoff = 1e-9)

# println(linkdims(mps))
# println(get_which_canonical_form(mps))
# println(get_which_canonical_form(mps))


mps = apply(Hz_mpo_noprime, mps; cutoff = cut_off, normalize = true)

# println(linkdims(mps))
# println(get_which_canonical_form(mps))
# println(get_which_canonical_form(mps))

mps = apply(Hz_mpo_primed, mps; cutoff = cut_off, normalize = true)

println(linkdims(mps))
# println(get_which_canonical_form(mps))

# H = get_Hamiltonian_noprime(sites, x, l_0, mg)
# println(inds(H))


# -- 

