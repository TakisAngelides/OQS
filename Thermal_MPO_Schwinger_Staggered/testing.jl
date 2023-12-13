using ITensors
using LinearAlgebra
using HDF5
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

# N = 8
# sites = siteinds("S=1/2", N)
# # state = [isodd(n) ? 1 : 2 for n=1:N]
# # mpo = get_initial_zero_charge_MPO(sites, state)
# # mpo = get_Hamiltonian(sites, 1.0, 1.0, 1.0)
# # mpo = mpo / tr(mpo)
# mpo = MPO(sites, "Id")
# mpo = mpo/tr(mpo)

# # println(commoninds(mpo[1], mpo[1]*mpo[2])...)
# # U, S, V = ITensors.svd(mpo[1]*mpo[2], commoninds(mpo[1], mpo[1]*mpo[2])..., lefttags = "Link,l=1", righttags = "Link,l=1")
# # println(U)
# # println(get_MPO_site_canonical_form(U, "first", 1))

# # @show mpo[1]
# # mpo[2] = S*V
# # println(inds(mpo))
# # normalize!(mpo)
# # orthogonalize!(mpo, 1)

# # println(get_MPO_canonical_form(mpo))
# x = 1.0
# a = 1.0im
# l_0 = 0.3
# mg = 0.2
# cutoff = 1e-9
# hj = x * op("S-", sites[1]) * op("S+", sites[2])
# Gj = exp(a * hj)
# # println(inds(Gj))

# idx = 1

# println(tr(mpo))
# println(get_MPO_canonical_form(mpo), " ", linkdims(mpo))
# Ho_mpo_list = get_exp_Ho_list(sites, a, x)
# # mpo = mpo/tr(mpo)
# apply_odd!(Ho_mpo_list, mpo)
# mpo = mpo/tr(mpo)
# println(tr(mpo))
# println(get_MPO_canonical_form(mpo), " ", linkdims(mpo))
# Hz_mpo = get_exp_Hz(sites, a, x, l_0, mg)
# println(tr(mpo))
# mpo = apply(apply(Hz_mpo, mpo), replaceprime(Hz_mpo', 2 => 0); cutoff = cutoff)
# println(linkdims(mpo))
# # mpo = mpo/tr(mpo)
# # println(tr(mpo))
# # # mpo = mpo/tr(mpo)
# # println(tr(mpo))
# # println(get_MPO_canonical_form(mpo), " ", linkdims(mpo))
# # He_mpo_list = get_exp_He_list(sites, a, x)
# # apply_He_mpo_list!(He_mpo_list, mpo)
# # println(get_MPO_canonical_form(mpo), " ", linkdims(mpo))


# # a = 1.0
# # tau = 1.0
# # x = 1.0
# # l_0 = 1.0
# # mg = 1.0
# # cut_off = 1e-16

# # gates = []
# # idx = 1
# # hj = x * op("S-", sites[idx]) * op("I", sites[idx+1]) * op("S+", sites[idx+2])
# # hj += x * op("S+", sites[idx]) * op("I", sites[idx+1]) * op("S-", sites[idx+2])
# # Gj = exp(a * hj)
# # left_idx, right_idx = idx, idx+2
# # push!(gates, (Gj, left_idx))

# # gate, left_idx = gates[1]

# # # mpo = MPO(gate, [sites[left_idx:right_idx]])
# # # println(linkdims(mpo))

# # println(inds(mps))
# # println(norm(mps))
# # println(get_which_canonical_form(mps), " ", linkdims(mps))
# # tmp = noprime(gate*mps[left_idx]*mps[left_idx+1]*mps[left_idx+2])
# # U1, S, V = ITensors.svd(tmp, commoninds(tmp, mps[left_idx]), cutoff = 1e-9, lefttags = "Link,l=$(idx)", righttags = "Link,l=$(idx)") 
# # R = S*V/norm(S)
# # println("Norm of S: ", norm(S), " ", diag(S))
# # U2, R = ITensors.qr(R, [commoninds(R, mps[left_idx + 1]), commoninds(R, U1)]) 
# # settags!(U2, "Link,l=$(idx+1)"; :tags => "qr")
# # settags!(R, "Link,l=$(idx+1)"; :tags => "qr")
# # mps[left_idx], mps[left_idx+1], mps[left_idx+2] = U1, U2, R
# # println(inds(mps))
# # println(norm(mps))
# # println(get_which_canonical_form(mps), " ", linkdims(mps))

# # println(linkdims(mps))
# # println(get_which_canonical_form(mps))

# # Ho_mpo_list = get_exp_Ho_list(sites, -tau*1im/2, x)
# # apply_Ho_mpo_list!(Ho_mpo_list, mpo; cutoff = 1e-9)

# # println(linkdims(mps))
# # println(get_which_canonical_form(mps))
# # for element in mps
# #     println(inds(element; :tags => "Link"))
# # end

# # Hz_mpo_noprime, Hz_mpo_primed = get_exp_Hz_list(sites, a, x, l_0, mg)

# # println(get_which_canonical_form(mps))

# # mps = apply(Hz_mpo_noprime, mps; cutoff = cut_off, normalize = true)

# # println(linkdims(mps))
# # println(get_which_canonical_form(mps))
# # println(get_which_canonical_form(mps))

# # mps = apply(Hz_mpo_primed, mps; cutoff = cut_off, normalize = true)

# # println(linkdims(mps))
# # println(get_which_canonical_form(mps))
# # println(get_which_canonical_form(mps))

# # println()
# # for element in mps
# #     println(inds(element; :tags => "Link"))
# # end

# # He_mpo_list = get_exp_He_list(sites, -tau*1im/2, x)
# # apply_He_mpo_list!(He_mpo_list, mps; cutoff = 1e-9)

# # println(linkdims(mps))
# # println(get_which_canonical_form(mps))
# # println(get_which_canonical_form(mps))


# # mps = apply(Hz_mpo_noprime, mps; cutoff = cut_off, normalize = true)

# # println(linkdims(mps))
# # println(get_which_canonical_form(mps))
# # println(get_which_canonical_form(mps))

# # mps = apply(Hz_mpo_primed, mps; cutoff = cut_off, normalize = true)

# # println(linkdims(mps))
# # println(get_which_canonical_form(mps))

# # H = get_Hamiltonian_noprime(sites, x, l_0, mg)
# # println(inds(H))


# -- 

# Prepare the ground state at l_0 = 0 as a rho MPO - this just needs the outer function of ITensors

# n = 6
# tol = 1e-14
# x = 1.0
# l_0 = 0.0
# mg = 0.1
# max_steps = 100
# cutoff = 1e-14
# sites = siteinds("S=1/2", n, conserve_qns = true) 
# H = get_Hamiltonian(sites, x, l_0, mg)
# state = [isodd(i) ? "0" : "1" for i = 1:n]
# psi0 = randomMPS(sites, state)
# obs = DMRGObserver(;energy_tol = tol)
# nsweeps = max_steps
# dmrg_energy, dmrg_state = dmrg(H, psi0; nsweeps, cutoff, observer=obs, outputlevel=1)

# mpo = outer(dmrg_state', dmrg_state; cutoff = 1e-20)
# println(inds(mpo))
# println(linkdims(dmrg_state))
# println(linkdims(mpo))
# println(dmrg_energy)
# println(inner(H, mpo))

# -- 

# let

#     N = 4
#     tau = 0.01
#     x = 1.0
#     l_0 = 0.2
#     mg = 0.3
#     cutoff = 1e-20
#     max_steps = 1000

#     sites = siteinds("S=1/2", N, conserve_qns = true) 
#     H = get_Hamiltonian(sites, x, l_0, mg)
#     state = [isodd(n) ? "0" : "1" for n = 1:N]
#     psi0 = randomMPS(sites, state, linkdims = 2)
#     obs = DMRGObserver(;energy_tol = tol)
#     nsweeps = max_steps
#     dmrg_energy, dmrg_state = dmrg(H, psi0; nsweeps, cutoff, observer=obs, outputlevel=1)

#     sites = siteinds("S=1/2", N, conserve_qns = true)
#     rho = MPO(sites, "Id")
#     orthogonalize!(rho, 1)
#     rho = rho/tr(rho)

#     odd_gates_4 = get_exp_Ho_list(sites, -tau/4, x) # odd/4
#     exp_Hz_mpo = get_exp_Hz(sites, -tau/4, x, l_0, mg) # 1+aH_z/4
#     even_gates_2 = get_exp_He_list(sites, -tau/2, x) # even/2
#     odd_gates_2 = get_exp_Ho_list(sites, -tau/2, x) # odd/2
#     H = get_Hamiltonian(sites, x, l_0, mg)

#     for step in 1:max_steps

#         println("Step is $step")
        
#         if step == 1

#             println(get_MPO_canonical_form(rho), " ", linkdims(rho), " ", inner(H, rho))
#             apply_odd!(odd_gates_4, rho; cutoff = cutoff)
#             rho = rho/tr(rho)
#             println(get_MPO_canonical_form(rho), " ", linkdims(rho), " ", inner(H, rho))
#             rho = apply(apply(exp_Hz_mpo, rho), replaceprime(dag(exp_Hz_mpo'), 2 => 0); cutoff = cutoff)
#             rho = rho/tr(rho)
#             println(get_MPO_canonical_form(rho), " ", linkdims(rho), " ", inner(H, rho))
#             apply_even!(even_gates_2, rho; cutoff = cutoff)
#             rho = rho/tr(rho)
#             println(get_MPO_canonical_form(rho), " ", linkdims(rho), " ", inner(H, rho))
#             rho = apply(apply(exp_Hz_mpo, rho), replaceprime(dag(exp_Hz_mpo'), 2 => 0); cutoff = cutoff)
#             rho = rho/tr(rho)
            
#         elseif step == max_steps
        
#             println(get_MPO_canonical_form(rho), " ", linkdims(rho), " ", inner(H, rho))
#             apply_odd!(odd_gates_2, rho; cutoff = cutoff)
#             rho = rho/tr(rho)
#             println(get_MPO_canonical_form(rho), " ", linkdims(rho), " ", inner(H, rho))
#             rho = apply(apply(exp_Hz_mpo, rho), replaceprime(dag(exp_Hz_mpo'), 2 => 0); cutoff = cutoff)
#             rho = rho/tr(rho)
#             println(get_MPO_canonical_form(rho), " ", linkdims(rho), " ", inner(H, rho))
#             apply_even!(even_gates_2, rho; cutoff = cutoff)
#             rho = rho/tr(rho)
#             println(get_MPO_canonical_form(rho), " ", linkdims(rho), " ", inner(H, rho))
#             rho = apply(apply(exp_Hz_mpo, rho), replaceprime(dag(exp_Hz_mpo'), 2 => 0); cutoff = cutoff)
#             rho = rho/tr(rho)
#             println(get_MPO_canonical_form(rho), " ", linkdims(rho), " ", inner(H, rho))
#             apply_odd!(odd_gates_4, rho; cutoff = cutoff)
#             rho = rho/tr(rho)
#             println(get_MPO_canonical_form(rho), " ", linkdims(rho), " ", inner(H, rho))
                
#         else

#             println(get_MPO_canonical_form(rho), " ", linkdims(rho), " ", inner(H, rho))
#             apply_odd!(odd_gates_2, rho; cutoff = cutoff)
#             rho = rho/tr(rho)
#             println(get_MPO_canonical_form(rho), " ", linkdims(rho), " ", inner(H, rho))
#             rho = apply(apply(exp_Hz_mpo, rho), replaceprime(dag(exp_Hz_mpo'), 2 => 0); cutoff = cutoff)
#             rho = rho/tr(rho)
#             println(get_MPO_canonical_form(rho), " ", linkdims(rho), " ", inner(H, rho))
#             apply_even!(even_gates_2, rho; cutoff = cutoff)
#             rho = rho/tr(rho)
#             println(get_MPO_canonical_form(rho), " ", linkdims(rho), " ", inner(H, rho))
#             rho = apply(apply(exp_Hz_mpo, rho), replaceprime(dag(exp_Hz_mpo'), 2 => 0); cutoff = cutoff)
#             rho = rho/tr(rho)

#         end

#     end

#     # println(inds(rho))
#     println(dmrg_energy)

# end

# -- 

# f = h5open("test1.h5", "w")
# sites = siteinds("S=1/2", 4)
# println(sites)
# mpo = randomMPO(sites)
# write(f, "mpo", mpo)
# close(f)
# f = h5open("test1.h5", "r")
# mpo_read = read(f, "mpo", MPO)
# sites_read = siteinds(mpo_read; :plev => 0)
# println("----------")
# println(reduce(vcat, siteinds(mpo_read; :plev => 0)))
# close(f)

# -- 
