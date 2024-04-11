using ITensors
using LinearAlgebra
using Plots
using Statistics
using DataFrames
using LsqFit
include("Utilities.jl")

function get_subtracted_gap(N, x, l_0, ma, lambda, D, max_sweeps, tol, sites, state)

    H = get_aH_Hamiltonian(sites, x, l_0, ma, lambda)
    gs_initial = randomMPS(sites, state)
    sweeps = Sweeps(max_sweeps, maxdim = D)
    observer = DMRGObserver(;energy_tol = tol) # This allows to track the energy during DMRG and if the change in energy is smaller than tol it stops
    println("Starting DMRG for GS with ma = $(ma)\n")
    gs_energy, gs = dmrg(H, gs_initial, sweeps; outputlevel = 1, observer = observer, ishermitian = true) # outputlevel = 1 prints out during DMRG otherwise set to 0, here one can include cutoff = cutoff instead of maxdim = D
    first_initial = randomMPS(sites, state)
    Ms = [gs]
    w = abs(real(gs_energy))
    println("\nStarting DMRG for 1st excited state with ma = $(ma)\n")
    first_energy, first = dmrg(H, Ms, first_initial, sweeps, weight = w, ishermitian = true, observer = observer, outputlevel = 1) # here one can include cutoff = cutoff instead of maxdim = D

    return (first_energy - gs_energy)*sqrt(x) - 1/sqrt(pi)

end

function get_max_pnd(N, x, l_0, ma, tau, max_steps, cutoff, sites, mps)
    
    Ho_mpo_list = get_exp_Ho(sites, -1im*tau/2) # tau*odd/2
    Hz_mpo = get_exp_Hz(sites, -1im*tau/2, x, l_0, ma) # 1-tau*H_z/2
    He_mpo_list = get_exp_He(sites, -1im*tau) # tau*even
    Ho_mpo_list_2 = get_exp_Ho(sites, -1im*tau) # tau*odd

    PND = get_particle_number_MPO(sites)/N
    particle_number_list = [real(inner(mps', PND, mps))]

    for step in 1:max_steps
        
        if step == 1
            apply_Ho_mpo_list!(Ho_mpo_list, mps; cutoff = cutoff)
            mps = apply(Hz_mpo, mps; cutoff = cutoff)
            apply_He_mpo_list!(He_mpo_list, mps; cutoff = cutoff)
            mps = apply(Hz_mpo, mps; cutoff = cutoff)
        elseif step == max_steps
            apply_Ho_mpo_list!(Ho_mpo_list_2, mps; cutoff = cutoff)
            mps = apply(Hz_mpo, mps; cutoff = cutoff)
            apply_He_mpo_list!(He_mpo_list, mps; cutoff = cutoff)
            mps = apply(Hz_mpo, mps; cutoff = cutoff)
            apply_Ho_mpo_list!(Ho_mpo_list, mps; cutoff = cutoff)    
        else
            apply_Ho_mpo_list!(Ho_mpo_list_2, mps; cutoff = cutoff)
            mps = apply(Hz_mpo, mps; cutoff = cutoff)
            apply_He_mpo_list!(He_mpo_list, mps; cutoff = cutoff)
            mps = apply(Hz_mpo, mps; cutoff = cutoff)
        end
        
        normalize!(mps)
        pnd = real(inner(mps', PND, mps))
        print("\rStep = $(step)")
        if (pnd < particle_number_list[end]) && (pnd > 1e-15)
            println("Max PND reached")
            return particle_number_list[end]
        end
        push!(particle_number_list, pnd)

    end

    println("Max PND not reached")
    return particle_number_list[end]

end

N = 20
sites = siteinds("S=1/2", N, conserve_qns = true)
state = [isodd(n) ? "0" : "1" for n = 1:N]
dirac_vacuum = get_dirac_vacuum_mps(sites)
x = 1
l_0 = 0
lambda = 0

tau = 0.01
max_steps = 2000

D = 20
cutoff = 1e-6
tol = 1e-6
max_sweeps = 1000

ma_list = LinRange(-0.4, -0.15, 10)

let

subtracted_gap_list = Float64[]
max_pnd_list = []

for ma in ma_list

    println("Starting ma = $(ma)\n")

    push!(subtracted_gap_list, get_subtracted_gap(N, x, l_0, ma, lambda, D, max_sweeps, tol, sites, state))
    push!(max_pnd_list, get_max_pnd(N, x, l_0, ma, tau, max_steps, cutoff, sites, dirac_vacuum))

    println("\nFinished ma = $(ma)\n")

end

p1 = plot(legend = :outertopright)
plot!(p1, ma_list, max_pnd_list, label = "Max-PND", color = :blue)
scatter!(p1, ma_list, max_pnd_list, label = false, color = :blue)
plot!(p1, ma_list, subtracted_gap_list, label = "Subtracted gap", color = :black)
scatter!(p1, ma_list, subtracted_gap_list, label = false, color = :black)

model(t, p) = p[1].+p[2].*t
p0 = [1.0, 1.0]
fit = curve_fit(model, ma_list, subtracted_gap_list, p0)
ms = -fit.param[1]/fit.param[2]

hline!(p1, [0], color = :red, label = false)
# vline!(p1, [-1/(8*sqrt(x))], color = :red, label = false)
# vline!(p1, [ms], color = :green, label = "Subtracted gap fit")
vline!(p1, [ma_list[argmax(max_pnd_list)]], color = :cyan, label = "Max-PND Max")

end
