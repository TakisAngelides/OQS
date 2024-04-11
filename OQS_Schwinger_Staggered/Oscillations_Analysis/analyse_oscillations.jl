using ITensors
using SparseArrays
using TupleTools
using Arpack
using Plots
using LinearAlgebra
using KrylovKit
using OpenQuantumTools
using HDF5
using LaTeXStrings
using Statistics
using LsqFit
include("/Users/takisangelides/Documents/PhD/Project_3_OQS/OQS/OQS_Schwinger_Staggered/Utilities.jl")
include("/Users/takisangelides/Documents/PhD/Project_3_OQS/OQS/Closed_Schwinger_Staggered/Utilities.jl")
ITensors.disable_warn_order()

function get_peaks(pnd)

    peak_indices = []
    upward = pnd[1] < pnd[2]
    for (idx, element) in enumerate(pnd)
        if (idx != length(pnd)) && (idx != 1)
            if ((element < pnd[idx+1]) && (element > pnd[idx-1]))
                upward = true
            end
            if (element > pnd[idx+1]) && (element < pnd[idx-1]) && upward
                push!(peak_indices, idx)
                upward = false
            end
        end
    end
    pnd_peaks = pnd[peak_indices]

    return peak_indices, pnd_peaks

end

function get_troughs(pnd)

    trough_indices = []
    downward = pnd[1] > pnd[2]
    for (idx, element) in enumerate(pnd)
        if (idx != length(pnd)) && (idx != 1)
            if ((element > pnd[idx+1]) && (element < pnd[idx-1]))
                downward = true
            end
            if (element < pnd[idx+1]) && (element > pnd[idx-1]) && downward
                push!(trough_indices, idx)
                downward = false
            end
        end
    end
    pnd_troughs = pnd[trough_indices]

    return trough_indices, pnd_troughs

end

let

N = 6
rhodim = binomial(N, div(N, 2))
e = 1
x = 1/(e^2)
volume = N/sqrt(x)
ma = 0.1
l_0 = 10
lambda = 0
env_corr_type = "delta"
sigma_over_a = 3.0
aD_0 = 0.1
beta = 1/0.001
aT = 1/beta
dt = 0.01
steps = 50000
at_list = dt*collect(0:steps)

L = get_Lindblad_reduced_sparse_matrix(N, x, ma, l_0, lambda, aD_0, sigma_over_a, aT, env_corr_type)
evolution_operator = exp(Matrix(L)*dt)

max_sweeps = 50
tol = 1e-16
D = 1000
sites = siteinds("S=1/2", N, conserve_qns = true)
state = [isodd(n) ? "0" : "1" for n = 1:N]
mps = randomMPS(sites, state)
H = get_aH_Hamiltonian(sites, x, 0, ma, lambda)
sweeps = Sweeps(max_sweeps, maxdim = D)
observer = DMRGObserver(;energy_tol = tol)
gs_energy, gs = dmrg(H, mps, sweeps; outputlevel = 0, observer = observer, ishermitian = true) 
rho = outer(gs', gs; cutoff = 0)
rho = project_zeroq(mpo_to_matrix(rho))

rho = rho/tr(rho)

# rho = get_dirac_vacuum_zeroq_density_matrix_sparse(N)

# L_tmp = get_Lindblad_reduced_sparse_matrix(N, x, ma, 0, lambda, aD_0, sigma_over_a, aT, env_corr_type)
# evolution_operator_tmp = exp(Matrix(L_tmp)*1000)
# rho_v_tmp = reshape(rho, rhodim*rhodim)
# for step in 1:steps
#     rho_v_tmp = evolution_operator_tmp*rho_v_tmp
# end
# rho = reshape(rho_v_tmp, rhodim, rhodim)

rho_v = reshape(rho, rhodim*rhodim)

pnd_op = project_zeroq(get_particle_number_operator_sparse(N))/N
pnd = [real(tr(rho * pnd_op))]

# charge_config = get_charge_config_from_zeroq_density_matrix(N, rho)
# println(charge_config)

for step in 1:steps

    rho_v = evolution_operator*rho_v
    rho = reshape(rho_v, rhodim, rhodim)
    push!(pnd, real(tr(rho * pnd_op)))

end

# charge_config = get_charge_config_from_zeroq_density_matrix(N, rho)
# println(charge_config)

p1 = plot(at_list, pnd, title = "N=$(N),aD=$(aD_0),ma=$(ma),l_0=$(l_0),beta=$(beta),tau=$(dt),e=$(e)", size=(800, 600), label = false)
xlabel!(p1, "at")
ylabel!(p1, "Particle number / N")

f = h5open("OQS_Schwinger_Staggered/Oscillations_Analysis/at_vs_pnd_$(aT).h5", "w")
write(f, "at_list", at_list)
write(f, "pnd", pnd)
close(f)

try

    peak_indices, pnd_peaks = get_peaks(pnd)
    trough_indices, pnd_troughs = get_troughs(pnd)

    if length(pnd_peaks) > length(pnd_troughs)
        peak_indices = peak_indices[1:length(pnd_troughs)]
        pnd_peaks = pnd_peaks[1:length(pnd_troughs)]
    else
        trough_indices = trough_indices[1:length(pnd_peaks)]
        pnd_troughs = pnd_troughs[1:length(pnd_peaks)]
    end

    at_peaks, at_troughs = at_list[peak_indices], at_list[trough_indices]
    periods_from_peaks = [at_peaks[idx]-at_peaks[idx-1] for idx in 2:length(at_peaks)]
    periods_from_troughs = [at_troughs[idx]-at_troughs[idx-1] for idx in 2:length(at_troughs)]
    scatter!(p1, at_peaks, pnd_peaks, label = "Peaks: Avg P = $(round(mean(periods_from_peaks), digits = 4))")
    scatter!(p1, at_troughs, pnd_troughs, label = "Troughs: Avg P = $(round(mean(periods_from_troughs), digits = 4))")

    oscillation_amplitudes = pnd_peaks - pnd_troughs
    at_oscillations = (at_peaks + at_troughs)/2
    p2 = scatter(at_oscillations, oscillation_amplitudes, title = "N=$(N),aD=$(aD_0),ma=$(ma),l_0=$(l_0),beta=$(beta),tau=$(dt),e=$(e)", label = false, size=(800, 600))
    xlabel!(p2, "Oscillation number")
    ylabel!(p2, "Oscillation amplitude")

    model(x, p) = p[1] * exp.(-p[2] .* x)
    p0 = [1.0, 0.1]
    fit = curve_fit(model, at_oscillations, oscillation_amplitudes, p0)
    params = coef(fit)
    plot!(p2, at_oscillations, model(at_oscillations, params), label = "Fitted exponential")
    display(p2)

catch ex

    println("Plot of oscillation amplitudes gave an exception.")

end

display(p1)

end
