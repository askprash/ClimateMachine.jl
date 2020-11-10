using Plots
using NCDatasets
using Statistics: mean
using FFTW

using ClimateMachine.ConfigTypes

using ClimateMachine.Spectra
using ClimateMachine.Spectra:
    compute_gaussian!,
    compute_legendre!,
    SpectralSphericalMesh,
    trans_grid_to_spherical!,
    compute_wave_numbers,
    power_spectrum_2d

include("test/Common/Spectra/spherical_helper_test.jl")

# select the file and variable to output
CLIMA_NETCDF = "/central/groups/esm/lenka/"; #location of .nc files
var_name = "v"
fnames = filter(x -> occursin(".nc", x), readdir( CLIMA_NETCDF ) );

# extract data
ds = NCDataset(CLIMA_NETCDF*fnames[1], "r");
lon = ds["long"][:];
lat = ds["lat"][:];
z = ds["level"][:] / 1e3; # height in kilometers
time = ds["time"][:];
var= ds[var_name][:];# lon, lat,lev, time
var=var[:, :, 4, 29];

nan_fill_value = 0.0
replace!(var, NaN=>nan_fill_value)

# set up grid vars
nlats = length(lat)
sinθ, wts = compute_gaussian!(nlats)
cosθ = sqrt.(1 .- sinθ .^2)
mass_weight = ones(Float64, 1);

# get the specktrum
spectrum, wave_numbers, spherical, mesh = power_spectrum_2d(AtmosGCMConfigType(),var, mass_weight)

# MAGNITUDE TEST
# Spectrum
m_ = collect(0:1:mesh.num_fourier)[:]
n_ = collect(0:1:mesh.num_spherical)[:]
contourf(m_,n_, (spectrum[:,:,1])', xlabel="m", ylabel ="n" )
plot(n_[2:end], sum(spectrum[:,1:end-1,1] .+ 0.00001 ,dims=1)', xaxis=:log, yaxis=:log)

# Check global magnitude
println(sum(0.5 * spectrum))

dθ = π / length(sinθ)
EKE = 0.5 .* var[:,:,1] .^ 2
#EKE = 0.5 .* (var[:,:,1] .- mean(var[:,:,1],dims=1)) .^ 2
println(sum( EKE[:,:] .* dθ .* reshape(cosθ, (1,length(cosθ))) * dθ / 4π )) # weighted average over Earth's sfc area (need in m2/s2)

# RECONSTRUCTION TEST
# grid to spherical to grid reconstruction
reconstruction = trans_spherical_to_grid!(mesh, spherical )

# Check visually
contourf(var[:,:,1])
contourf(reconstruction[:,:,1])
contourf(var[:,:,1] .- reconstruction[:,:,1])










"""
++++++++
# this is what the FFTW.fft does
len=1000
xarray = collect(1:1:len)
x = 1.0 * sin.(xarray / xarray[end] * 5.0 * 2π) + 2.0 * sin.(xarray / xarray[end] * 20.0 * 2π) +  convert(Array{Float64,1}, rand(len) )

fourier = fft(x)

N = size(x)[1]
n_ = collect(1:1:N) .- 1
fourier_slow = zeros(ComplexF64, N)
for kK in 0:(N-1)
    M = exp.(- im * 2 * π * n_ * kK / N)
    fourier_slow[kK+1] = sum(M .* x )
end

plot(n_[:],[real(fourier .* conj(fourier)), real(fourier_slow .* conj(fourier_slow))[:] , real(fourier_slow .* conj(fourier_slow))[:] - real(fourier .* conj(fourier))[:]])
plot(n_[:],[ real(fourier_slow .* conj(fourier_slow))[:] - real(fourier .* conj(fourier))[:]])

# check mgnitude
sum(0.5 .* (fourier_slow .* conj(fourier_slow))) / N
sum(0.5 .* (x .^2))

"""
