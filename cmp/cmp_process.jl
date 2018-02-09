#=


=#

# Load modules
using Preprocessing, radarIO

# Define filenames 
cmp_data = "jarvis.10mhz..2017.l2.cmp.transverse.h5"
v = 299792458 # m/s; speed of light

# Read in the preprocessing steps from 'preprocessing_steps.txt'
pro_steps = read_steps()


##
trace_mat, lle = h5cmp(cmp_data) 
m = size(trace_mat, 1)
n = size(trace_mat, 2) 


time_vec = ( collect(1:m ) - 1)./pro_steps.fs
x_vec = pro_steps.dxi + ( collect(1:n) - 1 ).*pro_steps.d_rxtx 

trace_pro = trace_mat



if pro_steps.filter_type != "none" 
	filthy = butt_design(pro_steps.filter_type, pro_steps.fc)
	for i in 1:n
		trace_pro[:,i] = filt(filthy, trace_pro[:,i])
	end

	println("Filtered")
end


if pro_steps.rms_norm
	trace_pro = mapslices(rms_norm, trace_pro, 1)

	println("RMS Normalization")
end


if pro_steps.stat_norm 
	for i in 1:n
		trace_pro[:,i] = running_stat(trace_pro[:,i], pro_steps.kwindow, pro_steps.stat_type)
	end
	println("Statistical Normalization")
end


if pro_steps.rm_bgrnd
	trace_pro, data_avg = background_removal(trace_pro)

	println("Background Removal")
end


# Time domain filter 

z = trace_pro[:, end] 
Z = fft(z).*conj(fft(z) )

omega = collect(0:(1/(length(z) -1)  ):1)

# find maxima
ind_max = find(real(Z) .== maximum( real(Z) ) )

R"""
	library(Hmisc)
	plot($omega[1:(length($z)/2)], Re($Z[1:(length($z)/2)]), type = "l")
	minor.tick(nx = 5)
	points($omega[$ind_max], Re($Z[$ind_max]))
"""

#=
# Apply time migration to account for reciever/transmitter spacing and zero all values before
t_adj = zeros(n, 1)

for i in 1:n
	trace_pro[:,i], t_adj[i] = time_adjust(trace_pro[:,i], x_vec[i], pro_steps.fs, v)
end
println("Time Adjusted")

# In the blue ice radar, the air wave is already adjusted for so we need to compensate for this. For example, the antenna spacing was input into the metadata to approximate the time zero for the time series.
air_wave_compensation = 62.5/v
t_adj_time = t_adj./pro_steps.fs + air_wave_compensation

t_adj_int = Int.(round(t_adj_time.*pro_steps.fs) )

for i in 1:n
	trace_pro[1:t_adj_int[i],i ] = 0
end


###### Take a look at the results

using RCall

R"""
library(RSEIS)
library(caTools)
source("~/.R/r_modules/plotting_functions.R")
plot_geophone_wiggles( $time_vec[1:800], $trace_pro[1:800,], c('red', 'blue'), FALSE, $t_adjust)
"""

=#


