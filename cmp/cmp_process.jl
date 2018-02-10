#=


=#

# Load modules
using Preprocessing, radarIO

# Define filenames 
cmp_data = "jarvis.10mhz..2017.l2.cmp.transverse.h5"

# The labeling of the file is unique so edit the following couple lines accordingly
title_id = split(cmp_data, ".")[end-3:end-1]
title_id =join(title_id, " - ")

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



#=

###############################################################################
# Since there is ringing in the signal let's create the time domain filter using the trace from the reciever at the furthest location because that is where the ringing becomes dominant 
ts = trace_pro[:,end]
z = [zeros(m); ts ]
Z = fft(z)#.*conj(fft(z) )

omega = collect(0:(length(Z) -1))./(4*m)

# find first and second maxima of the spectra
p = sortperm(abs( Z ) )

max_indices = p[(end-8):end]
Z_max = abs(Z[max_indices])

# Get the frequencies of those spectra
omega_maxima = omega[max_indices] 


# Zero out those frequencies then inverse ifft 
# pad then compute fft
for i in 1:n
	F = (fft( [zeros(m); trace_pro[:,i]] ) )
	F[ max_indices ] = 0 
	trace_pro[:,i] = real( ifft(F) )[(m+1):end]
end




# Let's do some basic fourier theory 
t = collect(1:length(omega) )
Y = sin(2*pi*omega[max1ind].*t) + coef.*sin(2*pi*omega[max2ind].*t)

# get the time lags
xcf = crosscor( Y, ts, collect( -( length(Y)-1):( length(Y)-1) ) )
lag = find(maximum(abs(xcf) ) .== abs(xcf) )
Y = Y.*sign( xcf[lag] )
tshift = (lag-m)[1]

# Create the time domain filter
time_domain_filter = sin(2*pi*omega_maxima[1].*(t+tshift) ) + coef.*sin(2*pi*omega_maxima[2].*(t+tshift) )
time_domain_filter = time_domain_filter./maximum(time_domain_filter)


# Remove the ringing from the traces 
multipliers = mapslices(maximum, trace_pro, 1) 

for i in 1:n
	# normalize the trace to 1 and filter 

	trace_pro[:,i] = multipliers[i].*( ( trace_pro[:,i]./multipliers[i] ) - time_domain_filter) 
end




R"""

	par( mfrow = c(2, 1), mai = c(1, 1, 0.1, 0.1)  )
	library(Hmisc)
	plot( $omega, Re( $Z ), type = "l", xlab = "Frequency (Hz)", ylab = "Amplitude" )
	minor.tick(nx = 5)
	points( $omega_maxima, $Z_maxima, col = "red", lwd = 1.5 )

	plot( $t, $time_domain_filter, lwd = 0.5, type = "l", col = "red", xaxt = "n", yaxt = "n")
	lines( $t, $ts/max( $ts ), lwd = 1.5 )
"""
###############################################################################

=#

if pro_steps.rm_bgrnd
	trace_pro, data_avg = background_removal(trace_pro)

	println("Background Removal")
end


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


#=
# Apply time migration to account for reciever/transmitter spacing and zero all values before
t_adj = zeros(n, 1)

for i in 1:n
	trace_pro[:,i], t_adj[i] = time_adjust(trace_pro[:,i], x_vec[i], pro_steps.fs, v)
end
println("Time Adjusted")

# In the blue ice radar, the air wave is already adjusted for so we need to compensate for this. For example, the antenna spacing was input into the metadata to approximate the time zero for the time series.
air_wave_compensation = 62.5/v
t_adj_time = t_adj./pro_steps.fs

t_adj_int = Int.(round(t_adj_time.*pro_steps.fs) )

for i in 1:n
	trace_pro[1:t_adj_int[i],i ] = 0
end

=#
# remove the air_wave_compensation
trace_pro[ Int( round(air_wave_compensation*pro_steps.fs) ):end, :]


###### Take a look at the results

using RCall
#=
R"""

"""
=#



R"""
M = 500
time_vec=$time_vec[1:M]
x_vec = $x_vec[1:($n-8)]
tp = $trace_pro[1:M, 1:($n-8)]

#col_pal = grey( seq(1, 0, length = 256) ) 
#dev.new(width = 5, height = 12, units = "in") 
#image(x_vec, time_vec, t(tp), col = col_pal, , ylim = c( max(time_vec), min(time_vec) ), xlab = "Rx-Tx Distance", ylab = "Two Way Travel Time", main = "T2-CMP" )


png(paste( "~/radar/cmp/", $title_id, ".png", sep = ""), height = 12, width = 3, units = "in", res = 100)
library(RSEIS)
library(caTools)
source("~/.R/r_modules/plotting_functions.R")
plot_geophone_wiggles( time_vec, -tp, c('blue', 'green'), TRUE, c(12, 3) )

axis(3, at = c( 1:length(x_vec) ), labels = x_vec, las = 3, cex.axis = 0.5, lwd = 2)
mtext($title_id, 1)
dev.off()
"""
