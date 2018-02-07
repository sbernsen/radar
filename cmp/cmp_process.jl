#=


=#

# Load modules
using Preprocessing, radarIO

# Define filenames 
cmp_data = "jarvis.10mhz..2017.l2.cmp.longitudinal.h5"

# Read in the preprocessing steps from 'preprocessing_steps.txt'
pro_steps = read_steps()


##
trace_mat, lle = h5cmp(cmp_data) 

if pro_steps.stat_norm 
	for i in 1:size(trace_pro,2)
		trace_pro[:,i] = running_stat(trace_mat[:,i], pro_steps.kwindow, pro_steps.stat_type)
	end
	println("Statistical Normalization")
end

if pro_steps.rms_norm
	trace_pro = mapslices(rms_norm, trace_pro, 1)

	println("RMS Normalization")
end

#=
if pro_steps.rm_bgrnd
	trace_pro, data_avg = background_removal(trace_pro)

	println("Background Removal")
end
=#

if pro_steps.filter_type != "none" 
	filthy = butt_design(pro_steps.filter_type, pro_steps.fc)
	for i in 1:size(trace_pro, 2)
		trace_pro[:,i] = filt(filthy, trace_pro[:,i])
	end

	println("Filtered")
end


using RCall
z = trace_pro[:,1]
zz = trace_mat[:,1]

R"""
dev.new()
par(mfrow = c(2, 1))
plot($z, type = "l")
plot($zz, type = "l")

dev.new()
image($trace_pro)
"""