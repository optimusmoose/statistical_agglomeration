################################################################################
#                      Statiscal Agglomeration
# This algorithm reduces a list of three dimensional direct injection
#     lipidomics output points (m/z, intensity, retention time) to a list
#     of two dimensional centroids (m/z, intensity). 
#
# USAGE: "ruby sa.rb filename" where filename contains the mass spec input
#     as a .csv with one line for each point, and each point in the format
#     <m/z>, <intensity>, <RT>. The default output is a list of centroided
#     peaks. That will produce a .csv where each line represents a centroid, 
#     with the first value being the m/z of that peak and the second value is 
#     the intensity of the peak. Alternatively, you can get each point by 
#     assignment by using: "ruby sa.rb filename verbose" This will produce
#     a .csv with one line per inputted point in the same format as the 
#     input file with the addition of a peak ID at the end of the line:
#     <m/z>, <intensity>, <RT>, <peak ID>. The peak ID "0" is reserved 
#     for noise points.
#
# REQUIREMENTS: None.
#
# LICENSE: Copyright 2012 Rob Smith. Free for non-commercial use provided this
#     header is included and the corresponding paper is cited in any resulting 
#     publication. Commercial usage must be licensed through the intellectual
#     property office at Brigham Young University. Email 2robsmith@gmail.com for
#     more information.
#
# NOTE ON FILE FORMATS: If you need to convert your mzML file to the csv format
#     used here, I recommend installing Mspire "gem install mspire". Run
#     "ruby mzML_to_csv.rb" included in this repository. You'll get a .csv file
#     in the right format.
################################################################################

require_relative 'bin-point'

#2 tail t-test. P vals are 0.1	0.05	0.02	0.01
$t_table = {2=>[2.92, 4.3027, 6.2054, 9.925], 3=>[2.3534, 3.1824, 4.1765, 5.8408], 4=>[2.1318, 2.7765, 3.4954, 4.6041], 5=>[2.015, 2.5706, 3.1634, 4.0321], 6=>[1.9432, 2.4469, 2.9687, 3.7074], 7=>[1.8946, 2.3646, 2.8412, 3.4995], 8=>[1.8595, 2.306, 2.7515, 3.3554], 9=>[1.8331, 2.2622, 2.685, 3.2498], 10=>[1.8125, 2.2281, 2.6338, 3.1693], 11=>[1.7959, 2.201, 2.5931, 3.1058], 12=>[1.7823, 2.1788, 2.56, 3.0545], 13=>[1.7709, 2.1604, 2.5326, 3.0123], 14=>[1.7613, 2.1448, 2.5096, 2.9768], 15=>[1.7531, 2.1315, 2.4899, 2.9467], 16=>[1.7459, 2.1199, 2.4729, 2.9208], 17=>[1.7396, 2.1098, 2.4581, 2.8982], 18=>[1.7341, 2.1009, 2.445, 2.8784], 19=>[1.7291, 2.093, 2.4334, 2.8609], 20=>[1.7247, 2.086, 2.4231, 2.8453], 21=>[1.7207, 2.0796, 2.4138, 2.8314], 22=>[1.7171, 2.0739, 2.4055, 2.8188], 23=>[1.7139, 2.0687, 2.3979, 2.8073], 24=>[1.7109, 2.0639, 2.391, 2.797], 25=>[1.7081, 2.0595, 2.3846, 2.7874], 26=>[1.7056, 2.0555, 2.3788, 2.7787], 27=>[1.7033, 2.0518, 2.3734, 2.7707], 28=>[1.7011, 2.0484, 2.3685, 2.7633], 29=>[1.6991, 2.0452, 2.3638, 2.7564], 30=>[1.6973, 2.0423, 2.3596, 2.75], 31=>[1.6955, 2.0395, 2.3556, 2.744], 32=>[1.6939, 2.0369, 2.3518, 2.7385], 33=>[1.6924, 2.0345, 2.3483, 2.7333], 34=>[1.6909, 2.0322, 2.3451, 2.7284], 35=>[1.6896, 2.0301, 2.342, 2.7238], 36=>[1.6883, 2.0281, 2.3391, 2.7195], 37=>[1.6871, 2.0262, 2.3363, 2.7154], 38=>[1.686, 2.0244, 2.3337, 2.7116], 39=>[1.6849, 2.0227, 2.3313, 2.7079], 40=>[1.6839, 2.0211, 2.3289, 2.7045], 41=>[1.6829, 2.0195, 2.3267, 2.7012], 42=>[1.682, 2.0181, 2.3246, 2.6981], 43=>[1.6811, 2.0167, 2.3226, 2.6951], 44=>[1.6802, 2.0154, 2.3207, 2.6923], 45=>[1.6794, 2.0141, 2.3189, 2.6896], 46=>[1.6787, 2.0129, 2.3172, 2.687], 47=>[1.6779, 2.0117, 2.3155, 2.6846], 48=>[1.6772, 2.0106, 2.3139, 2.6822], 49=>[1.6766, 2.0096, 2.3124, 2.68], 50=>[1.6759, 2.0086, 2.3109, 2.6778], 51=>[1.6753, 2.0076, 2.3095, 2.6757], 52=>[1.6747, 2.0066, 2.3082, 2.6737], 53=>[1.6741, 2.0057, 2.3069, 2.6718], 54=>[1.6736, 2.0049, 2.3056, 2.67], 55=>[1.673, 2.004, 2.3044, 2.6682], 56=>[1.6725, 2.0032, 2.3033, 2.6665], 57=>[1.672, 2.0025, 2.3022, 2.6649], 58=>[1.6716, 2.0017, 2.3011, 2.6633], 59=>[1.6711, 2.001, 2.3, 2.6618], 60=>[1.6706, 2.0003, 2.299, 2.6603], 61=>[1.6702, 1.9996, 2.2981, 2.6589], 62=>[1.6698, 1.999, 2.2971, 2.6575], 63=>[1.6694, 1.9983, 2.2962, 2.6561], 64=>[1.669, 1.9977, 2.2954, 2.6549], 65=>[1.6686, 1.9971, 2.2945, 2.6536], 66=>[1.6683, 1.9966, 2.2937, 2.6524], 67=>[1.6679, 1.996, 2.2929, 2.6512], 68=>[1.6676, 1.9955, 2.2921, 2.6501], 69=>[1.6672, 1.9949, 2.2914, 2.649], 70=>[1.6669, 1.9944, 2.2906, 2.6479], 71=>[1.6666, 1.9939, 2.2899, 2.6469], 72=>[1.6663, 1.9935, 2.2892, 2.6458], 73=>[1.666, 1.993, 2.2886, 2.6449], 74=>[1.6657, 1.9925, 2.2879, 2.6439], 75=>[1.6654, 1.9921, 2.2873, 2.643], 76=>[1.6652, 1.9917, 2.2867, 2.6421], 77=>[1.6649, 1.9913, 2.2861, 2.6412], 78=>[1.6646, 1.9908, 2.2855, 2.6403], 79=>[1.6644, 1.9905, 2.2849, 2.6395], 80=>[1.6641, 1.9901, 2.2844, 2.6387], 81=>[1.6639, 1.9897, 2.2838, 2.6379], 82=>[1.6636, 1.9893, 2.2833, 2.6371], 83=>[1.6634, 1.989, 2.2828, 2.6364], 84=>[1.6632, 1.9886, 2.2823, 2.6356], 85=>[1.663, 1.9883, 2.2818, 2.6349], 86=>[1.6628, 1.9879, 2.2813, 2.6342], 87=>[1.6626, 1.9876, 2.2809, 2.6335], 88=>[1.6624, 1.9873, 2.2804, 2.6329], 89=>[1.6622, 1.987, 2.28, 2.6322], 90=>[1.662, 1.9867, 2.2795, 2.6316], 91=>[1.6618, 1.9864, 2.2791, 2.6309], 92=>[1.6616, 1.9861, 2.2787, 2.6303], 93=>[1.6614, 1.9858, 2.2783, 2.6297], 94=>[1.6612, 1.9855, 2.2779, 2.6291], 95=>[1.6611, 1.9852, 2.2775, 2.6286], 96=>[1.6609, 1.985, 2.2771, 2.628], 97=>[1.6607, 1.9847, 2.2767, 2.6275], 98=>[1.6606, 1.9845, 2.2764, 2.6269], 99=>[1.6604, 1.9842, 2.276, 2.6264], 100=>[1.6602, 1.984, 2.2757, 2.6259], :inf=>[1.645, 1.960, 2.326, 2.576]}

def old_t_test(pop_mean, sample_mean, sample_size, sample_std_dev)
	#one-sample t-test
	t = (sample_mean.to_f - pop_mean)/(sample_std_dev / Math.sqrt(sample_size))
	df = sample_size-1

	t_table_row = []
	if df > 100
		t_table_row = $t_table[:inf]
	elsif df > 1
		t_table_row = $t_table[df]
	else
		return 1
	end

	p = 0
	t_table_row.each_with_index do |t_lookup, p_idx|
		p = p_idx + 1 if t >= t_lookup
	end
	p_s = [1,0.1,0.05,0.02,0.01]
	return p_s[p]
end

def variance_mean(data)
	mean = data.inject(0,:+) / data.size.to_f
	sum = 0
	data.each do |datum|
		sum += (datum - mean)**2
	end
	return [sum / data.size.to_f, mean]
end

#http://en.wikipedia.org/wiki/Welch%27s_t_test
def t_test(n1,d1,n2,d2)
	s1, mu1 = variance_mean(d1)
	s2, mu2 = variance_mean(d2)
	term = s1/d1.size.to_f + s2/d2.size.to_f # s1_over_n1_p_s2_over_n2

	t = (mu1 - mu2) / (Math.sqrt(term))
	zero_check = true if (n1 == 1 or n2 == 1)
	
	df = zero_check ? 0 : (term**2 / ((s1**2 / (n1**2 * (n1 - 1))) + (s2**2 / n2**2 * (n2 - 1)))).to_i
	t_table_row = []
	if df > 100
		t_table_row = $t_table[:inf]
	elsif df > 1
		t_table_row = $t_table[df]
	else
		return 0
	end

	p = 0
	t_table_row.each_with_index do |t_lookup, p_idx|
		p = p_idx + 1 if t >= t_lookup
	end
	p_s = [1,0.1,0.05,0.02,0.01]
	return p_s[p]
end

def t_test_helper(data1, data2)
	n1, d1, n2, d2 = nil
	if data1.size > data2.size
		n1 = data1.size
		d1 = data1.points
		n2 = data2.size
		d2 = data2.points
	else
		n1 = data2.size
		d1 = data2.points
		n2 = data1.size
		d2 = data1.points
	end
	p_mz = t_test(n1,d1.map{|point| point.mz},n2,d2.map{|point| point.mz})
	p_int = t_test(n1,d1.map{|point| point.intensity},n2,d2.map{|point| point.intensity})

	return [p_mz, p_int].max
end

def print_points(bins, filename)
	outfile = File.open(filename,"w")
	bins.each do |b|
		b.points.each do |pt|
			outfile.puts "#{pt.mz},#{pt.intensity},#{pt.rt},#{pt.cluster_id}"
		end
	end
	outfile.close
end

def print_bins(bins, filename)
	outfile = File.open(filename,"w")
	bins.each do |b|
		#mz centroid is calculated as an average of the point's mz's weighted
		#by each point's share of the intensity of the peak it belongs to
		mz = 0.0
		b.points.collect{|pt| mz += pt.mz * (pt.intensity/b.sum_intensity)}

		intensity = b.sum_intensity / b.points.size # average intensity

		outfile.puts "#{mz},#{intensity}"
	end
	outfile.close
end



#I/O
data_points = []

file = File.open(ARGV[0], "r")
while line = file.gets
	parts = line.chop.split(",")
	data_points << BIN_POINT::Point.new(parts[0].to_f, parts[1].to_f, parts[2].to_f)
end
file.close

data_points.sort!{|x,y| x.mz<=>y.mz}

dec_place = 100
last_mz = (data_points[0].mz * dec_place).to_i
bins = []
curr_bin = BIN_POINT::Bin.new
data_points.each do |datum|
	if (datum.mz * dec_place).to_i == last_mz
		curr_bin.add(datum)
	else
		bins << curr_bin
		curr_bin = BIN_POINT::Bin.new
		curr_bin.add(datum)
		last_mz = (datum.mz * dec_place).to_i
	end
end
bins << curr_bin

#Clean up the moved points
delete_list = []
(bins.size - 2).times do |b_idx|
	if bins[b_idx].size > 0 and bins[b_idx+1].size > 0
		if t_test_helper(bins[b_idx], bins[b_idx+1]) > 0.01
			bins[b_idx+1].combine(bins[b_idx])
			delete_list << b_idx
		end
	end
end
delete_list.reverse!
delete_list.each do |d_idx|
	bins.delete_at(d_idx)
end

#clean up the peak indices (no longer accurate due to agglomerations)
bctr = 1
bindex = 0
bins.each do |bin|
	if bin.size <= 1 #label all peaks of size less than 1 as noise
		bindex = 0
	else
		bindex = bctr
	end
	bin.points.each do |pt|
		pt.cluster_id = bindex
	end
	bctr += 1
end

if ARGV[1]=="verbose" #print each point and its peak assignment
	print_points(bins, "#{ARGV[0].gsub('.csv', '_sa.txt')}")
else 
	print_bins(bins, "#{ARGV[0].gsub('.csv', '_sa.txt')}")
end
