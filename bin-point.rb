module BIN_POINT
	class Bin
		attr_accessor :points, :min_int, :max_int, :sum_intensity
		@@num = 0

		def initialize(point=nil)
			@@num += 1
			@points = []
			@count = 0
			@sum = 0
			@min_int = nil
			@max_int = 0.0
			@sum_intensity = 0.0
			self.add(point) if point
		end

		def avg
			return @sum / @count
		end

		def add(point)
			@points << point
			@sum += point.mz
			@sum_intensity += point.intensity
			@count += 1
			@min_int ||= point.intensity
			@min_int = point.intensity if point.intensity < @min_int
			@max_int = point.intensity if point.intensity > @max_int
		end

		def std_dev
			sum_square_diff_mean = 0.0
			avg = self.avg
			@points.collect{|pt| sum_square_diff_mean += (pt.mz - avg)**2}
			return Math.sqrt(sum_square_diff_mean / @points.size)
		end

		def combine(bin)
			bin.points.each do |pt|
				self.add(pt)
			end
		end

		def size
			return @points.size
		end
	end

	class Point
		attr_accessor :mz, :intensity, :rt, :cluster_id, :true_cluster_id
	
		def initialize(mz, intensity, rt, cluster_id = nil)
			@mz = mz
			@intensity = intensity
			@rt = rt
			@cluster_id = cluster_id
			@true_cluster_id = nil
		end
	end
end
