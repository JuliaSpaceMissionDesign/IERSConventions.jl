

const COEFFS_CIO2000_SP = SVector(94.00, 3808.35, -119.94, -72574.09, 27.70, 15.61)

const COEFFS_CIO2000_S = [

	# j = 0
	[
		# 1-10
		PoissonSeries(-2640.73, 0.39, SVector(0, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-63.53, 0.02, SVector(0, 0, 0, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-11.75, -0.01, SVector(0, 0, 2, -2, 3), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-11.21, -0.01, SVector(0, 0, 2, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(4.57, 0.0, SVector(0, 0, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-2.02, 0.0, SVector(0, 0, 2, 0, 3), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1.98, 0.0, SVector(0, 0, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(1.72, 0.0, SVector(0, 0, 0, 0, 3), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(1.41, 0.01, SVector(0, 1, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(1.26, 0.01, SVector(0, 1, 0, 0, -1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 11-20
		PoissonSeries(0.63, 0.0, SVector(1, 0, 0, 0, -1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.63, 0.0, SVector(1, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.46, 0.0, SVector(0, 1, 2, -2, 3), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.45, 0.0, SVector(0, 1, 2, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.36, 0.0, SVector(0, 0, 4, -4, 4), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.24, 0.12, SVector(0, 0, 1, -1, 1), SVector(0, -8, 12, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.32, 0.0, SVector(0, 0, 2, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.28, 0.0, SVector(0, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.27, 0.0, SVector(1, 0, 2, 0, 3), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.26, 0.0, SVector(1, 0, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 21-30
		PoissonSeries(0.21, 0.0, SVector(0, 0, 2, -2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.19, 0.0, SVector(0, 1, -2, 2, -3), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.18, 0.0, SVector(0, 1, -2, 2, -1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.1, -0.05, SVector(0, 0, 0, 0, 0), SVector(0, 8, -13, 0, 0, 0, 0, 0, -1)),
		PoissonSeries(-0.15, 0.0, SVector(0, 0, 0, 2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.14, 0.0, SVector(2, 0, -2, 0, -1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.14, 0.0, SVector(0, 1, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.14, 0.0, SVector(1, 0, 0, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.14, 0.0, SVector(1, 0, 0, -2, -1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.13, 0.0, SVector(0, 0, 4, -2, 4), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 31-33
		PoissonSeries(0.11, 0.0, SVector(0, 0, 2, -2, 4), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.11, 0.0, SVector(1, 0, -2, 0, -3), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.11, 0.0, SVector(1, 0, -2, 0, -1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
	],

	# j = 1
	[
		# 34-36
		PoissonSeries(-0.07, 3.57, SVector(0, 0, 0, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(1.71, -0.03, SVector(0, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 0.48, SVector(0, 0, 2, -2, 3), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
	],

	# j = 2
	[
		# 37-46
		PoissonSeries(743.53, -0.17, SVector(0, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(56.91, 0.06, SVector(0, 0, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(9.84, -0.01, SVector(0, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-8.85, 0.01, SVector(0, 0, 0, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-6.38, -0.05, SVector(0, 1, 0, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-3.07, 0.0, SVector(1, 0, 0, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(2.23, 0.0, SVector(0, 1, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(1.67, 0.0, SVector(0, 0, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(1.3, 0.0, SVector(1, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.93, 0.0, SVector(0, 1, -2, 2, -2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 47-56
		PoissonSeries(0.68, 0.0, SVector(1, 0, 0, -2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.55, 0.0, SVector(0, 0, 2, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.53, 0.0, SVector(1, 0, -2, 0, -2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.27, 0.0, SVector(0, 0, 0, 2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.27, 0.0, SVector(1, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.26, 0.0, SVector(1, 0, -2, -2, -2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.25, 0.0, SVector(1, 0, 0, 0, -1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.22, 0.0, SVector(1, 0, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.21, 0.0, SVector(2, 0, 0, -2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.2, 0.0, SVector(2, 0, -2, 0, -1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 57-61
		PoissonSeries(0.17, 0.0, SVector(0, 0, 2, 2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.13, 0.0, SVector(2, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.13, 0.0, SVector(2, 0, 0, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.12, 0.0, SVector(1, 0, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.11, 0.0, SVector(0, 0, 2, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
	],

	# j = 3
	[
		# 62-65
		PoissonSeries(0.3, -23.51, SVector(0, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.03, -1.39, SVector(0, 0, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.01, -0.24, SVector(0, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 0.22, SVector(0, 0, 0, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
	],

	# j = 4
	[
		# 66-66
		PoissonSeries(-0.26, -0.01, SVector(0, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
	],

]
