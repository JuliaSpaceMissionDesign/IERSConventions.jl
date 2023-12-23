
const COEFFS_EEQ2000 = [

	# j = 0
	[
		# 1-10
		PoissonSeries(2640.96, -0.39, SVector(0, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(63.52, -0.02, SVector(0, 0, 0, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(11.75, 0.01, SVector(0, 0, 2, -2, 3), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(11.21, 0.01, SVector(0, 0, 2, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-4.55, 0.0, SVector(0, 0, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(2.02, 0.0, SVector(0, 0, 2, 0, 3), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(1.98, 0.0, SVector(0, 0, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1.72, 0.0, SVector(0, 0, 0, 0, 3), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1.41, -0.01, SVector(0, 1, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1.26, -0.01, SVector(0, 1, 0, 0, -1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 11-20
		PoissonSeries(-0.63, 0.0, SVector(1, 0, 0, 0, -1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.63, 0.0, SVector(1, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.46, 0.0, SVector(0, 1, 2, -2, 3), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.45, 0.0, SVector(0, 1, 2, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.36, 0.0, SVector(0, 0, 4, -4, 4), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.24, -0.12, SVector(0, 0, 1, -1, 1), SVector(0, -8, 12, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.32, 0.0, SVector(0, 0, 2, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.28, 0.0, SVector(0, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.27, 0.0, SVector(1, 0, 2, 0, 3), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.26, 0.0, SVector(1, 0, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 21-30
		PoissonSeries(-0.21, 0.0, SVector(0, 0, 2, -2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.19, 0.0, SVector(0, 1, -2, 2, -3), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.18, 0.0, SVector(0, 1, -2, 2, -1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.1, 0.05, SVector(0, 0, 0, 0, 0), SVector(0, 8, -13, 0, 0, 0, 0, 0, -1)),
		PoissonSeries(0.15, 0.0, SVector(0, 0, 0, 2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.14, 0.0, SVector(2, 0, -2, 0, -1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.14, 0.0, SVector(1, 0, 0, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.14, 0.0, SVector(0, 1, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.14, 0.0, SVector(1, 0, 0, -2, -1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.13, 0.0, SVector(0, 0, 4, -2, 4), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 31-33
		PoissonSeries(-0.11, 0.0, SVector(0, 0, 2, -2, 4), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.11, 0.0, SVector(1, 0, -2, 0, -3), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.11, 0.0, SVector(1, 0, -2, 0, -1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
	],

	# j = 1
	[
		# 34-34
		PoissonSeries(-0.87, 0.0, SVector(0, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
	],

]
