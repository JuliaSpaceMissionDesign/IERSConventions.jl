
const COEFFS_CPNC_XP = SVector(-17251.0, 2004191898.0, -429783.0, -198618.0)

const COEFFS_CPNC_X = [

	# j = 0
	[
		# 1-10
		PoissonSeries(-6.844318e6, 0.0, SVector(0, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(82169.0, 0.0, SVector(0, 0, 0, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(2521.0, 0.0, SVector(0, 0, 0, 2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(5096.0, 0.0, SVector(0, 0, 2, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-523908.0, 0.0, SVector(0, 0, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-15407.0, 0.0, SVector(0, 0, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-90552.0, 0.0, SVector(0, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-8585.0, 0.0, SVector(0, 1, -2, 2, -2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(58707.0, 0.0, SVector(0, 1, 0, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-20558.0, 0.0, SVector(0, 1, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 11-15
		PoissonSeries(-4911.0, 0.0, SVector(1, 0, -2, 0, -2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-6245.0, 0.0, SVector(1, 0, 0, -2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(28288.0, 0.0, SVector(1, 0, 0, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(2512.0, 0.0, SVector(1, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-11992.0, 0.0, SVector(1, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
	],

	# j = 1
	[
		# 16-17
		PoissonSeries(-3310.0, 205833.0, SVector(0, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 12814.0, SVector(0, 0, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
	],

]

const COEFFS_CPNC_YP = SVector(-5530.0, -25896.0, -22407275.0, 1901.0, 1113.0)

const COEFFS_CPNC_Y = [

	# j = 0
	[
		# 1-10
		PoissonSeries(0.0, 9.205236e6, SVector(0, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -89618.0, SVector(0, 0, 0, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -6918.0, SVector(0, 0, 2, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 573033.0, SVector(0, 0, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 20070.0, SVector(0, 0, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 97847.0, SVector(0, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -9593.0, SVector(0, 1, -2, 2, -2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 7387.0, SVector(0, 1, 0, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 22438.0, SVector(0, 1, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 2555.0, SVector(1, 0, -2, -2, -2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 11-15
		PoissonSeries(0.0, -5331.0, SVector(1, 0, -2, 0, -2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 3144.0, SVector(1, 0, 0, 0, -1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -3324.0, SVector(1, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 2636.0, SVector(1, 0, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 12903.0, SVector(1, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
	],

	# j = 1
	[
		# 16-17
		PoissonSeries(153042.0, 0.0, SVector(0, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(11714.0, 0.0, SVector(0, 0, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
	],

]
