
const COEFFS_CIP1980_X = [

	# j = 0
	[
		# 1-10
		PoissonSeries(-171996.0, 0.0, SVector(0, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-13187.0, 0.0, SVector(0, 0, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-2274.0, 0.0, SVector(0, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(2053.2, 0.0, SVector(0, 0, 0, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1426.0, 0.0, SVector(0, -1, 0, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(712.0, 0.0, SVector(1, 0, 0, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-517.0, 0.0, SVector(0, 1, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-386.0, 0.0, SVector(0, 0, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-301.0, 0.0, SVector(1, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(217.0, 0.0, SVector(0, -1, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 11-20
		PoissonSeries(158.0, 0.0, SVector(-1, 0, 0, 2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(129.0, 0.0, SVector(0, 0, 2, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(123.0, 0.0, SVector(-1, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(63.0, 0.0, SVector(1, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(63.0, 0.0, SVector(0, 0, 0, 2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-59.0, 0.0, SVector(-1, 0, 2, 2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-58.0, 0.0, SVector(-1, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-51.0, 0.0, SVector(1, 0, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-48.0, 0.0, SVector(-2, 0, 0, 2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(46.0, 0.0, SVector(-2, 0, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 21-30
		PoissonSeries(-38.0, 0.0, SVector(0, 0, 2, 2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-31.0, 0.0, SVector(2, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(29.0, 0.0, SVector(2, 0, 0, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(29.0, 0.0, SVector(1, 0, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(26.0, 0.0, SVector(0, 0, 2, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-22.0, 0.0, SVector(0, 0, 2, -2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(21.0, 0.0, SVector(-1, 0, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(17.0, 0.0, SVector(0, 2, 0, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-16.0, 0.0, SVector(0, 2, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(16.0, 0.0, SVector(-1, 0, 0, 2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 31-40
		PoissonSeries(-15.0, 0.0, SVector(0, 1, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-13.0, 0.0, SVector(1, 0, 0, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-12.0, 0.0, SVector(0, -1, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(11.0, 0.0, SVector(2, 0, -2, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-10.0, 0.0, SVector(-1, 0, 2, 2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-8.0, 0.0, SVector(1, 0, 2, 2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-7.0, 0.0, SVector(0, -1, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-7.0, 0.0, SVector(0, 0, 2, 2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-7.0, 0.0, SVector(1, 1, 0, -2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(7.0, 0.0, SVector(0, 1, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 41-50
		PoissonSeries(-6.0, 0.0, SVector(-2, 0, 0, 2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-6.0, 0.0, SVector(0, 0, 0, 2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(6.0, 0.0, SVector(2, 0, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(6.0, 0.0, SVector(1, 0, 0, 2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(6.0, 0.0, SVector(1, 0, 2, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-5.0, 0.0, SVector(0, 0, 0, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-5.0, 0.0, SVector(0, -1, 2, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-5.0, 0.0, SVector(2, 0, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(5.0, 0.0, SVector(1, -1, 0, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-4.0, 0.0, SVector(1, 0, 0, -1, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 51-60
		PoissonSeries(-4.0, 0.0, SVector(0, 0, 0, 1, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-4.0, 0.0, SVector(0, 1, 0, -2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(4.0, 0.0, SVector(1, 0, -2, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(4.0, 0.0, SVector(2, 0, 0, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(4.0, 0.0, SVector(0, 1, 2, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-3.0, 0.0, SVector(1, 1, 0, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-3.0, 0.0, SVector(1, -1, 0, -1, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-3.0, 0.0, SVector(-1, -1, 2, 2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-3.0, 0.0, SVector(0, -1, 2, 2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-3.0, 0.0, SVector(1, -1, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 61-70
		PoissonSeries(-3.0, 0.0, SVector(3, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-3.0, 0.0, SVector(-2, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(3.0, 0.0, SVector(1, 0, 2, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-2.0, 0.0, SVector(-1, 0, 2, 4, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-2.0, 0.0, SVector(1, 0, 0, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-2.0, 0.0, SVector(-1, 0, 2, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-2.0, 0.0, SVector(0, -2, 2, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-2.0, 0.0, SVector(-2, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(2.0, 0.0, SVector(2, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(2.0, 0.0, SVector(3, 0, 0, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 71-80
		PoissonSeries(2.0, 0.0, SVector(1, 1, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(2.0, 0.0, SVector(0, 0, 2, 1, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1.0, 0.0, SVector(1, 0, 0, 2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1.0, 0.0, SVector(1, 0, 2, 2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1.0, 0.0, SVector(1, 1, 0, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1.0, 0.0, SVector(0, 1, 0, 2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1.0, 0.0, SVector(0, 1, 2, -2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1.0, 0.0, SVector(0, 1, -2, 2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1.0, 0.0, SVector(1, 0, -2, 2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1.0, 0.0, SVector(1, 0, -2, -2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 81-90
		PoissonSeries(-1.0, 0.0, SVector(1, 0, 2, -2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1.0, 0.0, SVector(1, 0, 0, -4, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1.0, 0.0, SVector(2, 0, 0, -4, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1.0, 0.0, SVector(0, 0, 2, 4, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1.0, 0.0, SVector(0, 0, 2, -1, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1.0, 0.0, SVector(-2, 0, 2, 4, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1.0, 0.0, SVector(2, 0, 2, 2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1.0, 0.0, SVector(0, -1, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1.0, 0.0, SVector(0, 0, -2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(1.0, 0.0, SVector(0, 0, 4, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 91-100
		PoissonSeries(1.0, 0.0, SVector(0, 1, 0, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(1.0, 0.0, SVector(1, 1, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(1.0, 0.0, SVector(3, 0, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(1.0, 0.0, SVector(-2, 0, 2, 2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(1.0, 0.0, SVector(-1, 0, 0, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(1.0, 0.0, SVector(0, 0, -2, 2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(1.0, 0.0, SVector(0, 1, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(1.0, 0.0, SVector(-1, 0, 4, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(1.0, 0.0, SVector(2, 1, 0, -2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(1.0, 0.0, SVector(2, 0, 0, 2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 101-107
		PoissonSeries(1.0, 0.0, SVector(2, 0, 2, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(1.0, 0.0, SVector(2, 0, -2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(1.0, 0.0, SVector(1, -1, 0, -2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(1.0, 0.0, SVector(-1, 0, 0, 1, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(1.0, 0.0, SVector(-1, -1, 0, 2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(1.0, 0.0, SVector(0, 1, 0, 1, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1.2, 0.0, SVector(0, 0, 2, -2, 3), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
	],

	# j = 1
	[
		# 108-117
		PoissonSeries(-84.2, 5173.2, SVector(0, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(5.3, 322.2, SVector(0, 0, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(1.0, 54.8, SVector(0, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1.0, -50.5, SVector(0, 0, 0, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(4.3, 3.0, SVector(0, -1, 0, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.1, 0.0, SVector(1, 0, 0, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(1.5, 12.6, SVector(0, 1, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.4, 11.3, SVector(0, 0, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 7.3, SVector(1, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.5, -5.3, SVector(0, -1, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 118-126
		PoissonSeries(0.1, -4.0, SVector(0, 0, 2, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -3.0, SVector(-1, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.1, -1.8, SVector(1, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 1.5, SVector(-1, 0, 2, 2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.1, 1.8, SVector(-1, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 1.5, SVector(1, 0, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -1.3, SVector(-2, 0, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.1, 0.0, SVector(0, 2, 0, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.1, 0.0, SVector(0, 2, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
	],

]

const COEFFS_CIP1980_Y = [

	# j = 0
	[
		# 1-10
		PoissonSeries(0.0, 92025.0, SVector(0, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 5736.0, SVector(0, 0, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 977.0, SVector(0, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -893.7, SVector(0, 0, 0, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 54.0, SVector(0, -1, 0, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -7.0, SVector(1, 0, 0, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 224.0, SVector(0, 1, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 200.0, SVector(0, 0, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 129.0, SVector(1, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -95.0, SVector(0, -1, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 11-20
		PoissonSeries(0.0, -1.0, SVector(-1, 0, 0, 2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -70.0, SVector(0, 0, 2, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -53.0, SVector(-1, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -33.0, SVector(1, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -2.0, SVector(0, 0, 0, 2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 26.0, SVector(-1, 0, 2, 2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 32.0, SVector(-1, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 27.0, SVector(1, 0, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 1.0, SVector(-2, 0, 0, 2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -24.0, SVector(-2, 0, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 21-30
		PoissonSeries(0.0, 16.0, SVector(0, 0, 2, 2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 13.0, SVector(2, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -1.0, SVector(2, 0, 0, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -12.0, SVector(1, 0, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -1.0, SVector(0, 0, 2, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -10.0, SVector(-1, 0, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 7.0, SVector(0, 2, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -8.0, SVector(-1, 0, 0, 2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 9.0, SVector(0, 1, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 7.0, SVector(1, 0, 0, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 31-40
		PoissonSeries(0.0, 6.0, SVector(0, -1, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 5.0, SVector(-1, 0, 2, 2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 3.0, SVector(1, 0, 2, 2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 3.0, SVector(0, -1, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 3.0, SVector(0, 0, 2, 2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -3.0, SVector(0, 1, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 3.0, SVector(-2, 0, 0, 2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 3.0, SVector(0, 0, 0, 2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -3.0, SVector(2, 0, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -3.0, SVector(1, 0, 2, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 41-50
		PoissonSeries(0.0, 3.0, SVector(0, 0, 0, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 3.0, SVector(0, -1, 2, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 3.0, SVector(2, 0, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -2.0, SVector(2, 0, 0, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -2.0, SVector(0, 1, 2, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 1.0, SVector(-1, -1, 2, 2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 1.0, SVector(0, -1, 2, 2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 1.0, SVector(1, -1, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 1.0, SVector(3, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 1.0, SVector(-2, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 51-60
		PoissonSeries(0.0, 1.0, SVector(-1, 0, 2, 4, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 1.0, SVector(1, 0, 0, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 1.0, SVector(-1, 0, 2, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 1.0, SVector(0, -2, 2, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 1.0, SVector(-2, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -1.0, SVector(2, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -1.0, SVector(1, 1, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -1.0, SVector(0, 0, 2, 1, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 1.0, SVector(1, 0, 2, 2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, 1.0, SVector(-2, 0, 2, 4, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 61-64
		PoissonSeries(0.0, -1.0, SVector(1, 1, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -1.0, SVector(-2, 0, 2, 2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -1.0, SVector(-1, 0, 0, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.0, -1.0, SVector(2, 0, 2, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
	],

	# j = 1
	[
		# 65-74
		PoissonSeries(1529.9, 8.9, SVector(0, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(117.3, -3.1, SVector(0, 0, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(20.2, -0.5, SVector(0, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-18.3, 0.5, SVector(0, 0, 0, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(12.7, -0.1, SVector(0, -1, 0, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-6.3, 0.0, SVector(1, 0, 0, 0, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(4.6, -0.6, SVector(0, 1, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(3.4, 0.0, SVector(0, 0, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(2.7, -0.1, SVector(1, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1.9, 0.3, SVector(0, -1, 2, -2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),

		# 75-82
		PoissonSeries(-1.4, 0.0, SVector(-1, 0, 0, 2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1.2, 0.0, SVector(0, 0, 2, -2, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-1.1, 0.0, SVector(-1, 0, 2, 0, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.6, 0.0, SVector(1, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(-0.6, 0.0, SVector(0, 0, 0, 2, 0), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.5, 0.0, SVector(-1, 0, 2, 2, 2), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.5, 0.0, SVector(-1, 0, 0, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
		PoissonSeries(0.5, 0.0, SVector(1, 0, 2, 0, 1), SVector(0, 0, 0, 0, 0, 0, 0, 0, 0)),
	],

]
