import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.path import Path
from PostProcess import PostProcess, H, number_of_segregating_sites, proportion_pairwise_differences, distance_between

# from mpl_toolkits.mplot3d import Axes3D
# from collections import Counter
# from Cell import Cell
# from Util import Ellipse
# import itertools
import random
import numpy as np
import csv
import math
import time
import sys
# import pickle

# generates a coordinate somewhere in the tumor
def coordinate_generator():
	while True:
		yield np.random.random()*500 + 250

# generates an angle somewhere in the tumor
def angle_generator():
	while True:
		yield np.random.random() * 2* np.pi


class Statistics(object):
	def __init__(self, pp):
		"""
			Allows simple statistics collection
			@params:
				pp : a PostProcess module
		"""

		self.pp = pp

	# def random_sampling(self, radii, angle, x, y N):
	# 	"""
	# 		Samples the cells in the tumor randomly. The size of the patch
	# 		is radii
	# 		@params:
	# 			radii : array of two radi for major and minor axis
	# 			N: number of cells in the sample
	# 	"""
	#
	# 	final_sample = []
	# 	proportion_cancer = None
	# 	x_gen = coordinate_generator()
	# 	y_gen = coordinate_generator()
	# 	angle_gen = angle_generator()
	#
	#
	#
	# 		# sample from the patch
	# 		sample = self.pp.cells_in_ellipse_at(x, y, \
	# 			0 ,radii, type_restrictions = [ 1, 2, 3 ], rotate_by = angle )
	#
	#
	#
	# 	return final_sample, N, proportion_cancer


	def simple_stats(self, radii, iterations):
		"""
			Calculates statistics for a certains, e, N, t
			@params:
				e / 0 < float <= 1
					eccentricity of the sample
				N / int
					size of the sample
				t / 0 < float <=1
					the thresholding value
			@returns:
				tuple of:
				(ecc, N, T, S, SH, Epi, D, proportion_cancer)
				and corresponding standard deviations

				returns None if no sample could be generated
		"""

		to_be_returned = []

		x_gen = coordinate_generator()
		y_gen = coordinate_generator()
		angle_gen = angle_generator()

		tumor_COM = self.pp.tumor_COM()

		# number of samples to generate
		for i in xrange(iterations):
			x = x_gen.next()
			y = y_gen.next()
			angle = angle_gen.next()

			sample = self.pp.cells_in_ellipse_at(x, y, \
				0 , radii, type_restrictions = [ 1, 2, 3 ], rotate_by = angle )

			# distance between sample center and tumor center
			sample_COM = self.pp.COM(sample)

			# calculate statistics
			selected_cells = [ cell.id for cell in sample ]
			fa = self.pp.frequency_analyze( selected_cells )
			stats = self.get_stats(fa)
			N = len(sample)
			N_normal = int( np.sum( [ 1 if cell.type == 1 else 0 for cell in sample] ) )

			sample_statistics = {
				'distance': distance_between(tumor_COM, sample_COM),
				'N': N,
				'N_normal': N_normal,
				'N_cancer': N - N_normal,
				'proportion_cancer':(N-N_normal)/float(N),
				'a': radii[0],
				'b': radii[1],
				'x': x,
				'y': y,
				'angle': angle,
				'S':stats['S'],
				'D':stats['D'],
				'SH':stats['SH'],
				'Epi':stats['Epi']
			}
			# save to be returned
			to_be_returned.append( sample_statistics )
		# end for

		return to_be_returned

	def get_stats(self, fa):

		Epi = proportion_pairwise_differences(fa)
		S = number_of_segregating_sites(fa)
		# print 'fa[0]', fa[0]
		# print 'fa[1]',fa[1]

		SH = S/float(H(fa[1]-1))
		D = Epi-SH

		return {'S':S, 'SH':SH, 'Epi':Epi, 'D':D}
