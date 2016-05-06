import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.path import Path
from PostProcess import PostProcess, H, number_of_segregating_sites, proportion_pairwise_differences

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
def generate_coordinate():
	while True:
		yield np.random.random()*500 + 250

# generates an angle somewhere in the tumor
def generate_angle():
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

	def random_sampling(self, eccentricity=1, N=1000, N_samples=10, t=0.1):
		"""
			Samples the tumor randomly according to some
			threshold
			@params:
				N / int / 1000:
					the number of cells we want have totally
				eccentricity / int / 1:
					the shape of the ellipse we want to sample
				N_samples / int / 10:
					the number of samples we want to get to find an average
				t / float / 0.1:
					the threshold of cancer cells we want


		"""

		average_area = 70.0 #average area of a cell
		# from calculations of the radii of the sample
		# (area of patch = average_area * N)
		# (area of an ellipse = pi * a * b )
		# write a in terms of b and eccentricity
		# equate the two areas
		a = lambda N, ecc: np.sqrt( ( average_area * N ) / ( np.pi * np.sqrt( 1 - ecc**2 ) ) )
		b = lambda N, ecc: np.sqrt( ( ( average_area * N ) * np.sqrt( 1 - ecc**2 ) ) / np.pi )
		radii = [a(N, eccentricity), b(N, eccentricity)]

		# sample N_samples times
		threshold_met = False
		max_out = 0

		N_true = None # the true number of cells in this sample

		final_sample = []
		proportion_cancer = None
		x_gen = generate_coordinate()
		y_gen = generate_coordinate()
		angle_gen = generate_coordinate()

		while not threshold_met and max_out < 10:
			# generate coordinates
			x = x_gen.next()
			y = y_gen.next()
			angle = angle_gen.next()

			# sample from the patch
			normal_sample = self.pp.cells_in_ellipse_at(x, y, \
				0 ,radii, type_restrictions = [ 1 ], rotate_by = angle )
			cancer_sample = self.pp.cells_in_ellipse_at(x, y, 0, radii, \
				type_restrictions = [ 2 , 3 ], rotate_by = angle )
			# print 'normal sample:',normal_sample
			# print 'cancer sample:',cancer_sample
			N_true = len(normal_sample) + len(cancer_sample)

			if N_true < N:
				# we have fewer samples than required, skip this round
				print 'sample size ',N_true,' is too small compared to required ',N
				max_out +=1
				continue

			# the size of the cancer cells we must pick to get this
			# threshold:
			N_cancer_subsample = int(N*t)

			# we pick the minimum of the available and our required size
			N_subsample_cancer = min( len(cancer_sample), N_cancer_subsample )
			N_subsample_normal = min( len(normal_sample), N - N_subsample_cancer )

			normal_subsample = random.sample(normal_sample, N_subsample_normal)
			cancer_subsample = random.sample(cancer_sample, N_subsample_cancer)
			# print 'size of cancer_subsample:', len(cancer_subsample)
			# print 'size of normal_subsample:',len(normal_subsample)
			# # print 'type of subsample size', type(len(normal_subsample))
			# print 'total size of the sample:', float( len(cancer_subsample) + len(normal_subsample) )

			proportion_cancer =  len(cancer_subsample) / float( len(cancer_subsample) + len(normal_subsample) )

			if proportion_cancer >= t:
				threshold_met = True
				final_sample = normal_subsample + cancer_subsample
			else:
				max_out+=1
		#end while
		if not N_true or max_out == 9:
			print 'could not get a big enough sample'
		return final_sample, N, proportion_cancer


	def simple_stats(self, e, N, t):
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
		to_be_averaged = []
		for i in xrange(10):
			# repeat for averaging
			sample, N_true, proportion_cancer = self.random_sampling(\
				eccentricity=e, N=N, t=t )
			selected_cells = [ cell.id for cell in sample ]
			if not N_true:
				print 'failed to generate a sample'
				continue
			fa = self.pp.frequency_analyze( selected_cells )
			if fa[1] == 0:
				print 'failed to generate a sample'
				continue
			to_be_averaged.append( ( N_true, t ) + self.get_stats( fa ) + ( proportion_cancer, ) )
		# end repeat samplers
		if len(to_be_averaged) == 0:
			return False
		avgd = tuple( np.mean( np.array(to_be_averaged), axis=0) )
		sds = tuple( np.sd( np.array(to_be_averaged), axis=0) )
		return (e,) + avgd + sds 


	def get_stats(self, fa):

		Epi = proportion_pairwise_differences(fa)
		S = number_of_segregating_sites(fa)
		# print 'fa[0]', fa[0]
		# print 'fa[1]',fa[1]

		SH = S/H(fa[1]-1)
		D = Epi-SH

		return (S, SH, Epi, D)


