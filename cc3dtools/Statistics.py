import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.path import Path
import PostProcess

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
		radii = [a, b]

		# sample N_samples times
		threshold_met = False
		max_out = 0

		N_true = None # the true number of cells in this sample

		final_sample = []
		proportion_cancer = None

		while not threshold_met and max_out < 10:
			# generate coordinates
			x = generate_coordinate()
			y = generate_coordinate()
			angle = generate_angle()

			# sample from the patch
			normal_sample = self.pp.cells_in_ellipse_at(x, y, \
				0 ,radii, type_restrictions = [ 1 ], rotate_by = angle )
			cancer_sample = self.pp.cells_in_ellipse_at(x, y, 0, radii, \
				type_restrictions = [ 2 , 3 ], rotate_by = angle )
			
			N_true = len(normal_sample) + len(cancer_sample)

			if N_true < N:
				# we have fewer samples than required, skip this round
				print 'sample size ',N_true,' is too small compared to required ',N
				max_out +=1
				continue

			# the size of the cancer cells we must pick to get this
			# threshold:
			N_cancer_subsample = int(N_true*t)

			# we pick the minimum of the available and our required size
			N_subsample_cancer = min( len(cancer_sample), N_cancer_subsample )
			N_subsample_normal = min( len(normal_sample), N_true - subsample_cancer )

			normal_subsample = random.sample(normal_sample, N_subsample_normal)
			cancer_subsample = random.sample(cancer_sample, N_subsample_cancer)

			proportion_cancer =  len(cancer_subsample) / float( len(cancer_subsample) + len(normal_subsample) )

			if proportion_cancer >= t:
				threshold_met = True
				break	

			final_sample = normal_subsample + cancer_subsample
		#end while

		return final_sample, N_true, proportion_cancer


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
			to_be_averaged.append( ( N_true, t ) + self.get_stats( fa ) + ( proportion_cancer, ) )
		# end repeat samplers
		if len(to_be_averaged) == 0:
			return None
		avgd = tuple( np.mean( np.array(to_be_averaged), axis=0) )
		sds = tuple( np.sd( np.array(to_be_averaged), axis=0) )
		return (e,) + avgd + sds 


	def get_stats(self, fa):

		Epi = PostProcess.proportion_pairwise_differences(fa)
		S = PostProcess.number_of_segregating_sites(fa)
		SH = S/PostProcess.H(fa[1]-1)
		D = Epi-SH

		return (S, SH, Epi, D)


