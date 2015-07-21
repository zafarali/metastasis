from collections import namedtuple

Cell_ = namedtuple('Cell', 'id x y z type initial')
Cell_.__new__.__defaults__ = ( 0, 0, )

class Cell(Cell_):
	"""
		Generic holder for cell details.
		@params:
			id - unique cell identifier
			x, y, z - position of cell in space
			type(=0) - cel type 
			initial(=0) - non-zero if the cell is an initial cell
	"""

