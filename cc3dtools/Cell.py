from collections import namedtuple

# class Cell( namedtuple('Cell', 'x y z id type initial') ):
# 	def __new__( _cls , x , y , z , id , type , initial = 0 ):
# 		return super(Call, cls).__new__( _cls, ( x, y , id , type , inital ) ) 

Cell = namedtuple('Cell', 'id x y z type initial')
Cell.__new__.__defaults__ = ( 0, )
