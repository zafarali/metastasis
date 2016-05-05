import numpy as np
import matplotlib.pyplot as plt

class Ellipse():
    def __init__(self, a, x_0=0, y_0=0, ecc=0, rotate_by=0, b=None):
        """
            a, minor axis
            x_0 = 0, y_0 = 0, center of ellipse 
            ecc = 0, eccentricity of the ellipse
            rotate_by = 0 , the rotation of the ellipse
            b = 0, major axis if you do not want to define b through eccentricity
        """
        self.POINTS = np.linspace(0,2*np.pi, 150)
        
        self.b = b if b else np.sqrt( (1-ecc**2) * a**2 )
        self.a, self.x_0, self.y_0, self.ecc, self.rotate_by = np.float_(a), x_0, y_0, ecc, rotate_by
        
        x_s = np.cos(self.POINTS)
        y_s = np.sin(self.POINTS)
        
        self.cos_rotate_by = np.cos(rotate_by)
        self.sin_rotate_by = np.sin(rotate_by)
        self.x = self.x_0 + self.a * x_s * self.cos_rotate_by - self.b * y_s * self.sin_rotate_by
        self.y = self.y_0 + self.a * x_s * self.sin_rotate_by + self.b * y_s * self.cos_rotate_by
        
    def plot(self):
        """
            Plot this ellipse
        """
        plt.figure()
        plt.plot(self.x, self.y)
        plt.gca().set_aspect(1)
        plt.show()
    
    def is_inside(self, x_test, y_test):
        """
            checks if a test point is inside the defined ellipse.
        """    
        term_1 = ((self.cos_rotate_by * (x_test - self.x_0) + self.sin_rotate_by * (y_test - self.y_0)) / self.a )**2
        term_2 = ((self.sin_rotate_by * (x_test - self.x_0) + self.cos_rotate_by * (y_test - self.y_0)) / self.b )**2
        
        return ( term_1 + term_2 ) <= 1
        
    def is_inside_vectorized(self):
        raise NotImplementedError()
    
    