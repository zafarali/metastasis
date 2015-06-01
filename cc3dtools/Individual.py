### WRITTEN BY ZAFARALI AHMED
### MAY 2015

class Individual:
    def __init__ ( self , iid,  **kwargs ):
        """
            This class holds details of the individuals at different locations
            in the lineage

        """
        self.id = iid
        self.name = kwargs.get('name' , None )
        ## pre-processing to make it easier for us to get nodes
        # self.children = sorted( children, key = lambda x : ( x is None, x ) )
        self.more = kwargs
    
    # def __cmp__ ( self , other ):
    #     if self.id > other.id:
    #         return 1
    #     elif self.id < other.id:
    #         return -1

    def __repr__ ( self ):
        if self.name is None:
            return str(self.id)
        else: 
            return str(self.name)

    def num_descendants ( self ):
        return 0

    def plot ( self , center , generation , widith , final_points,  ):
        import matplotlib.pyplot as plt
        plt.plot( [ center , center ] , [ generation , generation - 1 ] , 'b' )
        plt.text( center , generation - 1.25 , str(self.name) if self.name is not None else str(self.id) , horizontalalignment='center' , rotation=90 , size='x-small')
        final_points.append( {'name':self.name if self.name is not None else int(self.id), 'location': ( center , generation - 1.25 ) } )
        return