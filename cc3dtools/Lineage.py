#### CODE SCAFFOLDING BY SIMON GRAVEL
#### MODIFICATIONS BY ZAFARALI AHMED
from Individual import Individual
import matplotlib.pyplot as plt

class Lineage:

    def __init__ ( self , sub1 = None , sub2 = None , time = 0 , isRoot = False ):
        """
            The lineage class is used to describe ancestral lineages. A lineage 
            either starts from a sample at time 0 (a leaf) or it is the result of the 
            coalescence of two lineages. In either case, it traces back in time 
            until it coalesces.
        
            attributes: 
                sub1
                    the first descendant lineage (default: None if the lineage starts at present)
                sub2
                    the second descendant lineage (default: None if the lineage starts at present)
                    If the lineage starts at present, both sub1 and sub2 are None. 
                If it is formed by coalescence, both are lineages. 
                T0
                    the time at which the lineage was first encountered (the most recent time)
                T1
                    the time at which the lineage coalesces (typically not known when the lineage is first created)
                nsamples
                    the number of samples that trace their ancestry to this lineage
                    length=T1-T0, the total times spent in the lineage. 
        """

        ## if we have only a single individual, this consists of the "root"
        self.isRoot = isRoot

        if isRoot:
            self.sub1 = sub1
            self.sub2 = None
            self.time = 0
            self.nsamples = 1
            return

        ## assertions for making sure we have useable information
        assert sub1 is not None, 'Sub1 cannot be none'
        assert sub2 is not None, 'sub2 cannot be none'
        assert isinstance( sub1 , ( Lineage , Individual ) ) , 'sub1 must either be a lineage or an individual'
        assert isinstance( sub2 , ( Lineage , Individual ) ) , 'sub2 must either be a lineage or an individual'
        assert time >= 0, 'time cannot be negative'

        self.sub1 = sub1 #the first descendant lineage
        self.sub2 = sub2 #the second descendant lineage 
        self.time = time #stores the time that the lineage came about

        if isinstance( sub1 , Individual ) and isinstance( sub2 , Individual ) :
            self.nsamples = 2 #since this is a leaf, there are only 2 individuals at this position
        elif isinstance( sub1 , Individual ) and isinstance( sub2 , Lineage ) :
            self.nsamples = sub2.nsamples + 1
        elif isinstance( sub1 , Lineage ) and isinstance( sub , Individual ):
            self.nsamples = sub1.nsamples + 1
        else:
            self.nsamples = sub1.nsamples + sub2.nsamples
        #     #then the lineage is a leaf, we store individuals as leaves or none
        #     assert isinstance( sub2, Individual ), 'Both sub1 and sub2 should be individuals'
            
        #     self.time = sub1

        #     

        # elif isinstance( sub1, Lineage ): 
        #     #we haven't defined T1 yet, but we assume that coalescing lineages will 
        #     #have a coalescing time:

        #     assert isinstance( sub2, Lineage ), 'sub2 must also be a lineage!'

        #     assert sub1.time == sub2.T1,"coalescing lineages should coalesce at the same time!"
            
        #     self.time = sub2.time #the founding time of the coalesced lineage is the coalescing time of 
        #     #the desending lineages.
        #     #Now compute the number of samples descending from the current lineage. 
        #     #Recursion makes this very simple! We can just add up the number of samples in each
        #     #sublineages
        #     self.nsamples = sub1.nsamples + sub2.nsamples #s/=.*/=.../
        # else:
        #     raise Exception(' sub1 is not a supported type, it must be either a lineage or an individual ')

    # def coalesce( self , lineage2 , time ):
    #     """
    #         coalesce the current lineage with a second lineage lineage2 at time t
    #         Updates the coalescence time in both lineages, and returns the coalesced lineage
    #     """
    #     self.T1 = time
    #     lineage2.T1 = time
    #     return Lineage ( self , lineage2 )

    def next_time( self ):
        return self.sub1.time

    def draw ( self, width ):
        final_points = []
        plt.figure()
        self.plot( 0 , self.num_descendants() , width , final_points )
        plt.xticks( [] )
        plt.xlabel( 'Cells' )
        plt.yticks([])
        plt.ylabel( 'Division' )
        plt.title('Tree')
        plt.show()
        return final_points

    def plot( self , center , generation , width , final_points ):
        """
            Plots the lineage and all its descending lineages. To avoid overlap when plotting
            multiple lineages, we specify the width and center of the lineage along the x-axis.
            Time is plotted along the y axis, and uses the lineage T0 and T1
        """

        # plot the vertical line 
        plt.plot( [ center , center ] , [ generation , generation - 1 ] , 'b' )

        n1 = self.sub1.num_descendants() or 1
        n2 = self.sub2.num_descendants() or 1
        w1 = n1 * 1. / ( n1 + n2 ) * width
        w2 = n2 * 1. / ( n1 + n2 ) * width
        mid1 = center - width / 2. + w1 / 2.
        mid2 = center + width / 2. - w2 / 2.

        #plot horizontal line
        plt.plot( [ mid1 , mid2 ] , [ generation - 1 , generation - 1 ] , 'b' )

        self.sub1.plot( mid1 , generation - 1 , w1 , final_points )
        self.sub2.plot( mid2 , generation - 1 , w2 , final_points )





        # check if this is a terminal lineage by seeing if both subchildren are leaves.
        # if isinstance( self.sub1 , Lineage ) and isinstance( self.sub2 , Lineage ) :
        #     #assign width proportional to the number of lineages in each sub-lineage
        #     n1=self.sub1.num_descendants()
        #     n2=self.sub2.num_descendants()
        #     w1=n1*1./(n1+n2)*width
        #     w2=n2*1./(n1+n2)*width
        #     mid1=center-width/2.+w1/2. #Find the center of each window
        #     mid2=center+width/2.-w2/2.
        #     plt.plot([mid1,mid2],[self.time,self.time],'b') #plot horizontal connector
        #     self.sub1.plot(w1,mid1) #plot descending lineages
        #     self.sub2.plot(w2,mid2)
        #     print 'terminated'

        # if isinstance( self.sub1 , Lineage ):
        #     n1 = self.sub1.num_descendants()
        #     n2 = 1
        #     w1 = n1*1./(n1+n2)*width
        #     w2 = n2*1./(n1+n2)*width
        #     mid1 = center-width/2.+w1/2. #Find the center of each window
        #     mid2 = center+width/2.-w2/2.

        #     plt.plot([mid1,mid2],[self.time,0],'b') #plot horizontal connector
        #     self.sub1.plot(w1,mid1) #plot descending lineages
        #     # add text to the position of sub2.

        # if isinstance ( self.sub2 , Lineage ):
        #     n2 = self.sub2.num_descendants()
        #     n1 = 1
        #     w1 = n1*1./(n1+n2)*width
        #     w2 = n2*1./(n1+n2)*width
        #     mid1 = center-width/2.+w1/2. #Find the center of each window
        #     mid2 = center+width/2.-w2/2.
        #     plt.plot([mid1,mid2],[self.time,0],'b') #plot horizontal connector
        #     self.sub2.plot(w2,mid2) #plot descending lineages
        #     # add text to the position of sub1



    # def get_length(self):
    #     """
    #         returns a tuple whose first element is the time spent in the current lineage, 
    #         and the second element is total length spent in the present lineage plus all descending lineages
    #     """
    #     self.length=self.T1-self.T0
    #     if self.sub1 is None: #then the time spent in the lineage and the total time are the same
    #         return (self.length,self.length)
    #     else:
    #         #to get the total time, use the recursion property again!
    #         #First compute the total length within the sub1 and sub2 lineages
    #         lengths1=self.sub1.get_length()[1]#s/=.*/.../
    #         lengths2=self.sub2.get_length()[1]#s/=.*/.../
            
    #         return (self.length,self.length+lengths1+lengths2) #s/\,self.*/,...)/
    
    def __repr__( self ):
        return '('+str(self.sub1)+', '+str(self.sub2)+')'

    def find ( self , key ):
        """
            Returns a tuple of (lineage, individual, LR) where the lineage contains 
            the individual we are searching for and LR(int) contains if its the first or second
            child or False otherwise
        """

        if isinstance( self.sub1 , Individual ):

            if self.sub1.id == key or self.sub1.name == key:

                return ( self , self.sub1 , 1 )


        if isinstance( self.sub2 , Individual ):

            if self.sub2.id == key or self.sub2.name == key:

                return ( self , self.sub2, 2 )

        ## we didn't find the individual, thus we continue searching downward
        sub1result = False
        if isinstance( self.sub1 , Lineage ):

            sub1result = self.sub1.find( key )

        sub2result = False
        if isinstance( self.sub2 , Lineage ):

            sub2result = self.sub2.find( key )

        ## we didn't find the individual even though we checked the subtrees, thus we return False

        if not sub1result and not sub2result:
            return False
        else:
            return sub1result or sub2result

    def divide ( self , parent , child1 , child2 = None , time = 0 ):
        """
            replaces an individual with children
        """


        search = self.find( parent )
        if search is False:
            return False
        container, replaced, lr = search  #lineage that contains the parent
        
        if child2 is None:
            child2 = replaced
        
        if self.isRoot:
            self.sub1 = child1
            self.sub2 = child2
            self.isRoot = False
            self.nsamples = 2
            return

        if lr == 1:
            container.sub1 = Lineage( sub1 = child1 , sub2 = child2 , time = time )
        elif lr == 2: 
            container.sub2 = Lineage( sub1 = child1 , sub2 = child2 , time = time )

        container.nsamples+=2

    def to_newick ( self ): 
        """
            convert this lineage to the newick standard for representing
            phylogeny for easy use later
        """
        return str(self)+';'

    def save_verbose ( self , format = ('t', 'p', 'c1', 'c2') ):
        pass

    def num_descendants ( self ):
        """
            Returns the number of descendants of the lineage
        """
        if isinstance( self.sub1 , Individual ) and isinstance( self.sub2 , Individual ) :
            return 2 #since this is a leaf, there are only 2 individuals at this position
        elif isinstance( self.sub1 , Individual ) and isinstance( self.sub2 , Lineage ) :
            return self.sub2.num_descendants() + 1
        elif isinstance( self.sub1 , Lineage ) and isinstance( self.sub2 , Individual ):
            return self.sub1.num_descendants() + 1
        else:
            return self.sub1.num_descendants() + self.sub2.num_descendants()

    # def descendants ( self ):
    #     if isinstance( self.sub1 , Individual ) and isinstance( self.sub2 , Individual ) :
    #         return [ self.sub1 , self.sub2 ] #since this is a leaf, there are only 2 individuals at this position
    #     elif isinstance( self.sub1 , Individual ) and isinstance( self.sub2 , Lineage ) :
    #         return [ self.sub1 ].extend( self.sub2.descendants() )
    #     elif isinstance( self.sub1 , Lineage ) and isinstance( self.sub2 , Individual ):
    #         return [ self.sub2 ].extend( self.sub1.descendants() )
    #     else:
    #         return self.sub1.descendants().extend( self.sub2.descendants() )
       

    @staticmethod
    def load_file( fileName = None , firstParentId = 1 , format = ( 't' , 'p' , 'c1' , 'c2' ) ):
        """
            This method loads a csv file with the given format and converts
            it into a lineage that can be manipulated or visualized.
            firstParent is the name of the root
        """

        import csv

        out = Lineage( sub1 = Individual( firstParentId , name = str( firstParentId ) ) , time = 0, isRoot = True )

        with open( fileName , 'r' ) as f :
            reader = csv.reader( f )
            for row in reader:
                info = dict( zip( format , row ) )
                out.divide( parent = info['p'] , \
                    child1 = Individual( int( info['c1'] ) , name = info['c1'] ) , \
                    child2 = Individual( int( info['c2'] ) , name = info['c2'] ) ,\
                    time = int ( info['t'] ) )
        return out


