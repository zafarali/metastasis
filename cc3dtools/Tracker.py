## AUTHOR: ZAFARALI AHMED (@ZAFARALI)
## GITHUB.COM/ZAFARALI


import csv

class Tracker:
    def __init__ ( self, fileName = '../divisionTrackerOutput.csv' ):
        self.fileName = fileName
        self.stash = []
        print 'cc3dtools.Tracker is no longer being maintained, use cc3dtools.Tracker2 and its extensions'
    ## this function saves a division event into the file
    ## in the format: time, parent, child1, child2
    def saveDivision ( self , time = 0 , parent = 1 , child1 = 2, child2 = 3 ):
        import csv
        with open ( self.fileName, 'a' ) as f:
            c = csv.writer( f )
            c.writerow( [ time , parent , child1, child2 ] )

    def stashDivision ( self , time = 0 , parent = 1 , child1 = 2, child2 =3 ):
        self.stash.append( [ time , parent , child1 , child2 ] )

    def saveStash( self ):
        import csv
        with open ( self.fileName, 'a' ) as f:
            c = csv.writer( f )
            for row in self.stash:
                c.writerow( row )

    def createLineage( self ):
        pass

    ## this function transfrms the output file
    ## into a newick standard tree format
    ## (A,B,(C,D),E);
    def outputToNewick ( self , newickFileName = None ):
    	if newickFileName is None:
    		newickFileName = self.fileName.split('csv')[0].join('newick')

class Tracker2( object ):
    def __init__ ( self , file_name = None ):
        """ 
            Allows us to store and track various events and properties
            @params:
                file_name / str [ mandatory ] 

            -> ( this class serves as a base class for most of the trackers used in cc3dtools )
        """
        assert file_name is not None , 'file_name must be supplied to use Tracker2'
        try:
            open(file_name, 'a').close()
        except Exception as e:
            print 'FAILED TO LOAD TRACKER2 MODULE'
            print 'ERROR\n'+str(e)
            raise e
        
        self.file_name = file_name
        self.internal_stash = [] # stores the rows for stashing

    def save_now( self , to_save ):
        """
            saves a list immediately to the file
            @params:
                to_save / list
                    list of entries to save immediately to file
        """
        assert type(to_save) is list , 'to_save must be a list'

        with open ( self.fileName, 'a' ) as f:
            c = csv.writer( f )
            c.writerow( to_save )

    def stash( self , to_stash ):
        """
            stores a list in the stash for future saving
            @params:
                to_stash / list
                    list of entries to stash
        """
        assert type( to_stash ) is list , 'to_stash must be a list'
        self.internal_stash.append( to_stash )

    def save_stash( self , flag = 'wb' ):
        """
            saves the current state of the stash to a file
            ( note this is an write function by default )
            @params:
                flag / str / 'w'
                    how to store the stash ('a' = append, 'w' = write )
        """
        assert type( flag ) is str , 'flag must be a string'
        assert flag != 'r' , 'cannot save with a read flag'

        try:
            with open( self.file_name , flag ) as f:
                c = csv.writer( f )
                for row in self.internal_stash:
                    c.writerow( row )

            print 'stash was saved to '+self.file_name
        except Exception as e:
            print 'EXCEPTION OCCURED WHEN TRYING TO SAVE TO ',self.file_name
            print 'The state of the stash is: \n',self.internal_stash
            print '\n Exception was: ',str(e)

