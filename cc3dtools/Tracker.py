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
    def __init__ ( self , file_name = None , template = None , preprocessor = None):
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

        self.internal_stash = [] # stores the rows for stashing

        if template is not None:
            with open( template , 'r' ) as f:
                reader = csv.reader( f )

                for row in reader:
                    self.internal_stash.append( row )

        self.file_name = file_name

        self.preprocessor = preprocessor
        

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
        if self.preprocessor:
            to_stash = self.preprocessor.preprocess(to_stash)

        # preprocess.preprocess returns FALSE if we should skip this stash object!
        if to_stash:
            self.internal_stash.append( to_stash )
        else:
            pass

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


# Generic TrackerPreprocessor
class TrackerPreprocessor(object):
    def __init__(self, fn):
        """
            A preprocessor class that has a function (TrackerPreprocessor.preprocess()) to be called before each Tracker.stash() event
            @params:
                fn / function [mandatory]
                    A function that takes in a list of data, does preprocessing on it.
                    If the data is valid, it returns the processed data
                    If the data in invalid, it returns false so that Tracker will skip it
        """
        self.fn = fn

    def preprocess(self, data):
        """
            calls the preprocessor
            @params:
                data / list
                    data that must be preproessed
        """
        return self.fn(data)


# this function creates a division preprocessor function according to our variables needed
# this might be an overly complicated way of doing it but maybe not?
def generate_divison_preprocessor( start_time = 0 , remove_roots = False ):
    """
        A generator function that returns a function based on a time shift in MCS
        @params:
            start_time / int / 0
                a timeshift to shift data[0] by
            remove_roots / bool / False
                will skip new roots i.e data[0] == 'R'
        @return:
            generated_function
                a function that takes in data and preprocesses it.
    """
    def generated_function(data):
        this = {
            'start_time' : start_time,
            'remove_roots' : remove_roots
        }

        if data[0] == 'R' and remove_roots:
            return False
        else:
            data[0] = data[0] + start_time
            return data
    return generated_function
#end generator

def generate_logger_preprocessor( start_time = 0 ):
    """
        A generator function that returns a function based on a time shift in MCS
        @params:
            start_time / int / 0
        @returns:
            a generated function
    """
    def generated_function( data ):
        this = {
            'start_time' : start_time
        }

        data[0] = data[0] + start_time
        return data
    return generated_function