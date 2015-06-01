## AUTHOR: ZAFARALI AHMED (@ZAFARALI)
## GITHUB.COM/ZAFARALI
## DIVISION TRACKER TO TRACK MITOTIC EVENTS
## ALSO ALLOWS US TO CONVERT BETWEEN MITOTIC EVENTS AND NEWICK FILEFORMAT


class Tracker:
    def __init__ ( self, fileName = '../divisionTrackerOutput.csv' ):
        self.fileName = fileName
        self.stash = []
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
