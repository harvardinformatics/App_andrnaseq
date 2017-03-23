import sys
import re
import unicodedata
import time
from bioservices import BioDBNet


def getAcc(fn, column):
    '''
    column: the name of the column to get
    '''
    with open(fn, 'r') as f:
        accL = []
        line = f.readline()
        ix = line.split('\t').index(column)
        ix += 1 # header has one less columns
        
        for line in f:
            line = line.split('\t')
            acc = line[ix]
            if not re.match('[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}', acc):
                print '%s is not an accession' % acc 
            else:
                accL.append(line[ix])                
    return accL


class Conversion:
    def __init__(self, ids, idtype, attr):
        self.s = BioDBNet()
        self.input_db = idtype # 'UniProt Accession'; use: s.getInputs()
        self.output_db = attr # 'Gene Symbol', 'UniProt Protein Name'; use: s.getOutputsForInputs('Uniprot Accession')
        self.inputValues = self.chunks(ids, 500) # ids must be a list

        print self.input_db
        print self.output_db

    def chunks(self, l, n):
        '''Yield successive n-sized chunks from l.
            
        code from stackoverflow.
        '''
        for i in range(0, len(l), n):
            yield l[i:i+n]

    def writeString(self, reqstr, filename):
        '''Write output from request to file with 2 columns.'''
        a = re.sub('\n', ',', reqstr).split(',')[1:-1]
        with open(filename, 'a') as of:
            for el in a:
                of.write(el + '\n')
            
    def submit(self): # vals is a generator
        l = []
        for vall in self.inputValues:
            vals = ','.join(vall)
            x = self.s.db2db(self.input_db, self.output_db, vals)
            x = unicodedata.normalize('NFKC', x).encode('ascii','ignore')
            l.append(x)
            time.sleep(10)
        return l

    def write2file(self, strlst, fn):
        for strng in strlst:
            self.writeString(strng, fn)


#s.db2db('UniProt Accession', ['Gene Symbol'], ['P78527, P15924, Q9NU22, P49327'])

if __name__ == '__main__':
    if len(sys.argv) < 3:
        usage = 'Usage: %s <file with UniProt accessions, e.g., x238.df> <e.g, \'Gene Symbol\', \'UnProt Accession\'>' % sys.argv[0]
        print usage
        print 'Example command line: %s' % "python getProtAnnotv2.py ../TempData/x243.txt 'Gene Symbol'"
        sys.exit(0)

    file = sys.argv[1]
    accL = getAcc(file, 'Accession') # python getProtAnnotv2.py TempData/x238.txt (tab delimited)

    if sys.argv[2] == 'Gene Symbol':
        cv = Conversion(accL, 'UniProt Accession', ['Gene Symbol'])
        res = cv.submit()
        cv.write2file(res, 'upAccession2GeneSymbol')
    elif sys.argv[2] == 'UniProt Protein Name':
        cv = Conversion(accL, 'UniProt Accession', ['UniProt Protein Name'])
        res = cv.submit()
        cv.write2file(res, 'upAccession2upProtName')
    else:
        print '%s not implemented, yet' % sys.argv[2]
    



            

