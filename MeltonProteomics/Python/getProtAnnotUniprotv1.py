import sys
import re
import unicodedata
import time
from bioservices import UniProt


def getAcc(column, *files):
    '''
    column: the name of the column to get
    '''
    accL = []
    for file in files:
        with open(file, 'r') as f:
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

    accS = set(accL)
    return list(accS)

def _getAcc(column, flst):
    '''
    column: the name of the column to get
    flst: a list of file names 
    '''
    accL = []
    for file in flst:
        with open(file, 'r') as f:
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

    accS = set(accL)
    return list(accS)

def chunks(l, n):
    '''Yield successive n-sized chunks from l.
            
    code from stackoverflow.
    '''
    for i in range(0, len(l), n):
        yield l[i:i+n]

class Conversion:
    def __init__(self, ids, *columns):
        self.u = UniProt(verbose=False)
        self.inputAccL = ids # ids is a list of UniProt accessions
        self.colstr = ', '.join(columns) # list of values, separated by commas

        print self.inputAccL[:5]
        print self.colstr

    def chunks(self, l, n):
        '''Yield successive n-sized chunks from l.
            
        code from stackoverflow.
        '''
        for i in range(0, len(l), n):
            yield l[i:i+n]

    def writeString(self, resstr, filename):
        '''Write output from request to file with len(columns) columns.'''
        #a = re.sub('\n', ',', reqstr).split(',')[1:-1]
        with open(filename, 'a') as of:
            of.write(resstr + '\n')
            
    def submit(self):
        l = []
        
        isheader = 1
        for acc in self.inputAccL:
            ret = self.u.search(acc + ', ' + 'organism:9606', columns=self.colstr, frmt='tab', limit=None)
            retl = ret.split('\n')
            
            if isheader:
                l.append(retl[0])
                isheader = 0

            if len(retl) > 1:
                l.append(retl[1])
                time.sleep(5)
            else:
                l.append('No entry for %s' % acc)
        return l

    def write2file(self, strlst, fn):
        for strng in strlst:
            self.writeString(strng, fn)

#u.search('Q9H9S0, organism:9606', columns='entry name, id, genes(PREFERRED), genes, protein names',
#         frmt='tab', limit=None)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        usage ='''Usage: %s <One or more UniProt acc files>
        Example: python getProtannot_v0UniProt.py ../TempData/x227.txt ../TempData/x238.txt etc
        ''' % sys.argv[0]
        print usage
        sys.exit(0)

    flst = sys.argv[1:]
    accL = _getAcc('Accession', flst)

    counter = 1
    for L in chunks(accL[299:], 50):
        cv = Conversion(L, 'id', 'genes(PREFERRED)', 'genes', 'protein names')
        res = cv.submit()
        print '%d accessions are done\n' % counter * 50
        counter = counter + 1
        time.sleep(30)
        cv.write2file(res, 'idAnnotWithUniProt')
    
