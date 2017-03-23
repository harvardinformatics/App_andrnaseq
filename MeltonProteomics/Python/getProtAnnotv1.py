import urllib
import urllib2
import sys
import re

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

def str2Dict(reqstr):
    '''Make dictionary from Request.submit() string.'''
    a = re.sub('\n', ',', reqstr).split(',')[1:-1]
    b = [tuple(x.split('\t')) for x in a]
    return dict(b)

def writeString(reqstr, filename):
    '''Write output from request to file with 2 columns.'''
    a = re.sub('\n', ',', reqstr).split(',')[1:-1]
    with open(filename, 'w') as of:
        for el in a:
            of.write(el + '\n')

class Request:
    def __init__(self, acc, ktype, attr):
        self.url = 'http://www.uniprot.org/uploadlists/'
        self.contact = "jstraubhaar@fas.harvard.edu"
        self.params =  {
            'from': ktype, # ACC
            'to':attr, # GENENAME, ID
            'format':'tab',
            'query':acc
        }

    def submit(self):
        data = urllib.urlencode(self.params)
        req = urllib2.Request(self.url, data)
        req.add_header = ('User-Agent', 'Python %s' % self.contact)
        resp = urllib2.urlopen(req)
        return resp.read()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        use = sys.argv[0] + ' <file name: sample file>'
        print use
        sys.exit(1)

    l = getAcc(sys.argv[1], 'Accession') # getProtAnnotv1.py TempData/x238.txt (tab delimited)
    s = ' '.join(l[:5])

    rq = Request(s, 'ACC', 'ID')
    res = rq.submit()
    
    resD = str2Dict(res)
    writeString(res, 'x238ProtDict')    

