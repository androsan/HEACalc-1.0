import re, itertools
from machine_learning_experimental import c

def CONVERT(raw):
    u=list(raw)
    for j in range(50):   
        try:
            a=u[j]
            b=u[j+1]
            if a.isalpha() and a.isupper() and b.isdigit():   # doda 'x' med npr. V in cifro
                u[j+1:j+1]=['x']
            elif a.isalpha() and a.isupper() and b.isalpha()and b.isupper():   # doda 'x1' med npr. V in npr. Cr
                u[j+1:j+1]=['x1']
            elif a.isalpha() and a.islower() and b.isalpha()and b.isupper():   # doda '1' med npr. Co in npr. Cr
                u[j+1:j+1]=['1']
        except IndexError:
            pass
    if u[-1].isalpha()and u[-1].isupper():
        u.append('x')
    o="".join(u)

    def split(s):
        res = []
        def replace(matchobj):
            res.append(matchobj.group(0))
            return ''

        letter = re.compile('^([A-Z]+[a-z]+)')
        number = re.compile('^(\.\d|\d+\.\d+|\d+)')

        if letter.match(s):
            c = itertools.cycle([letter, number])
        else:
            c = itertools.cycle([number, letter])

        for op in c:
            mods = op.sub(replace, s)
            if len(s) == len(mods):
                #print('mods: ',mods)
                return
            elif not mods:
                #print('res: ',res)
                return res
            s = mods

    res=split(o)

    keys=[]; vals=[]; d={}
    for i in range(len(res)):
        if i%2==0:
            keys.append(res[i])     
    for i in range(len(res)):
        if i%2!=0:
            vals.append(float(res[i]))
    if len(vals)< len(keys):
        vals.append(1)
    for i,j in zip(keys,vals):
        d[i]=j

    for i in d:
        if 'x' in i:
            try:
                f=i.replace('x','')
                d[f]=d.pop(i)
            except RuntimeError:
                pass
    return d

def data_from_literature():
    datki=[]
    for i in c:
        g=CONVERT(i[0])
        datki.append([g, i[1]])
    return datki

    
