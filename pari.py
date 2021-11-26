#['O','F','Mg','Al','Si','S','Ca']   # Vstopni podatki, vnesi elemente, katerih kombinacije te zanimajo..
insert= {'Co':1, 'Cr':1, 'Fe':1, 'Ni':1}

def COMBINATIONS(x):
    elementi=list(x.keys())
    test=['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'r', 's', 't', 'u', 'v', 'z']

    for i in range(len(test)-len(elementi)):
        del test[-1]

    def add_element(baza_A, baza_B):
        get=[]
        for i in baza_A:
            for j in baza_B:
                if i!=j:
                    get.append(i+j)
        return get
     
    def erase_duplicates(vrsta):
        y=[]
        for i in vrsta:
            y.append(set(i))
        y_new=[]
        for i in y:
            if i not in y_new:
                y_new.append(i)
        for j,i in enumerate(y_new):
            y_new[j]= "".join(i)

        return y_new

    def delete_same_letters(raw, n):    # n= število elementov
        y=[]
        for i in raw:
            y.append(set(i))
        y_new=[]
        for i in y:
            if len(list(i))==n:
                y_new.append(i)
        for j,i in enumerate(y_new):
            y_new[j]= "".join(i)
        return y_new

    rezultat=[]
    def find_combinations(baza):
        pairs_raw= add_element(baza,baza)
        pairs=erase_duplicates(pairs_raw)
        rezultat.append(pairs)
        t=-1
        while t<len(baza)-3:
            t=t+1
            xxx= add_element(rezultat[t],baza)
            xx=delete_same_letters(xxx,t+3)
            x=erase_duplicates(xx)
            rezultat.append(x)

    find_combinations(test)
    #print(rezultat)

    novo=[]
    for i in rezultat:
        novo_r=[]
        for j in i:
            qu= list(set(j))
            for kk, k in enumerate(qu):
                g=test.index(k)
                qu[kk] = elementi[g]
                m="".join(qu)
            novo_r.append(m)
        novo.append(novo_r)

    return novo


        




