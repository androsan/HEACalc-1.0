""" Podporni modul programa HEA-Calc 1.0, ki zajema iskanje korelacij med termodinamskimi
parametri za hitrejše odkrivanje novih HEA zlitin ---> analogija z brezdimenzijskimi števili"""
import math, statistics, pari, branje_csv, re
from hea_modeliranje_AI import normalize as normalize
from hea_modeliranje_AI import parske_entalpije as parske_entalpije
parske_kombinacije=branje_csv.parske_kombinacije
molske_mase=branje_csv.molske_data
density=branje_csv.density_data

def Rkvadrat(alfa,beta):
    a=parske_entalpije(alfa[beta][-2])
    h=[i[1] for i in a.items()]
    Hmix=sum(h)/len(h)
    Rsq=sum([(j-Hmix)**2 for j in h])
    Hfractions=[k*100/Hmix for k in h]
    print(alfa[beta][-3], Hfractions)
    return math.sqrt(Rsq)

def STD_parskih_funkcija(alfa,beta):
    c=parske_entalpije(alfa[beta][-2])
    h=[i[1] for i in c.items()]
    STD = round(statistics.stdev(h), 3)
    return STD

def gostota_funkcija(alfa,beta):
    vnos=normalize(alfa[beta][-2])
    im=sum([vnos[i]*molske_mase[i] for i in vnos])
    st=sum([vnos[i]*molske_mase[i]/density[i] for i in vnos])
    rho_mix=im/st
    return rho_mix
    

# funkcija I =========================>>> X-axis
def enka(baza, n):
    delta=baza[n][0] ; chi=baza[n][1] ; VEC=baza[n][2] ; dHmix=baza[n][3] ; dSmix=baza[n][4] ; Tm=baza[n][5] ; omega=baza[n][6] ; e_a=baza[n][7] ; molska=baza[n][8] ; d_mol=baza[n][9] ; Ssigma=baza[n][10] ; ksi=baza[n][11]
    try:
        #STD=STD_parskih_funkcija(baza,n)
        a='(delta+chi)**0.5/VEC'
        b='(delta+chi)**0.5/VEC  '    #'STD_parskih_funkcija*delta'
        #fi =  (dSmix/8.314-abs(dHmix)/(Tm+273))/abs(Ssigma)
        c='math.e**(-dHmix/(8.314*(Tm+273))+dSmix/8.314)'
        d='delta**2'
        
        string_x= a  

        X  =   eval(string_x) 
    except ZeroDivisionError:
        X=0
    return X, string_x

# funkcija II ========================>>> Y-axis
def dvojka(baza, n):
    delta=baza[n][0] ; chi=baza[n][1] ; VEC=baza[n][2] ; dHmix=baza[n][3] ; dSmix=baza[n][4] ; Tm=baza[n][5] ; omega=baza[n][6] ; e_a=baza[n][7] ; molska=baza[n][8]; d_mol=baza[n][9] ; Ssigma=baza[n][10] ; ksi=baza[n][11]
    try:
            gostota=gostota_funkcija(baza,n)
            #STD=STD_parskih_funkcija(baza,n)
            a= '(VEC+gostota)'
            b='(VEC+gostota)*dHmix/10'
            c= 'math.e**(-dHmix/(8.314*(Tm+273))+dSmix/8.314)'
            d= 'dHmix'

            string_y= a
            
            Y =   eval(string_y) 
            #Rkvadrat(baza,n)
            #(8.314*delta**2)/dSmix          

    except ZeroDivisionError or ValueError:
        Y=0
    return Y, string_y

