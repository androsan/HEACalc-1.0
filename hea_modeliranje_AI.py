""" Ta program je napisan za iskanje novih visokoentropijskih zlitin, skozi modeliranje termodinamskih
     parametrov in stabilnosti. Ime programa: HEA-CALC 1.0
     - - -                                                       - - - - -                                                - - -                              
     * * * Avtor: dr. Andraz Kocjan, Inštitut za kovinske materiale in tehnologije, december 2017. * * *
            - - -                                                - - - - -                                                - - -
---------------------------------------------------------------------------------------------------------------------------------
Seznam parametrov:

            COMMON:
            -----------------------------------------------
           delta    -  atomic size mismatch (%)
         kappa    -  electronegativity mismatch (%)
           VEC    -  valence electron concentration (-)
             e/a    -  number of valence electrons per number of atoms in the unit cell (-)
         dHmix    -  enthalpy of mixing (kJ/mol)
         dSmix    -  configurational entropy of mixing (J/mol)
             Tm    -  melting point (Kelvins)
        omega    - (-)
            -----------------------------------------------------
            Less Common:
            
          molska_pop   -  molar mass of an alloy (g/mol)      
     molske_sipanje   -  molar mass mismatch (%)
                dSsigma   -  mismatch entropy (J/mol)
            -----------------------------------------------------   
_________________________________________________
                                                                                                                                               """
import sys
sys.path.insert(0, '/path/to/C:/Users/akocjan/Desktop/Python programi/modeliranje HEA')
import math, re
import branje_csv, pari
import machine_learning_lesson
import matplotlib.pyplot as plt
from argparse import Namespace


def normalize(y):
    d={}
    for j in y:
        d[j]=y[j]
    vrsta= [d[i] for i in d]
    vsota= sum(vrsta)
    for i in d:
        d[i]= d[i]/vsota
    return d

"""-----------------------------------------------------------------
* * * B A Z A  podatkov o atomih (branje_csv.py oziroma podatki_o_atomih.csv) * * *
__________________________
                                                                             """
radii=branje_csv.radii
#kappa_data=branje_csv.kappa_data_Pauling
kappa_data=branje_csv.kappa_data_Allen
VEC_data=branje_csv.VEC_data
Tm_data=branje_csv.Tm_data
e_a_data=branje_csv.e_a_data
molske_data=branje_csv.molske_data
"""-----------------------------------------------------------------"""

#  Osnovne funkcije (delta, kappa, VEC, dHmix, dSmix, Tm)

def delta(x):                                                   # ------------------------------------- DELTA [%] 
    x=normalize(x)
    vrsta_avrg = [x[i]*radii[i] for i in x]
    avrg = sum(vrsta_avrg)
    vrsta_delta = [x[i]*(1-radii[i]/avrg)**2 for i in x]
    delta = 100*math.sqrt(sum(vrsta_delta))
    return round(delta,3)                                     

def kappa(x):                                                  # ------------------------------------ KAPPA [ %]
    x=normalize(x)
    vrsta_avrg = [x[i]*kappa_data[i] for i in x]
    avrg = sum(vrsta_avrg)
    vrsta_kappa = [x[i]*(1-kappa_data[i]/avrg)**2 for i in x]
    kappa = 100*math.sqrt(sum(vrsta_kappa))
    return round(kappa,3)

def VEC(x):                                                     # ------------------------------------ VEC [ -]
    x=normalize(x)
    vrsta_VEC=[x[i]*VEC_data[i] for i in x]
    VEC = sum(vrsta_VEC)
    return round(VEC,3)

def e_a(x):                                                     # ------------------------------------ e/a [ -]
    x=normalize(x)
    vrsta_e_a=[x[i]*e_a_data[i] for i in x]
    e_a = sum(vrsta_e_a)
    return round(e_a,3)

"""-----------------------------------------------------------------
* * * B A Z A  podatkov o parskih entalpijah mešanja
       (pari.py oziroma parske_entalpije.csv) * * *
__________________________
                                                                             """
parske_kombinacije=branje_csv.parske_kombinacije

"""-----------------------------------------------------------------"""

def dHmix(x):                                                   # ------------------------------------ dHmix [ kJ/mol]
    global wwww
    x=normalize(x)
    wwww=pari.COMBINATIONS(x)
    vrsta_Hmix = []
    for i in wwww[0]:
        aa= re.sub( r"([A-Z])", r" \1", i).split()
        h_mix = 4*parske_kombinacije[i]*x[aa[0]]*x[aa[1]]
        vrsta_Hmix.append(h_mix)
    dHmix = sum(vrsta_Hmix)
    return round(dHmix,3)


def dSmix(x):                                                   # ------------------------------------ dSmix [ J/mol]
    x=normalize(x)
    vrsta_Smix=[]
    for i in x:
        if x[i] != 0:
            vrsta_Smix.append(x[i]*math.log(x[i]))
    dSmix=-8.314*sum(vrsta_Smix)
    return round(dSmix,3)

def Tm(x):                                                        # ------------------------------------ Tm [ K]
    x=normalize(x)
    vrsta_Tm=[x[i]*Tm_data[i] for i in x]
    Tm = sum(vrsta_Tm)
    return round(Tm,)
    
def omega(x):                                                   # ------------------------------------ omega [ -]
    if dHmix(x)!= 0:
        omega=Tm(x)*dSmix(x)/(abs(dHmix(x))*1000)
        return round(omega,3)
    else:
        return 'dHmix is zero..'


def molske_sipanje(x):                                      # ------------------------------------ molske [%] 
    x=normalize(x)
    vrsta_avrg = [x[i]*molske_data[i] for i in x]
    avrg = sum(vrsta_avrg)
    vrsta_molske = [x[i]*(1-molske_data[i]/avrg)**2 for i in x]
    molske_sipanje = 100*math.sqrt(sum(vrsta_molske))
    return round(avrg,3), round(molske_sipanje,3)


def dSsigma(x, ksi_in):
    # ksi ===>>> atomic packing fraction, amrphous = 0.64, BCC=0.68, FCC&HCP=0.74
    x=normalize(x)
    zeta=1/(1-ksi_in)
    oooo=pari.COMBINATIONS(x)
    vrsta_y1 = []; vrsta_y2 = []
    for i in oooo[0]:
        aa= re.sub( r"([A-Z])", r" \1", i).split()
        y1_part = (radii[aa[0]]+radii[aa[1]])*(radii[aa[0]]-radii[aa[1]])**2*x[aa[0]]*x[aa[1]]
        vrsta_y1.append(y1_part)
        y2_part = radii[aa[0]]*radii[aa[1]]*(radii[aa[0]]-radii[aa[1]])**2*x[aa[0]]*x[aa[1]]
        vrsta_y2.append(y2_part)

    def sigma(p):
        sigma_vrsta=[x[i]*radii[i]**p for i in x]
        sigma=sum(sigma_vrsta)
        return sigma
     
    y1=1/sigma(3)* sum(vrsta_y1)
    y2=sigma(2)/(sigma(3))**2 * sum(vrsta_y2)
    y3=(sigma(2))**3/(sigma(3))**2
    
    dSs_kB = 1.5*(zeta**2-1)*y1 + 1.5*(zeta-1)**2*y2-(0.5*(zeta-1)*(zeta-3)+ math.log(zeta))*(1-y3)
    return round(dSs_kB,3)

def Fi(Sc, hmix, tm, Sexess):
    fi=(Sc/8.314-abs(hmix)/(tm+273))/Sexess
    return fi


""" GUI Tkinter === Uporabniski Graficni Vmesnik """
from tkinter import *

root=Tk()
root.geometry('1536x801+0+0')  # <<< ----  za en monitor
#root.geometry('1536x801+-1928+-8')# <<< ----  za dva monitorja
root.title('HEA-Calc 1.0 by Andro'+u'\u00AE')

# ksi (variable packing)is atomic packing fraction, i.e. 0.64 for amorphous, 0.68 for BCC and 0.74 for FCC
packing=DoubleVar()
packing.set(0.64)
Label(root, text='Atomic packing fraction ').grid(row=20, column=1, columnspan=2)
pck_ent=Entry(root, textvariable=packing, justify='center', width=6); pck_ent.grid(row=20, column=3)

def splosni_izracun(vnos, ksi):
        c= normalize(vnos)
        atomic_perc={}
        try:
            z_0= delta(c)
            z_1= kappa(c)
            z_2= VEC(c)
            z_3= dHmix(c)
            z_4= dSmix(c)
            z_5= Tm(c)
            z_6=omega(c)
            z_7=e_a(c)
            z_8=molske_sipanje(c)
            z_9=dSsigma(c, ksi)
            try:
                AI = machine_learning_lesson.insert_HEA_parameters(z_0, z_3, z_1, z_6, z_2)
                inteligence.set(AI[0])
            except ValueError:
                inteligence.set(' - ')
            a= list(c.keys())
            b= list(c.values())
            b=[u*100 for u in b]
            for i,j in zip(a,b):
                atomic_perc[i]=j
            mass_perc={}
            molske=[]
            for i in atomic_perc:
                val=atomic_perc[i]*molske_data[i]/100
                molske.append(val)
            ppp=sum(molske)
            for i in atomic_perc:
                val=atomic_perc[i]*molske_data[i]/ppp
                mass_perc[i]=val
            formula_at = " ".join([i+'-'+str(round(j,1))for i,j in zip(a,b)])
            am=list(mass_perc.keys())
            bm=list(mass_perc.values())
            formula_mass = " ".join([i+'-'+str(round(j,1))for i,j in zip(am,bm)])
            if fn.get()==0:
                res_formula.set('         '+formula_at+'    in at.%')
            elif fn.get()==1:
                res_formula.set('         '+formula_mass+'   in wt.%')
            res_delta.set(z_0); res_kappa.set(z_1)
            res_VEC.set(z_2); res_dHmix.set(z_3);res_dSmix.set(z_4)
            res_Tm.set(z_5); res_omega.set(z_6); res_e_a.set(z_7)
            if z_4 >= 12.5:
                S_R=z_4/8.314
                dSmix_comment.config(fg='#ff9933')
                op_dSmix.set('high entropy: '+"%.2f" % round(S_R, 2)+'R')
                # DELTA comments----------------------------
                if z_0 <= 4:
                    delta_comment.config(fg='green')
                    op_delta.set('HEA')
                elif z_0 > 4 and z_0 <= 6.6:
                    delta_comment.config(fg='blue')
                    op_delta.set('HEA-IM')
                elif z_0 > 6.6 and z_0 <= 12:
                    delta_comment.config(fg='blue')
                    op_delta.set('IM-BMG')
                elif z_0 > 12:
                    delta_comment.config(fg='red')
                    op_delta.set('BMG')
                # KAPPA comments----------------------------
                if z_1 <= 10.5:
                    kappa_comment.config(fg='green')
                    op_kappa.set('HEA')
                elif z_1 > 10.5:
                    kappa_comment.config(fg='blue')
                    op_kappa.set('IM-BMG')
                # VEC comments--------------------------------
                if z_2 <=6.88:
                    op_VEC.set('BCC')
                elif z_2 > 6.88 and z_2 <= 8:
                    op_VEC.set('BCC+FCC')
                elif z_2 > 8:
                    op_VEC.set('FCC')
                # dHmix comments-----------------------------
                if z_3 <= 5 and z_3 >= -5:
                    dHmix_comment.config(fg='green')
                    op_dHmix.set('HEA')
                elif z_3 > 5:
                    dHmix_comment.config(fg='red')
                    op_dHmix.set('seggregation')
                elif z_3 < -5 and z_3 >= -12:
                    dHmix_comment.config(fg='blue')
                    op_dHmix.set('HEA-IM')
                elif z_3 < -12 and z_3 >= -25:
                    dHmix_comment.config(fg='blue')
                    op_dHmix.set('IM-BMG')
                elif z_3 < -25:
                    dHmix_comment.config(fg='red')
                    op_dHmix.set('BMG')
                # OMEGA comments---------------------------
                try:
                    if z_6 >= 10:
                        omega_comment.config(fg='green')
                        op_omega.set('HEA ')
                    elif z_6 < 10 and z_6 >= 1.1 and z_0 <= 6.6:
                        omega_comment.config(fg='blue')
                        op_omega.set('HEA-IM')
                    elif z_6 < 10 and z_6 >= 1.1 and z_0 > 6.6:
                        omega_comment.config(fg='red')
                        op_omega.set('IM')
                    elif z_6 < 1.1:
                        omega_comment.config(fg='red')
                        op_omega.set('BMG')
                    elif z_6 =='dHmix is zero..':
                        omega_comment.config(fg='green')
                        op_omega.set('infinite')
                except TypeError:
                    pass
            elif z_4 < 12.5:
                dSmix_comment.config(fg='#ff80bf')
                op_dSmix.set('entropy too low..')
                op_delta.set(''); op_kappa.set(''); op_VEC.set(''); op_dHmix.set(''); op_omega.set('')
                
        except ValueError or ZeroDivisionError:
            pass

        return z_0 ,z_1 ,z_2,z_3 ,z_4 ,z_5 ,z_6 ,z_7,z_8[0],z_8[1],z_9, ksi, formula_at, atomic_perc, 'calculated'
            

e1var =StringVar(); e2var=StringVar() ; e3var=StringVar(); e4var=StringVar(); e5var=StringVar(); e6var=StringVar(); e7var=StringVar(); e8var=StringVar()
e1num=DoubleVar(); e2num=DoubleVar(); e3num=DoubleVar(); e4num=DoubleVar(); e5num=DoubleVar(); e6num=DoubleVar(); e7num=DoubleVar(); e8num=DoubleVar()

e_width=5
ee_width=5
font_e1='Arial 16 bold'
font_e11='Arial 16'
font_sub_buttons='Arial 12'

Label(root,text='           ').grid(row=0, column=0)
Label(root,text='Substitution   ').grid(row=1, column=1, sticky=E)
Label(root,text='Symbol   ').grid(row=2, column=1,sticky=E); Label(root,text='Fraction   ').grid(row=3, column=1, sticky=E)
button1x= Button(root, text='x', font=font_sub_buttons, command=lambda: pick_substitution(button1x,1),width=2);                   e1=Entry(root, textvariable=e1var, font=font_e1, width=e_width, justify='center'); e11=Entry(root, textvariable=e1num, font=font_e11,width=ee_width, justify='center')
button1y= Button(root, text='y', font=font_sub_buttons, command=lambda: pick_substitution(button1y,1),width=2);                   e2=Entry(root, textvariable=e2var, font=font_e1,width=e_width, justify='center'); e21=Entry(root, textvariable=e2num, font=font_e11,width=ee_width, justify='center')
button2x= Button(root, text='x', font=font_sub_buttons, command=lambda: pick_substitution(button2x,2),width=2);                   e3=Entry(root, textvariable=e3var, font=font_e1,width=e_width, justify='center'); e31=Entry(root, textvariable=e3num, font=font_e11,width=ee_width, justify='center')
button2y= Button(root, text='y', font=font_sub_buttons, command=lambda: pick_substitution(button2y,2),width=2);                   e4=Entry(root, textvariable=e4var, font=font_e1,width=e_width, justify='center'); e41=Entry(root, textvariable=e4num, font=font_e11,width=ee_width, justify='center')
button3x= Button(root, text='x', font=font_sub_buttons, command=lambda: pick_substitution(button3x,3),width=2);                   e5=Entry(root, textvariable=e5var, font=font_e1,width=e_width, justify='center'); e51=Entry(root, textvariable=e5num, font=font_e11,width=ee_width, justify='center')
button3y= Button(root, text='y', font=font_sub_buttons, command=lambda: pick_substitution(button3y,3),width=2);                   e6=Entry(root, textvariable=e6var, font=font_e1,width=e_width, justify='center'); e61=Entry(root, textvariable=e6num, font=font_e11,width=ee_width, justify='center')
button4x= Button(root, text='x', font=font_sub_buttons, command=lambda: pick_substitution(button4x,4),width=2);                   e7=Entry(root, textvariable=e7var, font=font_e1,width=e_width, justify='center'); e71=Entry(root, textvariable=e7num, font=font_e11,width=ee_width, justify='center')
button4y= Button(root, text='y', font=font_sub_buttons, command=lambda: pick_substitution(button4y,4),width=2);                   e8=Entry(root, textvariable=e8var, font=font_e1,width=e_width, justify='center'); e81=Entry(root, textvariable=e8num, font=font_e11,width=ee_width, justify='center')
button5x= Button(root, text='x', font=font_sub_buttons, command=lambda: pick_substitution(button5x,5),width=2);  button6x= Button(root, text='x', font=font_sub_buttons, command=lambda: pick_substitution(button6x,6),width=2)
button5y= Button(root, text='y', font=font_sub_buttons, command=lambda: pick_substitution(button5y,5),width=2);  button6y= Button(root, text='y', font=font_sub_buttons, command=lambda: pick_substitution(button6y,6),width=2)
button7x= Button(root, text='x', font=font_sub_buttons, command=lambda: pick_substitution(button7x,7),width=2);  button8x= Button(root, text='x', font=font_sub_buttons, command=lambda: pick_substitution(button8x,8),width=2)
button7y= Button(root, text='y', font=font_sub_buttons, command=lambda: pick_substitution(button7y,7),width=2);  button8y= Button(root, text='y', font=font_sub_buttons, command=lambda: pick_substitution(button8y,8),width=2)

X_sub={}; Y_sub={}

for i in range(1,9):
    X_sub["button{}x".format(i)]=False
    Y_sub["button{}y".format(i)]=False

def pick_substitution(a,cifra):
    global X_sub
    if a['relief']=='sunken':
        if a['text']=='x':
            X_sub["button{}x".format(cifra)]=False
        else:
            Y_sub["button{}y".format(cifra)]=False
        a.config(relief=RAISED, bg=root.cget('bg'))
    else:
        if a['text']=='x':
            X_sub["button{}x".format(cifra)]=True
            a.config(relief=SUNKEN, bg='#00cc99')
        else:
            Y_sub["button{}y".format(cifra)]=True
            a.config(relief=SUNKEN, bg='#cc6699')
        
    return X_sub, Y_sub

padeks=4
button1x.grid(column=2,row=1,sticky=W,padx=padeks); button1y.grid(column=2,row=1,sticky=E,padx=padeks)
button2x.grid(column=3,row=1,sticky=W,padx=padeks); button2y.grid(column=3,row=1,sticky=E,padx=padeks)
button3x.grid(column=4,row=1,sticky=W,padx=padeks); button3y.grid(column=4,row=1,sticky=E,padx=padeks)
button4x.grid(column=5,row=1,sticky=W,padx=padeks); button4y.grid(column=5,row=1,sticky=E,padx=padeks)
button5x.grid(column=6,row=1,sticky=W,padx=padeks); button5y.grid(column=6,row=1,sticky=E,padx=padeks)
button6x.grid(column=7,row=1,sticky=W,padx=padeks); button6y.grid(column=7,row=1,sticky=E,padx=padeks)
button7x.grid(column=8,row=1,sticky=W,padx=padeks); button7y.grid(column=8,row=1,sticky=E,padx=padeks)
button8x.grid(column=9,row=1,sticky=W,padx=padeks); button8y.grid(column=9,row=1,sticky=E,padx=padeks)

init_x=StringVar(); init_y=StringVar(); init_xy=StringVar()

init_frame = Frame(root, width=650)
init_frame.grid(row=1, column=10, columnspan=30, sticky=W, padx=30)
init_frame.config(borderwidth=1, relief=SUNKEN)
zulu=('Arial 18')
lab_initx = Label(init_frame, text='  Init X '); lab_inity = Label(init_frame,text='  Init Y '); lab_initxy = Label(init_frame,text='  Init XY ')
lab_initx.grid(row=0, column=0, sticky=E); lab_inity.grid(row=0, column=2, sticky=E); lab_initxy.grid(row=0, column=4, sticky=E)
initx_entry=Entry(init_frame ,textvariable=init_x, font=zulu, justify='center', width=5); inity_entry=Entry(init_frame,textvariable=init_y, font=zulu, justify='center', width=5); initxy_entry=Entry(init_frame,textvariable=init_xy, font=zulu, justify='center', width=5)
initx_entry.grid(row=0, column=1, sticky=W); inity_entry.grid(row=0, column=3, sticky=W); initxy_entry.grid(row=0, column=5, sticky=W); 

e1.grid(column=2, row=2); e11.grid(column=2, row=3)
e2.grid(column=3, row=2); e21.grid(column=3, row=3)
e3.grid(column=4, row=2); e31.grid(column=4, row=3)
e4.grid(column=5, row=2); e41.grid(column=5, row=3)
e5.grid(column=6, row=2); e51.grid(column=6, row=3)
e6.grid(column=7, row=2); e61.grid(column=7, row=3)
e7.grid(column=8, row=2); e71.grid(column=8, row=3)
e8.grid(column=9, row=2); e81.grid(column=9, row=3)

entries=[e1,e2,e3,e4,e5,e6,e7,e8]
numbers=[e1num,e2num,e3num,e4num,e5num,e6num,e7num,e8num]
             
def get_data():
    global k, varjable
    k={}
    varjable=[e1var.get(),e2var.get(),e3var.get(),e4var.get(),e5var.get(),e6var.get(),e7var.get(),e8var.get()]
    cifre=[e1num.get(),e2num.get(),e3num.get(),e4num.get(),e5num.get(),e6num.get(),e7num.get(),e8num.get()]
    for i,j in zip(varjable,cifre):
        try:
            if i[0].islower() is True:
                lj=i.replace(i[0],i[0].upper())
                k[lj]=j
            else:
                k[i]=j
        except IndexError:
            pass
    return k

def get_subs():
    input_x={}; input_y={}; input_xy={}; input_same={}
    varjable=[e1var.get(),e2var.get(),e3var.get(),e4var.get(),e5var.get(),e6var.get(),e7var.get(),e8var.get()]
    w=get_data()
    for i in varjable:
        if bool(i)==True and X_sub["button{}x".format(varjable.index(i)+1)] ==True and Y_sub["button{}y".format(varjable.index(i)+1)] ==False:
            try:
                input_x[i]=w[i]
            except KeyError:
                lj=i.replace(i[0],i[0].upper())
                input_x[lj]=w[lj]
                
        elif bool(i)==True and Y_sub["button{}y".format(varjable.index(i)+1)] ==True and X_sub["button{}x".format(varjable.index(i)+1)] ==False:
            try:
                input_y[i]=w[i]
            except KeyError:
                lj=i.replace(i[0],i[0].upper())
                input_y[lj]=w[lj]

        elif bool(i)==True and X_sub["button{}x".format(varjable.index(i)+1)] ==True and Y_sub["button{}y".format(varjable.index(i)+1)] ==True:
            try:
                input_xy[i]=w[i]
            except KeyError:
                lj=i.replace(i[0],i[0].upper())
                input_xy[lj]=w[lj]

        elif bool(i)==True and X_sub["button{}x".format(varjable.index(i)+1)] ==False and Y_sub["button{}y".format(varjable.index(i)+1)] ==False:
            try:
                input_same[i]=w[i]
            except KeyError:
                lj=i.replace(i[0],i[0].upper())
                input_same[lj]=w[lj]

    return input_x, input_y, input_xy, input_same

res_formula= StringVar();res_delta= StringVar();res_kappa=StringVar()
res_VEC=StringVar();res_dHmix=StringVar();res_dSmix=StringVar()
res_Tm=StringVar();res_omega=StringVar();res_e_a=StringVar()

op_delta=StringVar(); op_kappa=StringVar(); op_VEC=StringVar(); op_dHmix=StringVar(); op_omega=StringVar(); op_dSmix=StringVar()

fn= IntVar() # --- fn = function
fn.set(0)      # --- atomski %(default)

nosilec_podatkov=[]; dolzine=[]; one_shot=[]
nosilec_podatkov_status=StringVar()
nosilec_podatkov_status.set('NO DATA')

def stevilo_zlitin(numero):  # numero = stevilo shranjenih zlitin po posamezni iteraciji
    global dolzine
    dolzine.append(numero)
    stzlit=sum(dolzine)
    if stzlit==1:
        nosilec_podatkov_status.set(str(stzlit)+' alloy to show..')
    else:
        nosilec_podatkov_status.set(str(stzlit)+' alloys to show..')



def calculate():
    global mass_perc, one_shot
    w=get_data()
    odgovor=splosni_izracun(w, packing.get())
    one_shot.append([odgovor])
    nosilec_podatkov.append((one_shot, 'single_case'))
    stevilo_zlitin(1)
    nps.config(fg='lime')
    #print((one_shot, 'single_case'))
    #single_shot_plot((one_shot, 'single_case'))
    
Label(root,text='').grid(row=4, column=1)
formjula=Label(root,textvariable=res_formula, font=font_e1, fg='#1e90ff')
formjula.grid(row=3, column=10, sticky=W, columnspan=30)
Label(root,text='                                  '+u"\u0394".lower()+' =   ', font=font_e1).grid(row=6, column=10, sticky=E); Label(root,textvariable=res_delta, font=font_e11, width=7).grid(row=6, column=11, sticky=W)
Label(root,text='                                  '+u"\u03C7".lower()+' =   ', font=font_e1).grid(row=7, column=10, sticky=E); Label(root,textvariable=res_kappa, font=font_e11, width=7).grid(row=7, column=11, sticky=W)
Label(root,text='                            VEC =   ', font=font_e1).grid(row=8, column=10, sticky=E); Label(root,textvariable=res_VEC, font=font_e11, width=7).grid(row=8, column=11, sticky=W)
Label(root,text='                         '+u"\u0394"+'Hmix =   ', font=font_e1).grid(row=9, column=10, sticky=E); Label(root,textvariable=res_dHmix, font=font_e11, width=7).grid(row=9, column=11, sticky=W)
Label(root,text='                         '+u"\u0394"+'Smix =   ', font=font_e1).grid(row=10, column=10, sticky=E); Label(root,textvariable=res_dSmix, font=font_e11, width=7).grid(row=10, column=11, sticky=W)
Label(root,text='                               Tm =   ', font=font_e1).grid(row=11, column=10, sticky=E); Label(root,textvariable=res_Tm, font=font_e11, width=7).grid(row=11, column=11, sticky=W)
Label(root,text='                         '+u"\u03a9"+' =   ', font=font_e1).grid(row=12, column=10, sticky=E); Label(root,textvariable=res_omega, font=font_e11, width=7).grid(row=12, column=11, sticky=W)
Label(root,text='                         e/a =   ', font=font_e1).grid(row=13, column=10, sticky=E); Label(root,textvariable=res_e_a, font=font_e11, width=7).grid(row=13, column=11, sticky=W)

delta_comment=Label(root, textvariable=op_delta, font=font_e11); delta_comment.grid(row=6, column=12, sticky=W, padx=20)
kappa_comment=Label(root, textvariable=op_kappa, font=font_e11); kappa_comment.grid(row=7, column=12, sticky=W, padx=20)
VEC_comment=Label(root, textvariable=op_VEC, font=font_e11); VEC_comment.grid(row=8, column=12, sticky=W, padx=20)
dHmix_comment=Label(root, textvariable=op_dHmix, font=font_e11); dHmix_comment.grid(row=9, column=12, sticky=W, padx=20)
dSmix_comment=Label(root, textvariable=op_dSmix, font=font_e11, width=19); dSmix_comment.grid(row=10, column=12, sticky=W)
omega_comment=Label(root, textvariable=op_omega, font=font_e11); omega_comment.grid(row=12, column=12, sticky=W, padx=20)

# Machine-learning part
inteligence=StringVar()
intel=Label(root, text='AI predicted structure:   ', font='Arial 14', fg='grey'); intel.grid(row=14, column=10, sticky=E)
intel_com=Label(root, textvariable=inteligence, font='Arial 14 bold', fg='grey', width=7);intel_com.grid(row=14, column=11, sticky=W)

Label(root,text='').grid(row=19,column=2)    
Button(root, text='Calculate', command=calculate, font=('44')).grid(row=20, column=4, columnspan=2)


""" MANUAL CONTROL of ATOMIC FRACTIONS - sliders.. """

Label(root,text='To').grid(row=5, column=1)
Label(root,text='From').grid(row=18, column=1)

Label(root,text='STEP', fg='blue').grid(row=8,column=1, sticky=S)
div_fact=DoubleVar(); div_fact.set(0.01); razpon = 100
Entry(root, textvariable=div_fact, width=5, fg='blue', justify='center').grid(row=9, column=1, sticky=N)

def slajder_1(value):
    vr=abs(float(value))*div_fact.get()* eslfrom_1.get()/100
    e1num.set(round(vr,2))
    calculate()

def slajder_2(value):
    vr=abs(float(value))*div_fact.get()* eslfrom_2.get()/100
    e2num.set(round(vr,2))
    calculate()

def slajder_3(value):
    vr=abs(float(value))*div_fact.get()* eslfrom_3.get()/100
    e3num.set(round(vr,2))
    calculate()

def slajder_4(value):
    vr=abs(float(value))*div_fact.get()* eslfrom_4.get()/100
    e4num.set(round(vr,2))
    calculate()

def slajder_5(value):
    vr=abs(float(value))*div_fact.get()* eslfrom_5.get()/100
    e5num.set(round(vr,2))
    calculate()

def slajder_6(value):
    vr=abs(float(value))*div_fact.get()* eslfrom_6.get()/100
    e6num.set(round(vr,2))
    calculate()

def slajder_7(value):
    vr=abs(float(value))*div_fact.get()* eslfrom_7.get()/100
    e7num.set(round(vr,2))
    calculate()

def slajder_8(value):
    vr=abs(float(value))*div_fact.get()* eslfrom_8.get()/100
    e8num.set(round(vr,2))
    calculate()

eslfrom_1=DoubleVar(); eslto_1=DoubleVar(); eslfrom_1.set(razpon)
eslfrom_2=DoubleVar(); eslto_2=DoubleVar(); eslfrom_2.set(razpon)
eslfrom_3=DoubleVar(); eslto_3=DoubleVar(); eslfrom_3.set(razpon)
eslfrom_4=DoubleVar(); eslto_4=DoubleVar(); eslfrom_4.set(razpon)
eslfrom_5=DoubleVar(); eslto_5=DoubleVar(); eslfrom_5.set(razpon)
eslfrom_6=DoubleVar(); eslto_6=DoubleVar(); eslfrom_6.set(razpon)
eslfrom_7=DoubleVar(); eslto_7=DoubleVar(); eslfrom_7.set(razpon)
eslfrom_8=DoubleVar(); eslto_8=DoubleVar(); eslfrom_8.set(razpon)

e1_from=Entry(root, textvariable=eslfrom_1, font=('12'), width=4, justify='center'); e1_to=Entry(root, textvariable=eslto_1, font=('12'), width=4, justify='center')
e2_from=Entry(root, textvariable=eslfrom_2, font=('12'), width=4, justify='center'); e2_to=Entry(root, textvariable=eslto_2, font=('12'), width=4, justify='center')
e3_from=Entry(root, textvariable=eslfrom_3, font=('12'), width=4, justify='center'); e3_to=Entry(root, textvariable=eslto_3, font=('12'), width=4, justify='center')
e4_from=Entry(root, textvariable=eslfrom_4, font=('12'), width=4, justify='center'); e4_to=Entry(root, textvariable=eslto_4, font=('12'), width=4, justify='center')
e5_from=Entry(root, textvariable=eslfrom_5, font=('12'), width=4, justify='center'); e5_to=Entry(root, textvariable=eslto_5, font=('12'), width=4, justify='center')
e6_from=Entry(root, textvariable=eslfrom_6, font=('12'), width=4, justify='center'); e6_to=Entry(root, textvariable=eslto_6, font=('12'), width=4, justify='center')
e7_from=Entry(root, textvariable=eslfrom_7, font=('12'), width=4, justify='center'); e7_to=Entry(root, textvariable=eslto_7, font=('12'), width=4, justify='center')
e8_from=Entry(root, textvariable=eslfrom_8, font=('12'), width=4, justify='center'); e8_to=Entry(root, textvariable=eslto_8, font=('12'), width=4, justify='center')

e1_from.grid(row=5,column=2); e1_to.grid(row=18,column=2)
e2_from.grid(row=5,column=3); e2_to.grid(row=18,column=3)
e3_from.grid(row=5,column=4); e3_to.grid(row=18,column=4)
e4_from.grid(row=5,column=5); e4_to.grid(row=18,column=5)
e5_from.grid(row=5,column=6); e5_to.grid(row=18,column=6)
e6_from.grid(row=5,column=7); e6_to.grid(row=18,column=7)
e7_from.grid(row=5,column=8); e7_to.grid(row=18,column=8)
e8_from.grid(row=5,column=9); e8_to.grid(row=18,column=9)

def slajder_range():
    slider_1.config(from_=-eslfrom_1.get(), to=-eslto_1.get() )
    slider_2.config(from_=-eslfrom_2.get(), to=-eslto_2.get() )
    slider_3.config(from_=-eslfrom_3.get(), to=-eslto_3.get() )
    slider_4.config(from_=-eslfrom_4.get(), to=-eslto_4.get() )
    slider_5.config(from_=-eslfrom_5.get(), to=-eslto_5.get() )
    slider_6.config(from_=-eslfrom_6.get(), to=-eslto_6.get() )
    slider_7.config(from_=-eslfrom_7.get(), to=-eslto_7.get() )
    slider_8.config(from_=-eslfrom_8.get(), to=-eslto_8.get() )

Button(root, text='Set range', command=slajder_range, font=('44')).grid(row=20, column=6, columnspan=2)   

dolzina_slajderjev = 200

slider_1=Scale(root, from_=-eslfrom_1.get(), to=eslto_1.get(), orient=VERTICAL, length=dolzina_slajderjev, showvalue=0,command=slajder_1)
slider_2=Scale(root, from_=-eslfrom_2.get(), to=eslto_2.get(), orient=VERTICAL, length=dolzina_slajderjev, showvalue=0,command=slajder_2)
slider_3=Scale(root, from_=-eslfrom_3.get(), to=eslto_3.get(), orient=VERTICAL, length=dolzina_slajderjev, showvalue=0,command=slajder_3)
slider_4=Scale(root, from_=-eslfrom_4.get(), to=eslto_4.get(), orient=VERTICAL, length=dolzina_slajderjev, showvalue=0,command=slajder_4)
slider_5=Scale(root, from_=-eslfrom_5.get(), to=eslto_5.get(), orient=VERTICAL, length=dolzina_slajderjev, showvalue=0,command=slajder_5)
slider_6=Scale(root, from_=-eslfrom_6.get(), to=eslto_6.get(), orient=VERTICAL, length=dolzina_slajderjev, showvalue=0,command=slajder_6)
slider_7=Scale(root, from_=-eslfrom_7.get(), to=eslto_7.get(), orient=VERTICAL, length=dolzina_slajderjev, showvalue=0,command=slajder_7)
slider_8=Scale(root, from_=-eslfrom_8.get(), to=eslto_8.get(), orient=VERTICAL, length=dolzina_slajderjev, showvalue=0,command=slajder_8)

slider_1.grid(row=6, column=2, rowspan=10)
slider_2.grid(row=6, column=3, rowspan=10)
slider_3.grid(row=6, column=4, rowspan=10)
slider_4.grid(row=6, column=5, rowspan=10)
slider_5.grid(row=6, column=6, rowspan=10)
slider_6.grid(row=6, column=7, rowspan=10)
slider_7.grid(row=6, column=8, rowspan=10)
slider_8.grid(row=6, column=9, rowspan=10)

# CLICK SAVE
from tkinter.filedialog import asksaveasfile as savefile

def click_file_save():
        global cf
        cf= savefile(mode='a', defaultextension=".csv")    # shranjevanje rezultatov "na klik" v obliki .csv datoteke
        if cf is None:
                return
imt=0
def click_save():
        global imt
        imt=imt+1
        if imt==1:
                click_file_save()
        else:
                pass
        if bool(cf)==True:
                try:
                        with open(cf.name,cf.mode) as cg:
                                        cg.write(res_formula.get()+'\n')
                                        cg.write('   delta(%)   ,   chi(%)   ,   VEC   ,   dHmix(kJ/mol)   ,   dSmix(J/mol)   ,   Tm(deg.C)   ,   omega'+'\n')
                                        cg.write(res_delta.get()+','+res_kappa.get()+','+
                                                 res_VEC.get()+','+res_dHmix.get()+','+res_dSmix.get()+','+res_Tm.get()+','+
                                                 res_omega.get()+'\n')
                                        cg.write(op_delta.get()+','+op_kappa.get()+','+op_VEC.get()+','+op_dHmix.get()+',,,'+
                                                 op_omega.get()+'\n')
                                        cg.write('\n')           
                except AttributeError:
                        pass


# AUTO SAVE

def auto_file_save():
        global kol
        kol=savefile(mode='a', defaultextension=".csv", initialfile='auto_save_data.csv')    # shranjevanje rezultatov iteracij v obliki .csv datoteke
        if kol is None:
                return

display_num_points=StringVar()
num_points=-1
def auto_save():
        global num_points
        try:
            try:
                num_points+=1
                with open(kol.name,kol.mode) as han:
                    han.write(res_formula.get()+','+'calculated'+','+res_delta.get()+','+res_dHmix.get()+','+res_kappa.get()+','+res_omega.get()+','+res_VEC.get()+'\n')     
            except AttributeError:
                pass
        except NameError:
                auto_file_save()
                with open(kol.name,kol.mode) as huh:
                    huh.write('  formula  ,  Structure   ,   delta(%)   ,   Hmix kJ/mol   ,  chi(%)  ,   omega  ,  VEC  '+'\n')        

Button(root,text='Save',command=click_save, fg='green',font=('44')).grid(row=20, column=7, columnspan=2, sticky=E)

# Izbira atomskih ali masnih %

def izbira_funkcije():
        fn.get()
        if fn.get()==0:
            formjula.config(fg='#1e90ff')
            calculate()
            return 0     # - atomski %
        elif fn.get()==1:
            formjula.config(fg='#808080')
            calculate()
            return 1     # - masni %
        elif fn.get()==2:
            formjula.config(fg='red')
            mass_to_atomic(k)

Radiobutton(root, text='Atomic %',indicatoron = 0, variable=fn, value=0,fg='#1e90ff',font=('Arial','14'),command=izbira_funkcije, width=8).grid(row=1,column=12, sticky=E)
Radiobutton(root, text='Mass %',indicatoron = 0, variable=fn, value=1,fg='#808080',font=('Arial','14'),command=izbira_funkcije, width=8).grid(row=1,column=13, sticky=W)
Radiobutton(root, text=' Mass ==> At. % ',indicatoron = 0, variable=fn, value=2,fg='red',font=('Arial','14'),command=izbira_funkcije, width=14).grid(row=1,column=14, sticky=W)

def mass_to_atomic(pretvori):        #  PRETVORNIK   masnih % (npr. rezultati z EDS-a)v atomske %
    terna={}
    vrsta=[]
    for i in k:
        vrsta.append(k[i]/molske_data[i])
        iks=sum(vrsta)
    for j,i in zip(k,vrsta):
        terna[j]=round((100*i/iks),1)
    ama=list(terna.keys())
    bma=list(terna.values())
    formula_mass_to_atomic = " ".join([i+'-'+str(round(j,1))for i,j in zip(ama,bma)])
    res_formula.set('         '+formula_mass_to_atomic+' (at.%)')

Label(root, text='').grid(row=21, column=1)


# AVTOMATSKI ITERATOR oziroma ITERATOR MACHINE :)

""" Dodatek originalnemu programu, napisal Andraz Kocjan, v avgustu 2018, na Inštitutu za Materiale in tehnologije.
    Ta platforma omogoca iteracijo kemijskih sestav zlitin glede na izbrane elemente in intervale atomskih delezev ter
    izracun termodinamskih parametrov, ki jih nato (graficno)primerjamo z obstojecimi podatki iz literature. Oboje nam
    nudi hitrejso in lazjo pot pri iskanju novih visokoentropijskih zlitin."""

frame_bg='#b3b3cc'
frstc='#5c5c8a'
frame = Frame(root, bg=frame_bg, width=1500)
frame.grid(row=23, column=1, columnspan=30, sticky=W)
frame.config(borderwidth=3, relief=SUNKEN)

Label(frame, text='', bg=frame_bg).grid(column=1, row=23)
Label(frame, text=u"\u25cf"+' ITERATOR MACHINE '+u"\u25cf", fg='#ffffcc', font=('Arial','14','bold'), bg=frame_bg).grid(row=24, column=1, columnspan=1, sticky=NW)
Label(frame, text='', bg=frame_bg).grid(column=1, row=25)

# 1. ELEMENTS to VARY entries
Label(frame, text='Addition   ', fg='black', font=('Arial','11'), bg=frame_bg).grid(row=25, column=1, columnspan=1, sticky=E)
Label(frame, text='Elements to vary:  ', fg='black', font=('Arial','12', 'bold'), bg=frame_bg).grid(row=26, column=1, columnspan=1, sticky=E)
ai1var=StringVar(); ai2var=StringVar(); ai3var=StringVar(); ai4var=StringVar(); ai5var=StringVar(); ai6var=StringVar(); ai7var=StringVar(); ai8var=StringVar();
ai9var=StringVar(); ai10var=StringVar(); ai11var=StringVar(); ai12var=StringVar(); ai13var=StringVar(); ai14var=StringVar()
ai_bg='black'; ai_fg='#00ff00'; ai_cur='#00ff00'
ai1=Entry(frame, textvariable=ai1var, font=font_e1, width=e_width, bg=ai_bg, fg=ai_fg, justify='center'); ai1.grid(column=2, row=26, sticky=E); ai1.config(insertbackground=ai_cur)
ai2=Entry(frame, textvariable=ai2var, font=font_e1, width=e_width, bg=ai_bg, fg=ai_fg, justify='center'); ai2.grid(column=3, row=26, sticky=E); ai2.config(insertbackground=ai_cur)
ai3=Entry(frame, textvariable=ai3var, font=font_e1, width=e_width, bg=ai_bg, fg=ai_fg, justify='center'); ai3.grid(column=4, row=26, sticky=E); ai3.config(insertbackground=ai_cur)
ai4=Entry(frame, textvariable=ai4var, font=font_e1, width=e_width, bg=ai_bg, fg=ai_fg, justify='center'); ai4.grid(column=5, row=26, sticky=E); ai4.config(insertbackground=ai_cur)
ai5=Entry(frame, textvariable=ai5var, font=font_e1, width=e_width, bg=ai_bg, fg=ai_fg, justify='center'); ai5.grid(column=6, row=26, sticky=E); ai5.config(insertbackground=ai_cur)
ai6=Entry(frame, textvariable=ai6var, font=font_e1, width=e_width, bg=ai_bg, fg=ai_fg, justify='center'); ai6.grid(column=7, row=26, sticky=E); ai6.config(insertbackground=ai_cur)
ai7=Entry(frame, textvariable=ai7var, font=font_e1, width=e_width, bg=ai_bg, fg=ai_fg, justify='center'); ai7.grid(column=8, row=26, sticky=E); ai7.config(insertbackground=ai_cur)
ai8=Entry(frame, textvariable=ai8var, font=font_e1, width=e_width, bg=ai_bg, fg=ai_fg, justify='center'); ai8.grid(column=9, row=26, sticky=E); ai8.config(insertbackground=ai_cur)
ai9=Entry(frame, textvariable=ai9var, font=font_e1, width=e_width, bg=ai_bg, fg=ai_fg, justify='center'); ai9.grid(column=10, row=26, sticky=E); ai9.config(insertbackground=ai_cur)
ai10=Entry(frame, textvariable=ai10var, font=font_e1, width=e_width, bg=ai_bg, fg=ai_fg, justify='center'); ai10.grid(column=11, row=26, sticky=E); ai10.config(insertbackground=ai_cur)
ai11=Entry(frame, textvariable=ai11var, font=font_e1, width=e_width, bg=ai_bg, fg=ai_fg, justify='center'); ai11.grid(column=12, row=26, sticky=E); ai11.config(insertbackground=ai_cur)
ai12=Entry(frame, textvariable=ai12var, font=font_e1, width=e_width, bg=ai_bg, fg=ai_fg, justify='center'); ai12.grid(column=13, row=26, sticky=E); ai12.config(insertbackground=ai_cur)
ai13=Entry(frame, textvariable=ai13var, font=font_e1, width=e_width, bg=ai_bg, fg=ai_fg, justify='center'); ai13.grid(column=14, row=26, sticky=E); ai13.config(insertbackground=ai_cur)
ai14=Entry(frame, textvariable=ai14var, font=font_e1, width=e_width, bg=ai_bg, fg=ai_fg, justify='center'); ai14.grid(column=15, row=26, sticky=E); ai14.config(insertbackground=ai_cur)

# 2. START_RATIO entries
Label(frame, text='INITIAL at.%:  ', fg='black', font=('Arial','11'), bg=frame_bg).grid(row=27, column=1, columnspan=1, sticky=E)
im1var=IntVar(); im2var=IntVar(); im3var=IntVar(); im4var=IntVar(); im5var=IntVar(); im6var=IntVar(); im7var=IntVar(); im8var=IntVar(); im9var=IntVar();
im10var=IntVar(); im11var=IntVar(); im12var=IntVar(); im13var=IntVar(); im14var=IntVar()
im_bg='black'; im_fg='#99ccff'; im_cur='#99ccff'
im1=Entry(frame, textvariable=im1var, font=font_e1, width=e_width, bg=frstc, fg=im_fg, justify='center'); im1.grid(column=2, row=27, sticky=E); im1.config(insertbackground=im_cur)
im2=Entry(frame, textvariable=im2var, font=font_e1, width=e_width, bg=im_bg, fg=im_fg, justify='center'); im2.grid(column=3, row=27, sticky=E); im2.config(insertbackground=im_cur)
im3=Entry(frame, textvariable=im3var, font=font_e1, width=e_width, bg=im_bg, fg=im_fg, justify='center'); im3.grid(column=4, row=27, sticky=E); im3.config(insertbackground=im_cur)
im4=Entry(frame, textvariable=im4var, font=font_e1, width=e_width, bg=im_bg, fg=im_fg, justify='center'); im4.grid(column=5, row=27, sticky=E); im4.config(insertbackground=im_cur)
im5=Entry(frame, textvariable=im5var, font=font_e1, width=e_width, bg=im_bg, fg=im_fg, justify='center'); im5.grid(column=6, row=27, sticky=E); im5.config(insertbackground=im_cur)
im6=Entry(frame, textvariable=im6var, font=font_e1, width=e_width, bg=im_bg, fg=im_fg, justify='center'); im6.grid(column=7, row=27, sticky=E); im6.config(insertbackground=im_cur)
im7=Entry(frame, textvariable=im7var, font=font_e1, width=e_width, bg=im_bg, fg=im_fg, justify='center'); im7.grid(column=8, row=27, sticky=E); im7.config(insertbackground=im_cur)
im8=Entry(frame, textvariable=im8var, font=font_e1, width=e_width, bg=im_bg, fg=im_fg, justify='center'); im8.grid(column=9, row=27, sticky=E); im8.config(insertbackground=im_cur)
im9=Entry(frame, textvariable=im9var, font=font_e1, width=e_width, bg=im_bg, fg=im_fg, justify='center'); im9.grid(column=10, row=27, sticky=E); im9.config(insertbackground=im_cur)
im10=Entry(frame, textvariable=im10var, font=font_e1, width=e_width, bg=im_bg, fg=im_fg, justify='center'); im10.grid(column=11, row=27, sticky=E); im10.config(insertbackground=im_cur)
im11=Entry(frame, textvariable=im11var, font=font_e1, width=e_width, bg=im_bg, fg=im_fg, justify='center'); im11.grid(column=12, row=27, sticky=E); im11.config(insertbackground=im_cur)
im12=Entry(frame, textvariable=im12var, font=font_e1, width=e_width, bg=im_bg, fg=im_fg, justify='center'); im12.grid(column=13, row=27, sticky=E); im12.config(insertbackground=im_cur)
im13=Entry(frame, textvariable=im13var, font=font_e1, width=e_width, bg=im_bg, fg=im_fg, justify='center'); im13.grid(column=14, row=27, sticky=E); im13.config(insertbackground=im_cur)
im14=Entry(frame, textvariable=im14var, font=font_e1, width=e_width, bg=im_bg, fg=im_fg, justify='center'); im14.grid(column=15, row=27, sticky=E); im14.config(insertbackground=im_cur)

# 3. STOP_RATIO entries
Label(frame, text='FINAL at.%:  ', fg='black', font=('Arial','11'), bg=frame_bg).grid(row=28, column=1, columnspan=1, sticky=E)
rs1var=IntVar(); rs2var=IntVar(); rs3var=IntVar(); rs4var=IntVar(); rs5var=IntVar(); rs6var=IntVar(); rs7var=IntVar(); rs8var=IntVar(); rs9var=IntVar();
rs10var=IntVar(); rs11var=IntVar(); rs12var=IntVar(); rs13var=IntVar(); rs14var=IntVar()
rs_bg_end='black'; rs_fg_end='#cc99ff'; rs_cur_end='#cc99ff'
rs1=Entry(frame, textvariable=rs1var, font=font_e1, width=e_width, bg=frstc, fg=rs_fg_end, justify='center'); rs1.grid(column=2, row=28, sticky=E); rs1.config(insertbackground=rs_cur_end)
rs2=Entry(frame, textvariable=rs2var, font=font_e1, width=e_width, bg=rs_bg_end, fg=rs_fg_end, justify='center'); rs2.grid(column=3, row=28, sticky=E); rs2.config(insertbackground=rs_cur_end)
rs3=Entry(frame, textvariable=rs3var, font=font_e1, width=e_width, bg=rs_bg_end, fg=rs_fg_end, justify='center'); rs3.grid(column=4, row=28, sticky=E); rs3.config(insertbackground=rs_cur_end)
rs4=Entry(frame, textvariable=rs4var, font=font_e1, width=e_width, bg=rs_bg_end, fg=rs_fg_end, justify='center'); rs4.grid(column=5, row=28, sticky=E); rs4.config(insertbackground=rs_cur_end)
rs5=Entry(frame, textvariable=rs5var, font=font_e1, width=e_width, bg=rs_bg_end, fg=rs_fg_end, justify='center'); rs5.grid(column=6, row=28, sticky=E); rs5.config(insertbackground=rs_cur_end)
rs6=Entry(frame, textvariable=rs6var, font=font_e1, width=e_width, bg=rs_bg_end, fg=rs_fg_end, justify='center'); rs6.grid(column=7, row=28, sticky=E); rs6.config(insertbackground=rs_cur_end)
rs7=Entry(frame, textvariable=rs7var, font=font_e1, width=e_width, bg=rs_bg_end, fg=rs_fg_end, justify='center'); rs7.grid(column=8, row=28, sticky=E); rs7.config(insertbackground=rs_cur_end)
rs8=Entry(frame, textvariable=rs8var, font=font_e1, width=e_width, bg=rs_bg_end, fg=rs_fg_end, justify='center'); rs8.grid(column=9, row=28, sticky=E); rs8.config(insertbackground=rs_cur_end)
rs9=Entry(frame, textvariable=rs9var, font=font_e1, width=e_width, bg=rs_bg_end, fg=rs_fg_end, justify='center'); rs9.grid(column=10, row=28, sticky=E); rs9.config(insertbackground=rs_cur_end)
rs10=Entry(frame, textvariable=rs10var, font=font_e1, width=e_width, bg=rs_bg_end, fg=rs_fg_end, justify='center'); rs10.grid(column=11, row=28, sticky=E); rs10.config(insertbackground=rs_cur_end)
rs11=Entry(frame, textvariable=rs11var, font=font_e1, width=e_width, bg=rs_bg_end, fg=rs_fg_end, justify='center'); rs11.grid(column=12, row=28, sticky=E); rs11.config(insertbackground=rs_cur_end)
rs12=Entry(frame, textvariable=rs12var, font=font_e1, width=e_width, bg=rs_bg_end, fg=rs_fg_end, justify='center'); rs12.grid(column=13, row=28, sticky=E); rs12.config(insertbackground=rs_cur_end)
rs13=Entry(frame, textvariable=rs13var, font=font_e1, width=e_width, bg=rs_bg_end, fg=rs_fg_end, justify='center'); rs13.grid(column=14, row=28, sticky=E); rs13.config(insertbackground=rs_cur_end)
rs14=Entry(frame, textvariable=rs14var, font=font_e1, width=e_width, bg=rs_bg_end, fg=rs_fg_end, justify='center'); rs14.grid(column=15, row=28, sticky=E); rs14.config(insertbackground=rs_cur_end)

# 4. STEP_SIZE entries
Label(frame, text='STEP:  ', fg='black', font=('Arial','11'), bg=frame_bg).grid(row=29, column=1, columnspan=1, sticky=E)
ss1var=DoubleVar(); ss2var=DoubleVar(); ss3var=DoubleVar(); ss4var=DoubleVar(); ss5var=DoubleVar(); ss6var=DoubleVar(); ss7var=DoubleVar(); ss8var=DoubleVar(); ss9var=DoubleVar();
ss10var=DoubleVar(); ss11var=DoubleVar(); ss12var=DoubleVar(); ss13var=DoubleVar(); ss14var=DoubleVar()
ss_bg='black'; ss_fg='#ffff99'; ss_cur='#ffff99'
ss1=Entry(frame, textvariable=ss1var, font=font_e1, width=e_width, bg=frstc, fg=ss_fg, justify='center'); ss1.grid(column=2, row=29, sticky=E); ss1.config(insertbackground=ss_cur)
ss2=Entry(frame, textvariable=ss2var, font=font_e1, width=e_width, bg=ss_bg, fg=ss_fg, justify='center'); ss2.grid(column=3, row=29, sticky=E); ss2.config(insertbackground=ss_cur)
ss3=Entry(frame, textvariable=ss3var, font=font_e1, width=e_width, bg=ss_bg, fg=ss_fg, justify='center'); ss3.grid(column=4, row=29, sticky=E); ss3.config(insertbackground=ss_cur)
ss4=Entry(frame, textvariable=ss4var, font=font_e1, width=e_width, bg=ss_bg, fg=ss_fg, justify='center'); ss4.grid(column=5, row=29, sticky=E); ss4.config(insertbackground=ss_cur)
ss5=Entry(frame, textvariable=ss5var, font=font_e1, width=e_width, bg=ss_bg, fg=ss_fg, justify='center'); ss5.grid(column=6, row=29, sticky=E); ss5.config(insertbackground=ss_cur)
ss6=Entry(frame, textvariable=ss6var, font=font_e1, width=e_width, bg=ss_bg, fg=ss_fg, justify='center'); ss6.grid(column=7, row=29, sticky=E); ss6.config(insertbackground=ss_cur)
ss7=Entry(frame, textvariable=ss7var, font=font_e1, width=e_width, bg=ss_bg, fg=ss_fg, justify='center'); ss7.grid(column=8, row=29, sticky=E); ss7.config(insertbackground=ss_cur)
ss8=Entry(frame, textvariable=ss8var, font=font_e1, width=e_width, bg=ss_bg, fg=ss_fg, justify='center'); ss8.grid(column=9, row=29, sticky=E); ss8.config(insertbackground=ss_cur)
ss9=Entry(frame, textvariable=ss9var, font=font_e1, width=e_width, bg=ss_bg, fg=ss_fg, justify='center'); ss9.grid(column=10, row=29, sticky=E); ss9.config(insertbackground=ss_cur)
ss10=Entry(frame, textvariable=ss10var, font=font_e1, width=e_width, bg=ss_bg, fg=ss_fg, justify='center'); ss10.grid(column=11, row=29, sticky=E); ss10.config(insertbackground=ss_cur)
ss11=Entry(frame, textvariable=ss11var, font=font_e1, width=e_width, bg=ss_bg, fg=ss_fg, justify='center'); ss11.grid(column=12, row=29, sticky=E); ss11.config(insertbackground=ss_cur)
ss12=Entry(frame, textvariable=ss12var, font=font_e1, width=e_width, bg=ss_bg, fg=ss_fg, justify='center'); ss12.grid(column=13, row=29, sticky=E); ss12.config(insertbackground=ss_cur)
ss13=Entry(frame, textvariable=ss13var, font=font_e1, width=e_width, bg=ss_bg, fg=ss_fg, justify='center'); ss13.grid(column=14, row=29, sticky=E); ss13.config(insertbackground=ss_cur)
ss14=Entry(frame, textvariable=ss14var, font=font_e1, width=e_width, bg=ss_bg, fg=ss_fg, justify='center'); ss14.grid(column=15, row=29, sticky=E); ss14.config(insertbackground=ss_cur)

# 5. DECIMALS entries
Label(frame, text='COMBINATIONS at.%:  ', fg='black', font=('Arial','11'), bg=frame_bg).grid(row=30, column=1, columnspan=1, sticky=E)
dd1var=DoubleVar(); dd2var=DoubleVar(); dd3var=DoubleVar(); dd4var=DoubleVar(); dd5var=DoubleVar(); dd6var=DoubleVar(); dd7var=DoubleVar(); dd8var=DoubleVar(); dd9var=DoubleVar();
dd10var=DoubleVar(); dd11var=DoubleVar(); dd12var=DoubleVar(); dd13var=DoubleVar(); dd14var=DoubleVar()
dd_bg='black'; dd_fg='#ffffff'; dd_cur='#ffffff'
dd1=Entry(frame, textvariable=dd1var, font=font_e1, width=e_width, bg=frstc, fg=dd_fg, justify='center'); dd1.grid(column=2, row=30, sticky=E); dd1.config(insertbackground=dd_cur)
dd2=Entry(frame, textvariable=dd2var, font=font_e1, width=e_width, bg=dd_bg, fg=dd_fg, justify='center'); dd2.grid(column=3, row=30, sticky=E); dd2.config(insertbackground=dd_cur)
dd3=Entry(frame, textvariable=dd3var, font=font_e1, width=e_width, bg=dd_bg, fg=dd_fg, justify='center'); dd3.grid(column=4, row=30, sticky=E); dd3.config(insertbackground=dd_cur)
dd4=Entry(frame, textvariable=dd4var, font=font_e1, width=e_width, bg=dd_bg, fg=dd_fg, justify='center'); dd4.grid(column=5, row=30, sticky=E); dd4.config(insertbackground=dd_cur)
dd5=Entry(frame, textvariable=dd5var, font=font_e1, width=e_width, bg=dd_bg, fg=dd_fg, justify='center'); dd5.grid(column=6, row=30, sticky=E); dd5.config(insertbackground=dd_cur)
dd6=Entry(frame, textvariable=dd6var, font=font_e1, width=e_width, bg=dd_bg, fg=dd_fg, justify='center'); dd6.grid(column=7, row=30, sticky=E); dd6.config(insertbackground=dd_cur)
dd7=Entry(frame, textvariable=dd7var, font=font_e1, width=e_width, bg=dd_bg, fg=dd_fg, justify='center'); dd7.grid(column=8, row=30, sticky=E); dd7.config(insertbackground=dd_cur)
dd8=Entry(frame, textvariable=dd8var, font=font_e1, width=e_width, bg=dd_bg, fg=dd_fg, justify='center'); dd8.grid(column=9, row=30, sticky=E); dd8.config(insertbackground=dd_cur)
dd9=Entry(frame, textvariable=dd9var, font=font_e1, width=e_width, bg=dd_bg, fg=dd_fg, justify='center'); dd9.grid(column=10, row=30, sticky=E); dd9.config(insertbackground=dd_cur)
dd10=Entry(frame, textvariable=dd10var, font=font_e1, width=e_width, bg=dd_bg, fg=dd_fg, justify='center'); dd10.grid(column=11, row=30, sticky=E); dd10.config(insertbackground=dd_cur)
dd11=Entry(frame, textvariable=dd11var, font=font_e1, width=e_width, bg=dd_bg, fg=dd_fg, justify='center'); dd11.grid(column=12, row=30, sticky=E); dd11.config(insertbackground=dd_cur)
dd12=Entry(frame, textvariable=dd12var, font=font_e1, width=e_width, bg=dd_bg, fg=dd_fg, justify='center'); dd12.grid(column=13, row=30, sticky=E); dd12.config(insertbackground=dd_cur)
dd13=Entry(frame, textvariable=dd13var, font=font_e1, width=e_width, bg=dd_bg, fg=dd_fg, justify='center'); dd13.grid(column=14, row=30, sticky=E); dd13.config(insertbackground=dd_cur)
dd14=Entry(frame, textvariable=dd14var, font=font_e1, width=e_width, bg=dd_bg, fg=dd_fg, justify='center'); dd14.grid(column=15, row=30, sticky=E); dd14.config(insertbackground=dd_cur)

# OSNOVA fiksno
e1var.set('Fe'); e2var.set('Cr'); e3var.set('Ni')
e1num.set(1); e2num.set(1); e3num.set(1)

# SPREMENLJIVKE fiksno
ai1var.set('Ti')#;  ai2var.set('Zr')        # elementi
im1var.set(0)#;  im2var.set(0)            # zacetni atomski %
rs1var.set(25)#; im22var.set(25)      # koncni atomski %
ss1var.set(5.0)#; ss2var.set(5)             # korak (at.%)
dd1var.set(1.0)#; dd2var.set(1)             # scaling

nps= Label(frame, textvariable=nosilec_podatkov_status, width=23, bg='#47476b', fg='#ff3300', font=('Calibry Body', '10', 'bold'))
nps.grid(row=25, column=13, columnspan=3)

# ADDITION of elements for Substitution iteration - using radiobuttons
font_add_buttons=('Arial 11')
mikromodes={}
for i in range(1,12):
    mikromodes["add_mode{}".format(i)]=StringVar()
MODES=Namespace(**mikromodes)
for i in range(1,12):
    exec("MODES.add_mode{}.set('x')".format(i))
urnebes={}
longines=[i for i in range(1,12)]
modes_dict=MODES.__dict__
for n,m in zip(longines, modes_dict):
    urnebes[n]=modes_dict[m]
rbtn={}
for i in range(1,12):
    rbtn["button{}x_add".format(i)] = Radiobutton (frame, text='x', font=font_add_buttons, selectcolor='#00cc99', variable=urnebes[i], value='x', indicatoron=0, width=2)
    rbtn["button{}y_add".format(i)] = Radiobutton (frame, text='y', font=font_add_buttons, selectcolor='#cc6699', variable=urnebes[i], value='y', indicatoron=0, width=2)
gma=Namespace(**rbtn)

for i in range(1,12):
    exec("gma.button{0}x_add.grid(column={1},row=25,sticky=W,padx=3)".format(i,i+1))
    exec("gma.button{0}y_add.grid(column={1},row=25,sticky=E,padx=3)".format(i,i+1))

Master_formula=StringVar()
Label(root, text='').grid(row=16, rowspan=3, column=10, columnspan=6)
mstf=Label(root, textvariable=Master_formula, fg='#336699', font=('Arial 24'))
mstf.grid(row=17, rowspan=5, column=10, columnspan=6)

def pick_addition():
    X_add={}; Y_add={}
    selected_adds=[urnebes[i].get()for i in urnebes]
    itd=get_iterator_data()
    try:
        x_range=itd[0][1:4]
        y_range=itd[1][1:4]
    except IndexError:
        y_range=[]
    for i in range(len(itd)):
        try:
            if selected_adds[i]=='x':
                X_add[itd[i][0]]=itd[i][4]
            else:
                Y_add[itd[i][0]]=itd[i][4]
        except IndexError:
            pass
    return X_add, Y_add, x_range, y_range

# Implicite formula of substitution iteration
def implicitna_formula(sub, add):
    all_subs=normalize_100(get_data())
    sub_x={}; sub_y={}; sub_xy={}; sub_same={}
    add_x={}; add_y={}; add_xy={}; add_same={}
    for i in sub[0]:#............................................... (AaBb)1-x  ==> SUB_X
        sub_x[i]=all_subs[i]
    if bool(init_x.get())==False:
        x_init=float(round(sum(list(sub_x.values())),1))
        if (x_init).is_integer():
            x_init=int(x_init)
    else:
        x_init=float(init_x.get())
    sub_xnorm_raw=normalize(sub_x); sub_xnorm={}
    for x in sub_xnorm_raw:
        sub_xnorm[x]=round(sub_xnorm_raw[x]/min(list(sub_xnorm_raw.values())),1)
    n=[]
    for i,j in zip(list(sub_xnorm.keys()),list(sub_xnorm.values())):
        if all(k ==list(sub_xnorm.values())[0] for k in list(sub_xnorm.values())) or j==1.0:
            n.append(i)
        else:
            if (j).is_integer():
                n.append(i+str(int(j)))
            else:
                n.append(i+str(round(j,2)))
    if bool("".join(n)):
        sub_xstring='('+"".join(n)+')'+str(x_init)+ '-x'
    else:
        sub_xstring=''

    for i in sub[1]:#............................................... (CcDd)1-y  ==> SUB_Y
        sub_y[i]=all_subs[i]
    if bool(init_y.get())==False:
        y_init=float(round(sum(list(sub_y.values())),1))
        if (y_init).is_integer():
            y_init=int(y_init)
    else:
        y_init=float(init_y.get())
    sub_ynorm_raw=normalize(sub_y); sub_ynorm={}
    for y in sub_ynorm_raw:
        sub_ynorm[y]=round(sub_ynorm_raw[y]/min(list(sub_ynorm_raw.values())),1)
    n=[]
    for i,j in zip(list(sub_ynorm.keys()),list(sub_ynorm.values())):
        if all(k ==list(sub_ynorm.values())[0] for k in list(sub_ynorm.values())) or j==1.0:
            n.append(i)
        else:
            if (j).is_integer():
                n.append(i+str(int(j)))
            else:
                n.append(i+str(round(j,2)))
    if bool("".join(n)):
        sub_ystring='('+"".join(n)+')'+str(y_init)+ '-m'
    else:
        sub_ystring=''
        
    same=100-x_init-y_init   #.................................................... GgHh ==> SUB_SAME
    sub3_norm=normalize(sub[3])
    for i in sub[3]:
        sub_same[i]=sub3_norm[i]*same
    n=[]    
    for i,j in zip(list(sub_same.keys()),list(sub_same.values())):
        if (j).is_integer():
            n.append(i+str(int(j)))
        else:
            n.append(i+str(round(j,1)))
    sub_samestring="".join(n)

    add_x_raw=normalize(add[0])#.......................................... (IiJj)x ==> ADD_X
    n=[]
    for x in add_x_raw:
        add_x[x]=round(add_x_raw[x]/min(list(add_x_raw.values())),1)
    for i,j in zip(list(add_x.keys()),list(add_x.values())):
        if all(k ==list(add_x.values())[0] for k in list(add_x.values())) or j==1.0:
            n.append(i)
        else:
            if (j).is_integer():
                n.append(i+str(int(j)))
            else:
                n.append(i+str(round(j,2)))
    if bool("".join(n)):
        add_xstring='('+"".join(n)+')'+ 'x'
    else:
        add_xstring=''

    add_y_raw=normalize(add[1])#................................ (KkLl)x ==> ADD_Y
    n=[]
    for y in add_y_raw:
        add_y[y]=round(add_y_raw[y]/min(list(add_y_raw.values())),1)
    for i,j in zip(list(add_y.keys()),list(add_y.values())):
        if all(k ==list(add_y.values())[0] for k in list(add_y.values())) or j==1.0:
            n.append(i)
        else:
            if (j).is_integer():
                n.append(i+str(int(j)))
            else:
                n.append(i+str(round(j,2)))
    if bool("".join(n)):
        add_ystring='('+"".join(n)+')'+ 'm'
    else:
        add_ystring=''

    if switch.get()=='small_iterator_data':
        prikaz_formule_raw=sub_xstring+' '+sub_samestring+' Mx'
        prikaz_formule_subscripts=prikaz_formule_raw.translate(Subscripts)
        Master_formula.set(prikaz_formule_subscripts)
    elif switch.get()=='substitution_iterator_data':
        prikaz_formule_raw=sub_xstring+' '+add_xstring+' '+sub_ystring+' '+add_ystring+' '+sub_samestring
        prikaz_formule_subscripts=prikaz_formule_raw.translate(Subscripts)
        Master_formula.set(prikaz_formule_subscripts)
        
    return x_init, sub_xstring, sub_samestring, add_xstring, y_init, sub_ystring, add_ystring
    

# Grabs the data for all iteration mechanisms
def get_iterator_data():
    big_hub=[]
    elements_entries=[ai1var.get(), ai2var.get(), ai3var.get(), ai4var.get(), ai5var.get(), ai6var.get(), ai7var.get(), ai8var.get(), ai9var.get(), ai10var.get(), ai11var.get(), ai12var.get(), ai13var.get(), ai14var.get()]
    initial_entries=[im1var.get(), im2var.get(), im3var.get(), im4var.get(), im5var.get(), im6var.get(), im7var.get(), im8var.get(), im9var.get(), im10var.get(), im11var.get(), im12var.get(), im13var.get(), im14var.get()]
    final_entries=[rs1var.get(), rs2var.get(), rs3var.get(), rs4var.get(), rs5var.get(), rs6var.get(), rs7var.get(), rs8var.get(), rs9var.get(), rs10var.get(), rs11var.get(), rs12var.get(), rs13var.get(), rs14var.get()]
    step_entries=[ss1var.get(),ss2var.get(),ss3var.get(),ss4var.get(),ss5var.get(),ss6var.get(),ss7var.get(),ss8var.get(),ss9var.get(),ss10var.get(),ss11var.get(),ss12var.get(),ss13var.get(),ss14var.get()]
    decimals_entries=[dd1var.get(), dd2var.get(), dd3var.get(), dd4var.get(), dd5var.get(), dd6var.get(), dd7var.get(), dd8var.get(), dd9var.get(), dd10var.get(), dd11var.get(), dd12var.get(), dd13var.get(), dd14var.get()]
    for a,b,c,d,e in zip(elements_entries, initial_entries, final_entries, step_entries, decimals_entries):
        big_hub.append([a,b,c,d,e])
    for u,t in enumerate(big_hub):
        if bool(t[0])== True:
            cr = t[0][0]
            if cr.islower():
                rk=cr.replace(cr,cr.upper())
                try:
                    big_hub[u][0]=rk+big_hub[u][0][1]
                except IndexError:
                    big_hub[u][0]=rk
    big_hub=[iks for iks in big_hub if bool(iks[0])==True]
    return big_hub


""" ALGORITMI za ITERATOR MACHINE """
import itertools
# LIMITS entries ; vnosi kriterijev shranjevanja:
t=-1
limit_mode=StringVar()
def settings():     
    global kriteriji, t
    settings = Toplevel()
    settings.geometry('560x400+600+200')
    settings.title('Nastavitve kriterijev shranjevanja')
    frame_bg='#669999'
    settings.config(bg=frame_bg)
    entry_width= 7; font_entry = ('Arial', '22'); bg_entry='black' ; fg_entry='#ffffff'
    font_label = ('Arial', '16', 'bold'); fg_label='black'
    font_title=('Arial','14')
    
    Label(settings, text='      ', bg=frame_bg, height=5, width=10).grid(row=0, column=0)
    Label(settings, text='      ', bg=frame_bg, height=2).grid(row=6, column=0)
    Label(settings, text='MIN', bg=frame_bg, font=font_title).grid(row=5, column=6)
    Label(settings, text='parameter', bg=frame_bg, font=font_title, width=10).grid(row=5, column=7)
    Label(settings, text='MAX', bg=frame_bg, font=font_title).grid(row=5, column=8)

    lim_del_min=DoubleVar(); lim_del_max=DoubleVar(); lim_H_min=DoubleVar(); lim_H_max=DoubleVar()
    lim_S_min=DoubleVar(); lim_S_max=DoubleVar()

    Label(settings, text=' <  '+u"\u0394"+'Smix  <    ', fg=fg_label, font=font_label, bg=frame_bg).grid(column=7, row=9)    
    # delta
    Label(settings, text=' <  '+u"\u0394".lower()+'  <    ', fg=fg_label, font=font_label, bg=frame_bg).grid(column=7, row=7)
    lim1=Entry(settings, textvariable=lim_del_min, font=font_entry, width=entry_width,bg=bg_entry, fg=fg_entry, justify='center'); lim1.grid(column=6, row=7, sticky=E); lim1.config(insertbackground=fg_entry)
    lim2=Entry(settings, textvariable=lim_del_max, font=font_entry, width=entry_width,bg=bg_entry, fg=fg_entry, justify='center'); lim2.grid(column=8, row=7, sticky=E); lim2.config(insertbackground=fg_entry)
    # Hmix
    Label(settings, text=' <  '+u"\u0394"+'Hmix  <    ', fg=fg_label, font=font_label, bg=frame_bg).grid(column=7, row=8)
    lim3=Entry(settings, textvariable=lim_H_min, font=font_entry, width=entry_width,bg=bg_entry, fg=fg_entry, justify='center'); lim3.grid(column=6, row=8, sticky=E); lim3.config(insertbackground=fg_entry)
    lim4=Entry(settings, textvariable=lim_H_max, font=font_entry, width=entry_width,bg=bg_entry, fg=fg_entry, justify='center'); lim4.grid(column=8, row=8, sticky=E); lim4.config(insertbackground=fg_entry)
    # Smix    
    Label(settings, text=' <  '+u"\u0394"+'Smix  <    ', fg=fg_label, font=font_label, bg=frame_bg).grid(column=7, row=9)
    lim5=Entry(settings, textvariable=lim_S_min, font=font_entry, width=entry_width,bg=bg_entry, fg=fg_entry, justify='center'); lim5.grid(column=6, row=9, sticky=E); lim5.config(insertbackground=fg_entry)
    lim6=Entry(settings, textvariable=lim_S_max, font=font_entry, width=entry_width,bg=bg_entry, fg=fg_entry, justify='center'); lim6.grid(column=8, row=9, sticky=E); lim6.config(insertbackground=fg_entry)    

    # Default values of limits
    def_del_min=0 ;  def_del_max=6.6 ;  def_H_min=-12 ;  def_H_max=-5 ;  def_S_min=12.5 ;  def_S_max=100
    # Unlimited values
    unlim_del_min=0 ;  unlim_del_max=100 ;  unlim_H_min=-100 ;  unlim_H_max=100 ;  unlim_S_min=0 ;  unlim_S_max=100
    t +=1
    def set_default():
        global kriteriji
        lim_del_min.set(def_del_min); lim_del_max.set(def_del_max); lim_H_min.set(def_H_min); lim_H_max.set(def_H_max); lim_S_min.set(def_S_min); lim_S_max.set(def_S_max)
        kriteriji=[lim_del_min.get(), lim_del_max.get(), lim_H_min.get(), lim_H_max.get(), lim_S_min.get(), lim_S_max.get()]
        print(kriteriji)
        limit_mode.set('default limits')
        settings.withdraw()
        return kriteriji

    def set_unlimited():
        global kriteriji
        lim_del_min.set(unlim_del_min); lim_del_max.set(unlim_del_max); lim_H_min.set(unlim_H_min); lim_H_max.set(unlim_H_max); lim_S_min.set(unlim_S_min); lim_S_max.set(unlim_S_max)
        kriteriji=[lim_del_min.get(), lim_del_max.get(), lim_H_min.get(), lim_H_max.get(), lim_S_min.get(), lim_S_max.get()]
        print(kriteriji)
        limit_mode.set('unlimited')
        settings.withdraw()
        return kriteriji

    if t==0:
        kriteriji=set_unlimited()
    else:
        lim_del_min.set(kriteriji[0]); lim_del_max.set(kriteriji[1]); lim_H_min.set(kriteriji[2]); lim_H_max.set(kriteriji[3]); lim_S_min.set(kriteriji[4]); lim_S_max.set(kriteriji[5])
        
    def set_kriteriji():
        global kriteriji
        kriteriji=[lim_del_min.get(), lim_del_max.get(), lim_H_min.get(), lim_H_max.get(), lim_S_min.get(), lim_S_max.get()]
        print(kriteriji)
        limit_mode.set('custom limits')
        settings.withdraw()
        return kriteriji

    Label(settings, text='', fg=fg_label, font=font_label, bg=frame_bg).grid(column=7, row=13)
    Label(settings, text='', fg=fg_label, font=font_label, bg=frame_bg).grid(column=7, row=14)
    Button(settings, text='         Set        ', justify=CENTER, width=20, command=set_kriteriji).grid(row=15, column=7, columnspan=3, sticky=W)
    Button(settings, text='         Default        ', justify=CENTER, width=20, command=set_default).grid(row=15, column=6, columnspan=3, sticky=W)
    Button(settings, text='         No limits        ', justify=CENTER, width=20, command=set_unlimited).grid(row=15, column=8, columnspan=3, sticky=W)
    return kriteriji

kriteriji=settings()

def normalize_100(insert):
    podatki=normalize(insert)
    bbb=list(podatki.values())
    bbb=[u*100 for u in bbb]
    for i,j in zip(podatki,bbb):
            podatki[i]=j
    return podatki

def stevc():
    global stevec
    stevec = 0
    return stevec

stevc()

def vary_conc(vnos, element, start_ratio, stop_ratio, step_size, criteria):
    global stevec
    rezultat=[]
    sub=get_subs(); add=pick_addition(); input_x=sub[0]; input_same=sub[3]; x_init=implicitna_formula(sub,add)[0]
    old_x=normalize_100(input_x)
    for i in old_x:
        old_x[i]=old_x[i]*x_init/100
    try:
        same=normalize(input_same)
        same_init=100-x_init   # začetni delež elementov, ki se ne zamenjujejo
        for i in same:
                same[i]=same[i]*same_init
        vnos={**old_x, **same}
    except NameError:
        vnos=old_x.copy()
    vnos_0=vnos.copy()
    old_x=normalize(old_x)
    #vnos=normalize_100(vnos)
    #vnos_0=vnos.copy()
    od = int(start_ratio/step_size)
    do = int(stop_ratio/step_size)
    for i in range(od, do+1):
        stevec+=1
        atomic_perc={}
        vnos[element]=i*step_size
        for x in old_x:
            #vnos[x]=round(vnos_0[x]*(1-(i*step_size)/100), 1)
            vnos[x]= old_x[x]*(x_init-i*step_size)
        try:
            mega=splosni_izracun(vnos, packing.get())
        except UnboundLocalError:
                pass
        # kriterij shranjevanja v Python list
        if mega[4] > criteria[4] and mega[3] > criteria[2] and mega[3] < criteria [3] and mega[0] > criteria[0] and mega[0] < criteria[1]:
            
            if choice_mode.get()=='basic' and switch.get()!='big_iterator_data':
                rezultat.append(('M='+element, mega[-2][element],) + mega)
            elif choice_mode.get()=='modified' or switch.get()=='big_iterator_data':
                rezultat.append(mega)
        # tu doloci kriterije za shranjevanje
        #auto_save()
    return rezultat


""" UNDER CONSTRUCTION """

def vary_nova(vnos, element, start_ratio, stop_ratio, step_size):
    sub=get_subs(); add=pick_addition(); input_x=sub[0]; input_same=sub[3]; x_init=implicitna_formula(sub,add)[0]
    old_x=normalize_100(input_x)
    for i in old_x:
        old_x[i]=old_x[i]*x_init/100
    try:
        same=normalize(input_same)
        same_init=100-x_init   # začetni delež elementov, ki se ne zamenjujejo
        for i in same:
                same[i]=same[i]*same_init
        vnos={**old_x, **same}
    except NameError:
        vnos=old_x.copy()
    vnos_0=vnos.copy()
    old_x=normalize(old_x)

    od = int(start_ratio/step_size)
    do = int(stop_ratio/step_size)
    for i in range(od, do+1):
        vnos[element]=i*step_size
        for x in old_x:
            vnos[x]= old_x[x]*(x_init-i*step_size)
        try:
            mega=splosni_izracun(vnos, packing.get())
        except UnboundLocalError:
                pass
        print(mega[12])

        


""" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - """



def vary_conc_BIG(vnos, element, start_ratio, stop_ratio, step_size, criteria):
    global stevec
    rezultat=[]
    vnos=normalize_100(vnos)
    vnos_0=vnos.copy()
    od = int(start_ratio/step_size)
    do = int(stop_ratio/step_size)
    for i in range(od, do+1):
        stevec+=1
        atomic_perc={}
        vnos[element]=i*step_size
        for x in vnos_0:
            vnos[x]=round(vnos_0[x]*(1-(i*step_size)/100), 1)
        try:
            mega=splosni_izracun(vnos, packing.get())
        except UnboundLocalError:
                pass
        # kriterij shranjevanja v Python list
        if mega[4] > criteria[4] and mega[3] > criteria[2] and mega[3] < criteria [3] and mega[0] > criteria[0] and mega[0] < criteria[1]:
                rezultat.append(mega)
        # tu doloci kriterije za shranjevanje
        #auto_save()
    return rezultat


def vary_subs(criteria):
    global stevec, rezultat_substitution
    stevc()
    sub=get_subs(); add=pick_addition()
    input_x=sub[0]; input_y=sub[1]; input_xy=sub[2]; input_same=sub[3]
    rezultat_substitution=[]
    
    x_init=implicitna_formula(sub,add)[0]
    y_init=implicitna_formula(sub,add)[4]
    
    old_x=normalize_100(input_x)
    old_y=normalize_100(input_y)
    for i in old_x:
        old_x[i]=old_x[i]*x_init/100
    for i in old_y:
        old_y[i]=old_y[i]*y_init/100
    try:
        same=normalize(input_same)
        same_init=100-x_init-y_init   # začetni delež elementov, ki se ne zamenjujejo
        for i in same:
                same[i]=same[i]*same_init
        vnos={**old_x, **old_y, **same}
    except NameError:
        vnos=old_x.copy()

    # ++++++++ enako sledi za input_xy..
  
    x_od = int(add[2][0]/add[2][2]); x_do = int(add[2][1]/add[2][2])
    old_x=normalize(old_x); add_norm=normalize(add[0])
    if bool(sub[1])== True or y_init!=0 or bool(add[1])==True:
        old_y=normalize(old_y); add_norm_y=normalize(add[1])
        y_od = int(add[3][0]/add[3][2]); y_do = int(add[3][1]/add[3][2])
        for j in range(y_od, y_do+1):   #  addition of  Y
            for y in add[1]:
                vnos[y]= j*add[3][2]*add_norm_y[y]
            for y in old_y:
                vnos[y]= old_y[y]*(y_init-j*add[3][2])
            rez_sub=[]		    
            for i in range(x_od, x_do+1):   #  addition of  X
                stevec+=1
                for x in add[0]:
                    vnos[x]= i*add[2][2]*add_norm[x]
                for x in old_x:
                    vnos[x]= old_x[x]*(x_init-i*add[2][2])
                try:
                    mega=splosni_izracun(vnos, packing.get())
                except UnboundLocalError:
                        pass
                # kriterij shranjevanja v Python list
                if mega[4] > criteria[4] and mega[3] > criteria[2] and mega[3] < criteria [3] and mega[0] > criteria[0] and mega[0] < criteria[1]:
                    if choice_mode.get()=='basic':
                        rez_sub.append(('m='+str(j*add[3][2]),i*add[2][2],) + mega)
                    else:
                        rez_sub.append(mega)
                # tu doloci kriterije za shranjevanje
                #auto_save()
            rezultat_substitution.append(rez_sub)
    else:
        rez_sub=[]		    
        for i in range(x_od, x_do+1):   #  addition of  X
            stevec+=1
            for x in add[0]:
                vnos[x]= i*add[2][2]*add_norm[x]
            for x in old_x:
                vnos[x]= old_x[x]*(x_init-i*add[2][2])
            try:
                mega=splosni_izracun(vnos, packing.get())
            except UnboundLocalError:
                    pass
            # kriterij shranjevanja v Python list
            if mega[4] > criteria[4] and mega[3] > criteria[2] and mega[3] < criteria [3] and mega[0] > criteria[0] and mega[0] < criteria[1]:
                if choice_mode.get()=='basic':
                    rez_sub.append((implicitna_formula(sub,add)[1],i*add[2][2],) + mega)
                else:
                    rez_sub.append(mega)
            # tu doloci kriterije za shranjevanje
            #auto_save()
        rezultat_substitution.append(rez_sub)
        
    lift=[r for r in rezultat_substitution if r!=[]]
    ma=list(itertools.chain(*lift))
    counter=len(ma)
    procent = round(100*counter/stevec, 1)
    counter_display.set('Substitution Iterator report:   "Got '+str(counter)+' alloys! That is '+str(procent)+' % out of '+str(stevec)+' iterations!"')
    print(' ===>>>  Got ',counter,' alloys! That is ',procent,' % out of ',stevec,' iterations! <<<===')
    switch.set('substitution_iterator_data')
    nosilec_podatkov.append((rezultat_substitution, 'rezultat_substitution'))
    stevilo_zlitin(counter)
    nps.config(fg='lime')
    implicitna_formula(sub,add)
    return rezultat_substitution

switch = StringVar(); counter_display = StringVar()

def BIG_iterator(criteria):
    global BIG_data
    switch.set('big_iterator_data')
    stevc(); BIG_data=[]; w=get_data(); rotor=get_iterator_data()
    vnos=normalize_100(w)
    if len(rotor)==1:
        first_el=rotor[0][0]
        mm=rotor[0]
        frka=vary_conc_BIG(vnos, first_el, mm[1], mm[2], mm[3], criteria)
        BIG_data.append(frka)

    if len(rotor)==2:
        first_el=rotor[0][0]; second_el=rotor[1][0]
        vnos[first_el]=0
        g0=w.keys()
        vnos_0=vnos.copy()
        od1 = int(rotor[0][1]/rotor[0][3]); do1 = int(rotor[0][2]/rotor[0][3])
        for i in range(od1, do1+1):
            vnos[first_el]=i*rotor[0][3]
            for el0 in g0:
                vnos[el0]= vnos_0[el0]*(1-(i*rotor[0][3])/100) #""" To vrstico popravi tako kot je v osnovnem vary_conc iteratorju. """
            mm=rotor[1]
            frka=vary_conc_BIG(vnos, second_el, mm[1], mm[2], mm[3], criteria)
            BIG_data.append(frka)

    if len(rotor)==3:
        first_el=rotor[0][0]; second_el=rotor[1][0]; third_el=rotor[2][0]
        vnos[first_el]=0; vnos[second_el]=0
        g0=w.keys()
        g1=list(g0)+[first_el]
        vnos_0=vnos.copy()
        od1 = int(rotor[0][1]/rotor[0][3]); do1 = int(rotor[0][2]/rotor[0][3])
        od2 = int(rotor[1][1]/rotor[1][3]); do2 = int(rotor[1][2]/rotor[1][3])
        for i in range(od1, do1+1):
            vnos[first_el]=i*rotor[0][3]
            for el0 in g0:
                vnos[el0]= vnos_0[el0]*(1-(i*rotor[0][3])/100)
            vnos_1=vnos.copy()
            for j in range(od2, do2+1):
                vnos[second_el]=j*rotor[1][3]
                for el1 in g1:
                    vnos[el1]= vnos_1[el1]*(1-(j*rotor[1][3])/100)
                mm=rotor[2]
                frka=vary_conc_BIG(vnos, third_el, mm[1], mm[2], mm[3], criteria)
                BIG_data.append(frka)        

    if len(rotor)==4:
        first_el=rotor[0][0]; second_el=rotor[1][0]; third_el=rotor[2][0]; fourth_el=rotor[3][0]
        vnos[first_el]=0; vnos[second_el]=0; vnos[third_el]=0
        g0=w.keys()
        g1=list(g0)+[first_el]
        g2=g1+[second_el]
        vnos_0=vnos.copy()
        od1 = int(rotor[0][1]/rotor[0][3]); od2 = int(rotor[1][1]/rotor[1][3]); od3 = int(rotor[2][1]/rotor[2][3])
        do1 = int(rotor[0][2]/rotor[0][3]); do2 = int(rotor[1][2]/rotor[1][3]); do3 = int(rotor[2][2]/rotor[2][3])
 
        for i in range(od1, do1+1):
            vnos[first_el]=i*rotor[0][3]
            for el0 in g0:
                vnos[el0]= vnos_0[el0]*(1-(i*rotor[0][3])/100)
            vnos_1=vnos.copy()
            for j in range(od2, do2+1):
                vnos[second_el]=j*rotor[1][3]
                for el1 in g1:
                    vnos[el1]= vnos_1[el1]*(1-(j*rotor[1][3])/100)
                vnos_2=vnos.copy()
                for k in range(od3, do3+1):
                    vnos[third_el]=k*rotor[2][3]
                    for el2 in g2:
                        vnos[el2]=vnos_2[el2]*(1-(k*rotor[2][3])/100)
                    vnos[fourth_el]=0
                    mm=rotor[3]
                    frka=vary_conc_BIG(vnos, fourth_el, mm[1], mm[2], mm[3], criteria)
                    #frka.append(second_el+'= '+str(i*rotor[1][3]))
                    BIG_data.append(frka)

    if len(rotor)==5:
        first_el=rotor[0][0]; second_el=rotor[1][0]; third_el=rotor[2][0]; fourth_el=rotor[3][0]; fifth_el=rotor[4][0]
        vnos[first_el]=0; vnos[second_el]=0; vnos[third_el]=0; vnos[fourth_el]=0
        g0=w.keys()
        g1=list(g0)+[first_el]
        g2=g1+[second_el]
        g3=g2+[third_el]
        vnos_0=vnos.copy()
        od1 = int(rotor[0][1]/rotor[0][3]); od2 = int(rotor[1][1]/rotor[1][3]); od3 = int(rotor[2][1]/rotor[2][3]); od4 = int(rotor[3][1]/rotor[3][3])
        do1 = int(rotor[0][2]/rotor[0][3]); do2 = int(rotor[1][2]/rotor[1][3]); do3 = int(rotor[2][2]/rotor[2][3]); do4 = int(rotor[3][2]/rotor[3][3])
 
        for i in range(od1, do1+1):  # 1st element
            vnos[first_el]=i*rotor[0][3]
            for el0 in g0:
                vnos[el0]= vnos_0[el0]*(1-(i*rotor[0][3])/100)
            vnos_1=vnos.copy()
            
            for j in range(od2, do2+1):  # 2nd element
                vnos[second_el]=j*rotor[1][3]
                for el1 in g1:
                    vnos[el1]= vnos_1[el1]*(1-(j*rotor[1][3])/100)
                vnos_2=vnos.copy()
                
                for k in range(od3, do3+1):  # 3rd element
                    vnos[third_el]=k*rotor[2][3]
                    for el2 in g2:
                        vnos[el2]=vnos_2[el2]*(1-(k*rotor[2][3])/100)
                    vnos_3=vnos.copy()
                    
                    for l in range(od4, do4+1):  # 4th element
                        vnos[fourth_el]=l*rotor[3][3]
                        for el3 in g3:
                            vnos[el3]=vnos_3[el3]*(1-(l*rotor[3][3])/100)
                            
                        vnos[fifth_el]=0
                        mm=rotor[4]
                        frka=vary_conc_BIG(vnos, fifth_el, mm[1], mm[2], mm[3], criteria)# 5th element
                        #frka.append(second_el+'= '+str(i*rotor[1][3]))
                        BIG_data.append(frka)

    if len(rotor)==6:
        first_el=rotor[0][0]; second_el=rotor[1][0]; third_el=rotor[2][0]; fourth_el=rotor[3][0]; fifth_el=rotor[4][0]; sixth_el=rotor[5][0]
        vnos[first_el]=0; vnos[second_el]=0; vnos[third_el]=0; vnos[fourth_el]=0; vnos[fifth_el]=0
        g0=w.keys()
        g1=list(g0)+[first_el]
        g2=g1+[second_el]
        g3=g2+[third_el]
        g4=g3+[fourth_el]
        vnos_0=vnos.copy()
        od1 = int(rotor[0][1]/rotor[0][3]); od2 = int(rotor[1][1]/rotor[1][3]); od3 = int(rotor[2][1]/rotor[2][3]); od4 = int(rotor[3][1]/rotor[3][3]); od5 = int(rotor[4][1]/rotor[4][3])
        do1 = int(rotor[0][2]/rotor[0][3]); do2 = int(rotor[1][2]/rotor[1][3]); do3 = int(rotor[2][2]/rotor[2][3]); do4 = int(rotor[3][2]/rotor[3][3]); do5 = int(rotor[4][2]/rotor[4][3])
 
        for i in range(od1, do1+1):  # 1st element
            vnos[first_el]=i*rotor[0][3]
            for el0 in g0:
                vnos[el0]= vnos_0[el0]*(1-(i*rotor[0][3])/100)
            vnos_1=vnos.copy()
            
            for j in range(od2, do2+1):  # 2nd element
                vnos[second_el]=j*rotor[1][3]
                for el1 in g1:
                    vnos[el1]= vnos_1[el1]*(1-(j*rotor[1][3])/100)
                vnos_2=vnos.copy()
                
                for k in range(od3, do3+1):  # 3rd element
                    vnos[third_el]=k*rotor[2][3]
                    for el2 in g2:
                        vnos[el2]=vnos_2[el2]*(1-(k*rotor[2][3])/100)
                    vnos_3=vnos.copy()
                    
                    for l in range(od4, do4+1):  # 4th element
                        vnos[fourth_el]=l*rotor[3][3]
                        for el3 in g3:
                            vnos[el3]=vnos_3[el3]*(1-(l*rotor[3][3])/100)
                            vnos_4=vnos.copy()

                            for m in range(od5, do5+1):  # 5th element
                                vnos[fifth_el]=m*rotor[4][3]
                                for el4 in g4:
                                    vnos[el4]=vnos_4[el4]*(1-(m*rotor[4][3])/100)
                            
                        vnos[sixth_el]=0
                        mm=rotor[5]
                        frka=vary_conc_BIG(vnos, sixth_el, mm[1], mm[2], mm[3], criteria)# 6th element
                        #frka.append(second_el+'= '+str(i*rotor[1][3]))
                        BIG_data.append(frka)
            
    lift=[r for r in BIG_data if r!=[]]
    ma=list(itertools.chain(*lift))
    counter=len(ma)
    procent = round(100*counter/stevec, 1)
    counter_display.set('BIG Iterator report:   "Got '+str(counter)+' alloys! That is '+str(procent)+' % out of '+str(stevec)+' iterations!"')
    print(' ===>>>  Got ',counter,' alloys! That is ',procent,' % out of ',stevec,' iterations! <<<===')
    nosilec_podatkov.append((BIG_data, 'BIG_data'))
    stevilo_zlitin(counter)
    nps.config(fg='lime')
    return BIG_data

def small_iterator(criteria):
    global results
    stevc()
    results=[]  
    w=get_data()
    small_elems=get_iterator_data()
    for sit in small_elems:
        w[sit[0]]=0
        try:
            frka=vary_conc(w,sit[0], sit[1], sit[2], sit[3], criteria)
            #frka.append('variing element = '+sit[0])
            results.append(frka)
        except UnboundLocalError:
            pass
        del w[sit[0]]
        
    lift=[r for r in results if r!=[]]
    ma=list(itertools.chain(*lift))
    counter=len(ma)
    procent = round(100*counter/stevec, 1)
    counter_display.set('Small Iterator report:   "Got '+str(counter)+' alloys! That is '+str(procent)+' % out of '+str(stevec)+' iterations!"')
    print(' ===>>>  Got ',counter,' alloys! That is ',procent,' % out of ',stevec,' iterations! <<<===')
    switch.set('small_iterator_data')
    nosilec_podatkov.append((results, 'results'))
    stevilo_zlitin(counter)
    nps.config(fg='lime'); plotSmall_btn.config(state="normal", fg='black')
    return results


from pari import COMBINATIONS as fff

no_el=IntVar()
no_el.set(5)

def kombinatorika(criteria):
    global results_comb
    results_comb=[]
    v=get_data()
    nu=get_iterator_data()
    q={}
    for p in nu:
        q[p[0]]=1
    w=fff(q)
    for o in w[no_el.get()-2]:
        e={}
        m=re.findall('[A-Z][a-z]*',o)
        for j in m:
            for c in nu:
                if c[0]==j:
                    e[j]=c[-1]            
        fu={**v,**e}
        try:
            odg=splosni_izracun(fu, packing.get())
        except UnboundLocalError:
            pass
        # kriterij shranjevanja v Python list
        if odg[4] > criteria[4] and odg[3] > criteria[2] and odg[3] < criteria [3] and odg[0] > criteria[0] and odg[0] < criteria[1]:
            results_comb.append(odg)
        # tu doloci kriterije shranjevanja v csv.file
        #auto_save()
    moznih=len(w[no_el.get()-2])
    shranjenih=len(results_comb)
    print(moznih,'  moznih zlitin. Po filtriranju shranjenih ',shranjenih,' sestav.')
    procent = round(100*shranjenih/moznih, 1)
    counter_display.set('Combinations report:   "Got '+str(shranjenih)+' alloys! That is '+str(procent)+' % out of '+str(moznih)+' iterations!"')
    #switch.set('combination_iterator_data')
    nosilec_podatkov.append((results_comb, 'results_comb'))
    stevilo_zlitin(shranjenih)
    nps.config(fg='lime')
    Master_formula.set('')
    return results_comb
        
aaa=Button(frame, text='BIG Iter', command=lambda: BIG_iterator(kriteriji)); aaa.grid(row=24, column=3, columnspan=2, sticky=W)#aaa.grid(row=24, column=2, columnspan=3)
bbb=Button(frame, text='Small Iter', command=lambda: small_iterator(kriteriji)); bbb.grid(row=24, column=3, columnspan=3, sticky=E)
ccc=Button(frame, text='Substitution', command=lambda:vary_subs(kriteriji)); ccc.grid(row=24, column=6, columnspan=2, sticky=W)
ddd=Button(frame, text='Combinations', command=lambda: kombinatorika(kriteriji)); ddd.grid(row=24, column=6, columnspan=3, sticky=E)
Label(frame, text='No.- of elements', font=('Arial','10'), fg='#0066ff', bg=frame_bg).grid(column=10, row=24, columnspan=2, sticky=W)
num_komb=Entry(frame, textvariable=no_el, font=font_e1, width=e_width, bg='#003366', fg='#00ffff', justify='center')
num_komb.grid(column=9, row=24, columnspan=1, sticky=E); num_komb.config(insertbackground='#00ffff')
Button(frame, text='Save As', command=auto_file_save).grid(row=24, column=12, columnspan=2, sticky=W)

aaa.config(width=12); bbb.config(width=12); ccc.config(width=12); ddd.config(width=12)

Label(frame, text='', bg=frame_bg).grid(row=26, column=16)
all_elements=['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Al', 'Si', 'Zr', 'Nb', 'Mo', 'Ta','W', 'Hf','Sn','Y']     # A list of 14 transition metal candidates for HEAs

def input_elements():
    wow=get_data()
    in_el=[]
    for i in all_elements:
        if i not in wow.keys():
            in_el.append(i)
    for u,rr in enumerate(in_el):
        exec("ai{0}var.set('{1}')".format(u+1,rr))
    input_values()
    return

def input_values():
    p=get_iterator_data()
    for i in range(len(p)-1):
        exec("im{0}var.set('{1}')".format(i+2,p[0][1]))
        exec("rs{0}var.set('{1}')".format(i+2,p[0][2]))
        exec("ss{0}var.set('{1}')".format(i+2,p[0][3]))
        exec("dd{0}var.set('{1}')".format(i+2,p[0][4]))
    return
    
def clear_all():
    for sss in range(1,15):
        exec("ai{}var.set('')".format(sss))
Button(frame, text='Insert Elements', width=15, command=input_elements).grid(row=24, column=17)
Button(frame, text='Input Values', width=15,command=input_values).grid(row=25, column=17)
Button(frame, text='Clear Elements', width=15,command=clear_all).grid(row=26, column=17)
Label(frame, text=' ', bg=frame_bg).grid(row=26, column=19)

x_varka=StringVar(); y_varka=StringVar()
Label(frame, text='X - axis', bg=frame_bg, font=('Arial','13')).grid(row=28, column=17, sticky=S); Label(frame, text='Y - axis', bg=frame_bg, font=('Arial','13')).grid(row=28, column=18, sticky=S)
OPTIONS=["Atomic %", '    '+u"\u0394".lower()+'    ', '    '+u"\u03C7".lower()+'    ', "   VEC   ", '  '+u"\u0394"+"Hmix  ", '  '+u"\u0394"+"Smix  ", "    Tm   ", '    '+u"\u03a9"+'   ', "  e/a   "]
x_pick = OptionMenu(frame, x_varka, *OPTIONS); x_pick.grid(row=29, column=17, rowspan=2, sticky=N)
y_pick = OptionMenu(frame, y_varka, *OPTIONS); y_pick.grid(row=29, column=18, rowspan=2, sticky=N)
x_pick.config(width=7, font=font_e1, justify='center', fg='#3399ff', activebackground='#ffff00')
y_pick.config(width=7, font=font_e1, justify='center', fg='#ff66cc', activebackground='black', activeforeground='#ffff00')
x_varka.set(OPTIONS[1]); y_varka.set(OPTIONS[4])

def nearest(array, e,f):
    x=[]; y=[]; dist=[]; gt=(1,2,3)
    if choice_mode.get()=='modified':
        merged_array = list(itertools.chain(*array))
        x = [n[0] for n in merged_array]   
        y = [m[1] for m in merged_array]  
    elif choice_mode.get()=='basic':
        if switch.get()=='small_iterator_data' or switch.get()=='substitution_iterator_data':
            d={}
            for j,i in enumerate(OPTIONS):
                    d[i]=j+1
        elif switch.get()=='big_iterator_data':
            d={}
            for j,i in enumerate(OPTIONS):
                if i !='Atomic %':
                    d[i]=j-1
        merged_array = list(itertools.chain(*array))
        for h,k in enumerate(merged_array):
            if isinstance(k,type(gt))==False:
                del merged_array[h]
        x+= [n[d[x_varka.get()]] for n in merged_array]
        y+= [m[d[y_varka.get()]] for m in merged_array]
    tx=abs(max(x)- min(x))
    ty=abs(max(y)- min(y)) 
    r=ty/tx
    y=[i/r for i in y]
    z=[(i,j)for i,j in zip(x,y)]
    for i in z:
        d=math.sqrt((i[0]-e)**2 + (i[1]-f/r)**2)
        dist.append(d)
    n=dist.index(min(dist))
    return merged_array[n]

def nearest_comb(array, e,f):
    d={}; dist=[]
    for j,i in enumerate(OPTIONS):
        if i !='Atomic %':
            d[i]=j-1
    if choice_mode.get()=='basic':
        x=[n[d[x_varka.get()]] for n in array]
        y=[m[d[y_varka.get()]] for m in array]
    else:
        x=[n[0] for n in array]
        y=[m[1] for m in array]
    tx=abs(max(x)- min(x))
    ty=abs(max(y)- min(y))
    r=ty/tx
    y=[i/r for i in y]
    z=[(i,j)for i,j in zip(x,y)]
    for i in z:
        d=math.sqrt((i[0]-e)**2 + (i[1]-f/r)**2)
        dist.append(d)
    n=dist.index(min(dist))
    return array[n]

def onclick_comb(event):
    global good_job_comb, pozicija_detajlov
    afna= (event.xdata, event.ydata)
    if choice_mode.get()=='basic':
        good_job_comb = nearest_comb(results_comb, afna[0], afna[1])
    else:
        good_job_comb = nearest_comb(advanced, afna[0], afna[1])
    #print(good_job_comb)
    try:
        pozicija_detajlov=details.winfo_geometry()
        details.destroy()
        get_details()
    except  NameError:
        pass

import operator, statistics

bz=-1
def get_details():
    global details, bz
    bz+=1
    on = switch.get()
    if on == 'combination_iterator_data':
        data = good_job_comb[:-2]
        a = parske_entalpije(good_job_comb[-2])
        struktura=good_job_comb[-1]
        
    elif on == 'small_iterator_data':
        data = good_job_small[:-2]
        a = parske_entalpije(good_job_small[-2])
        struktura=good_job_small[-1]

    elif on == 'substitution_iterator_data':
        data = good_job_substitution[:-2]
        a = parske_entalpije(good_job_substitution[-2])
        struktura=good_job_substitution[-1]
        
    elif on == 'big_iterator_data':
        data = good_job_big[:-2]
        a = parske_entalpije(good_job_big[-2])
        struktura=good_job_big[-1]

    elif on=='total_plot' or on=='corr_plot':
        data = good_job_lit[:-2] 
        a = parske_entalpije(good_job_lit[-2])
        struktura=good_job_lit[-1]  # ==> podatki o kristalni strukturi (FCC, BCC, BMG,...)oz. načinu izračuna
  
    background='#ffcc99'
    details = Toplevel()
    details.title('Details about selected alloy')
    details.configure(background=background)
    font_1 = ('Arial', '13', 'bold'); fg_1='#4040bf' ; fg_2='#ff3399' ; font_2 =('Arial','14'); font_0=('Arial', '11')
    Label(details, text='      ', bg=background, height=2, width=2).grid(row=0, column=0, columnspan=10)
    Label(details, text='      ', bg=background).grid(row=6, column=0, columnspan=10)
    Label(details,text='Detailed data for alloy (in at.%): ', fg='#4040bf', bg=background, font=('Arial','18')).grid(row=5, column=0, columnspan=10)
    Label(details,text='    '+data[-1], fg='#4040bf', bg=background, font=('Arial','22')).grid(row=7, column=0, sticky=E, columnspan=10)
    Label(details, text='                                                           ', bg=background).grid(row=8, column=1, columnspan=1)
    Label(details, text='     Pair mixing enthalpies:  ', bg=background, font=font_2, fg=fg_1, width=20).grid(row=9, column=0, columnspan=4, sticky=NW)
    Label(details, text='     ', bg=background, font=font_2, fg=fg_1).grid(row=9, column=4, columnspan=1, sticky=NE) 
    b=sorted(a.items(), key=operator.itemgetter(1))
    for j,i in enumerate(b):
        Label(details, text=i[0]+' :  ', bg=background, fg=fg_1, font=font_0).grid(row=10+j, column=0, sticky=E)
        Label(details, text=str(i[1])+'   kJ/ mol', bg=background, fg='#9933ff', font=font_0).grid(row=10+j, column=1, sticky=W)

    c=[i[1] for i in a.items()]
    std = round(statistics.stdev(c), 3)
    Label(details, text='             -------------------------', bg=background, fg=fg_1, font=font_1).grid(row=11+j, column=0, columnspan=2, sticky=W)
    Label(details, text='STD'+' :  ', bg=background, fg=fg_1, font=font_1).grid(row=12+j, column=0, sticky=E)
    Label(details, text=str(std)+'  %', bg=background, fg='#9933ff', font=font_1).grid(row=12+j, column=1, sticky=W)
    def draw_pair_enthalpies():
        plt.ion()
        fig,ax = plt.subplots()
        plt.bar(a.keys(), a.values())
        ax.set_ylim(-10, 5)
        plt.title(data[-1]+'   '+text_strukture[struktura])
    Label(details,text='', bg=background).grid(row=13+j, column=1)
    Button(details, text='Plot pair enthalpies' , command=draw_pair_enthalpies).grid(row=14+j, column=0, sticky=E)
    font_str=("Arial", "30", "bold")
    barve_strukture={'FCC':'aqua','BCC':'lime','FCC+BCC':'#b366ff','BMG':'#cc0000','IM':'grey','calculated':'magenta'}
    font_strukture={'FCC':font_str,'BCC':font_str,'FCC+BCC':font_str,'BMG':font_str,'IM':font_str,'calculated':("Arial", "22", "italic")}
    text_strukture={'FCC':'FCC','BCC':'BCC','FCC+BCC':'FCC+BCC','BMG':'BMG','IM':'IM','calculated':'unknown'}
    Label(details, text='Structure:   ', bg=background, font=font_2, fg=fg_1).grid(row=9, column=3, columnspan=3, sticky=N)
    try:
        struktura_= Label(details, text=text_strukture[struktura], font=font_strukture[struktura], bg=background, fg=barve_strukture[struktura])
    except KeyError:
        pass  

    if struktura !='BCC' and struktura != 'FCC+BCC' and 'BCC' in struktura and 'FCC' not in struktura:
        struktura_= Label(details, bg=background, text=struktura, fg='green', font=font_str)
    elif struktura !='FCC' and struktura != 'FCC+BCC' and 'FCC' in struktura and 'BCC' not in struktura:
        struktura_= Label(details, bg=background, text=struktura, fg='#003399', font=font_str)
    elif 'FCC' in struktura and 'BCC' in struktura and struktura != 'FCC+BCC' :
        struktura_= Label(details, bg=background, text=struktura, fg='#6600cc', font=("Arial", "24"))

    struktura_.grid(row=10, column=3, columnspan=3, sticky=N, rowspan=2)
    o=4
    Label(details, text='     ', bg=background, font=font_2, fg=fg_1).grid(row=9+o, column=4, columnspan=3, sticky=NE) 
    Label(details, text='Stability parameters:     ', bg=background, font=font_2, fg=fg_1).grid(row=9+o, column=3, columnspan=3, sticky=NE)

    f={}
    enote=['  %', '  %', '  -', '  kJ/mol', '  J/mol K', '  '+u"\u00b0"+'C', '  -', '  -']
    data=list(data)
    while len(data)> 13:
        del data[0]
    for h,k, l in zip(data[:-1],enote,OPTIONS[1:]):
            f[l]=(h,k)
    g=f.items()
    for j,i in enumerate(g):
        Label(details, text=i[0]+' :  ', bg=background, fg='#0099ff', font=font_1).grid(row=10+j+o, column=3, sticky=E)
        Label(details, text=str(i[1][0])+i[1][1], bg=background, fg=fg_2, font=font_1).grid(row=10+j+o, column=4, sticky=W)

    if bz==0:
        details.geometry('+900+20')
    else:
        cut=pozicija_detajlov.index('+')
        details.geometry(pozicija_detajlov[cut:])
        details.deiconify()
    
def onclick_BIG(event):
    global good_job_big, pozicija_detajlov
    afna= (event.xdata, event.ydata)
    if choice_mode.get()=='basic':
        good_job_big = nearest(BIG_data, afna[0], afna[1])
    else:
        good_job_big = nearest(advanced, afna[0], afna[1])
    #print(good_job_big)
    try:
        pozicija_detajlov=details.winfo_geometry()
        details.destroy()
        get_details()
    except  NameError:
        pass

def onclick_small(event):
    global good_job_small, pozicija_detajlov
    afna= (event.xdata, event.ydata)
    if choice_mode.get()=='basic':
        good_job_small = nearest(results, afna[0], afna[1])
    else:
        good_job_small = nearest(advanced, afna[0], afna[1])
    #print(good_job_small)
    try:
        pozicija_detajlov=details.winfo_geometry()
        details.destroy()
        get_details()
    except  NameError:
        pass

def onclick_substitution(event):
    global good_job_substitution, pozicija_detajlov
    afna= (event.xdata, event.ydata)
    if choice_mode.get()=='basic':
        good_job_substitution = nearest(rezultat_substitution, afna[0], afna[1])
    else:
        good_job_substitution = nearest(advanced, afna[0], afna[1])
    #print(good_job_small)
    try:
        pozicija_detajlov=details.winfo_geometry()
        details.destroy()
        get_details()
    except  NameError:
        pass

import CORRELATIONS
from importlib import reload

Subscripts = str.maketrans("0123456789x-m", '₀₁₂₃₄₅₆₇₈₉ₓ₋ₘ')

def Plot_2D(e):
    switch.set('small_iterator_data')
    #e = [x for x in e if x != []]
    d={}
    for j,i in enumerate(OPTIONS):
        d[i]=j+1
    """
    try:
        plt.close()
    except NameError:
        pass
    """
    barve=['blue','red', 'magenta', '#00ff00', '#006600', '#00e6e6', '#ff9900', '#666699', '#003366', 'yellow', '#cc0000', '#006666', '#9900cc', '#ff9999']
    markerji=['s', 's', 's', '^', 'o', 'o', 'o', '.', 'v', 'v', '>', '>', '<', 's','^','.']
    lojzek={}
    sub=get_subs(); add=pick_addition()
    sf=implicitna_formula(sub,add)
    w_string=sf[1]+sf[2]+'Mx'
    w_string=w_string.translate(Subscripts)
    if choice_mode.get()=='basic':
        plt.ion()
        fig,ax = plt.subplots()
        for i in range(len(e)):
            lojzek["x{}".format(i)]=[n[d[x_varka.get()]] for n in e[i]]
            lojzek["y{}".format(i)]=[m[d[y_varka.get()]] for m in e[i]]
            lojzek["colors"]=barve[i]
            lojzek["markers"]=markerji[i]
            amg=Namespace(**lojzek)
            exec("ax.plot(amg.x{0}, amg.y{0}, label=e[{0}][0][0], c=amg.colors, marker=amg.markers)".format(i))
        enote=['  %', '  %', '  -', '  kJ/mol', '  J/mol K', '  '+u"\u00b0"+'C', '  -', '  -']
        plt.xlabel(x_varka.get()); plt.ylabel(y_varka.get())
    else:
        reload(CORRELATIONS)
        from CORRELATIONS import enka as fun_x
        from CORRELATIONS import dvojka as fun_y
        plt.ion()
        fig,ax = plt.subplots()
        for i in range(len(e)):
            zaloga_x=[fun_x(e[i],n)[0] for n in range(len(e[i])) ]
            zaloga_y=[fun_y(e[i],m)[0] for m in range(len(e[i])) ]
            lojzek["x{}".format(i)]=[fun_x(e[i],n)[0] for n in range(len(e[i])) ]
            lojzek["y{}".format(i)]=[fun_y(e[i],m)[0] for m in range(len(e[i])) ]
            lojzek["colors"]=barve[i]
            lojzek["markers"]=markerji[i]
            amg=Namespace(**lojzek)
            exec("ax.plot(amg.x{0}, amg.y{0}, label='M='+list(e[{0}][0][-2].keys())[-1], c=amg.colors, marker=amg.markers)".format(i))
            n=-1
            for u,v in zip(zaloga_x, zaloga_y):
                n+=1
                e[i][n]=(v,)+ e[i][n]
                e[i][n]=(u,)+ e[i][n]
        #plotSmall_btn.config(state="disabled", fg='light grey')
        plt.xlabel(fun_x(e[0],0)[1], fontsize=15); plt.ylabel(fun_y(e[0],0)[1], fontsize=15)
                
    plt.title(w_string, fontsize=20); plt.legend(loc='best')
    fig.canvas.mpl_connect('button_press_event', onclick_small); plt.show()


def Substitution_plot(e):
    switch.set('substitution_iterator_data')
    e = [x for x in e if x != []]
    d={}
    for j,i in enumerate(OPTIONS):   # Check the options ===>>> Replace Atomic% with X
        d[i]=j+1
    """
    try:
        plt.close()
    except NameError:
        pass
    """
    plt.ion()
    fig,ax = plt.subplots()
    barve=['blue','red', 'magenta', '#00ff00', '#006600', '#00e6e6', '#ff9900', '#666699', '#003366', 'yellow', '#cc0000', '#006666', '#9900cc', '#ff9999']
    markerji=['s', 's', 's', '^', 'o', 'o', 'o', '.', 'v', 'v', '>', '>', '<', 's']
    lojzek={}
    add=pick_addition()
    if bool(add[3]):
        m=[i*add[3][2] for i in range(int(add[3][0]/add[3][2]), int(add[3][1]/add[3][2])+1)]
    if choice_mode.get()=='basic':
        for i in range(len(e)):
            lojzek["x{}".format(i)]=[n[d[x_varka.get()]] for n in e[i]]
            lojzek["y{}".format(i)]=[m[d[y_varka.get()]] for m in e[i]]
            lojzek["colors"]=barve[i]
            lojzek["markers"]=markerji[i]
            amg=Namespace(**lojzek)
            exec("ax.plot(amg.x{0}, amg.y{0}, label=e[{0}][0][0], c=amg.colors, marker=amg.markers)".format(i))
        enote=['  %', '  %', '  -', '  kJ/mol', '  J/mol K', '  '+u"\u00b0"+'C', '  -', '  -']
        plt.xlabel(x_varka.get()); plt.ylabel(y_varka.get())
    else:
        reload(CORRELATIONS)
        from CORRELATIONS import enka as fun_x
        from CORRELATIONS import dvojka as fun_y
        for i in range(len(e)):
            zaloga_x=[fun_x(e[i],n)[0] for n in range(len(e[i])) ]
            zaloga_y=[fun_y(e[i],mm)[0] for mm in range(len(e[i])) ]
            lojzek["x{}".format(i)]=[fun_x(e[i],n)[0] for n in range(len(e[i])) ]
            lojzek["y{}".format(i)]=[fun_y(e[i],m)[0] for m in range(len(e[i])) ]
            lojzek["colors"]=barve[i]
            lojzek["markers"]=markerji[i]
            lojzek["m"]=m[i]
            amg=Namespace(**lojzek)
            if bool(add[3]):
                exec("ax.plot(amg.x{0}, amg.y{0}, label='m='+str(amg.m), c=amg.colors, marker=amg.markers)".format(i))
            else:
                exec("ax.plot(amg.x{0}, amg.y{0}, label='m=0', c=amg.colors, marker=amg.markers)".format(i))
            n=-1
            for u,v in zip(zaloga_x, zaloga_y):
                n+=1
                e[i][n]=(v,)+ e[i][n]
                e[i][n]=(u,)+ e[i][n]
        plt.xlabel(fun_x(e[0],0)[1], fontsize=15); plt.ylabel(fun_y(e[0],0)[1], fontsize=15)
       
    plt.legend(loc='best'); plt.title(Master_formula.get(), fontsize=20)
    fig.canvas.mpl_connect('button_press_event', onclick_substitution)
    plt.show()

def BIG_plot(e):
    switch.set('big_iterator_data')
    e = [x for x in e if x != []]
    d={}
    for j,i in enumerate(OPTIONS):
        if i !='Atomic %':
            d[i]=j-1
    """
    try:
        plt.close()
    except NameError:
        pass
    """
    plt.ion()
    fig,ax = plt.subplots()
    barve=['blue','red', 'magenta', '#00ff00', '#006600', '#00e6e6', '#ff9900', '#666699', '#003366', 'yellow', '#cc0000', '#006666', '#9900cc', '#ff9999']
    markerji=['s', 's', 's', '^', 'o', 'o', 'o', '.', 'v', 'v', '>', '>', '<', 's']
    if choice_mode.get()=='basic':
        iksi=[]; ipsiloni=[]
        for i in range(len(e)):
            iksi += [n[d[x_varka.get()]] for n in e[i]]
            ipsiloni +=[m[d[y_varka.get()]] for m in e[i]]
        ax.scatter(iksi, ipsiloni, c='b', marker='s')
        enote=['  %', '  %', '  -', '  kJ/mol', '  J/mol K', '  '+u"\u00b0"+'C', '  -', '  -']
        enota_x= enote[d[x_varka.get()]]; enota_y= enote[d[y_varka.get()]]
        plt.xlabel(x_varka.get()+'('+enota_x+')'); plt.ylabel(y_varka.get()+'('+enota_y+')'); plt.title(y_varka.get()+'('+enota_y+')  '+' vs. '+x_varka.get()+'('+enota_x+')', fontsize=18)
    else:
        reload(CORRELATIONS)
        from CORRELATIONS import enka as fun_x
        from CORRELATIONS import dvojka as fun_y
        zaloga_x=[]; zaloga_y=[]
        for i in range(len(e)):
            zaloga_x+=[fun_x(e[i],n)[0] for n in range(len(e[i])) ]
            zaloga_y+=[fun_y(e[i],m)[0] for m in range(len(e[i])) ]
        chx = [zaloga_x[x:x+6] for x in range(0, len(zaloga_x), 6)]
        chy = [zaloga_y[x:x+6] for x in range(0, len(zaloga_y), 6)]
        for i in range(len(e)):
            n=-1
            for u,v in zip(chx[i], chy[i]):
                n+=1
                e[i][n]=(v,)+ e[i][n]
                e[i][n]=(u,)+ e[i][n]
        ax.scatter(zaloga_x, zaloga_y, s=20, c='blue', marker='s')
        plt.xlabel(fun_x(e[0],0)[1], fontsize=15); plt.ylabel(fun_y(e[0],0)[1], fontsize=15)
        plt.title('Advanced Plot of Big Iteration', fontsize=20)
        #plotBIG_btn.config(state="disabled", fg='light grey')     
    fig.canvas.mpl_connect('button_press_event', onclick_BIG)
    plt.show()

def Plot_COMB(e):
    switch.set('combination_iterator_data')
    d={}
    for j,i in enumerate(OPTIONS):
        if i !='Atomic %':
            d[i]=j-1
    """
    try:
        plt.close()
    except NameError:
        pass
    """
    
    if choice_mode.get()=='basic':
        plt.ion()
        fig,ax = plt.subplots()
        iksi=[n[d[x_varka.get()]] for n in e]
        ipsiloni=[m[d[y_varka.get()]] for m in e]
        ax.scatter(iksi, ipsiloni, s=20, c='green', marker='o')
        plt.xlabel(x_varka.get()); plt.ylabel(y_varka.get()); plt.title(y_varka.get()+' vs. '+x_varka.get(), fontsize=20)
    else:
        reload(CORRELATIONS)
        from CORRELATIONS import enka as fun_x
        from CORRELATIONS import dvojka as fun_y
        plt.ion()
        fig,ax = plt.subplots()
        zaloga_x=[] ; zaloga_y=[]
        for n in range(len(e)):
            zaloga_x.append(fun_x(e,n)[0])
            zaloga_y.append(fun_y(e,n)[0])
        ax.scatter(zaloga_x, zaloga_y, s=20, c='blue', marker='o')
        n=-1
        for u,v in zip(zaloga_x, zaloga_y):
            n+=1
            e[n]=(v,)+ e[n]
            e[n]=(u,)+ e[n]
        #plotComb_btn.config(state="disabled", fg='light grey')
        plt.xlabel(fun_x(e,0)[1], fontsize=15); plt.ylabel(fun_y(e,0)[1], fontsize=15)
        plt.title('Advanced Plot of Combinations', fontsize=20)
    fig.canvas.mpl_connect('button_press_event', onclick_comb)
    plt.show()

def Plot_2D_together():
    if choice_mode.get()=='basic':
        Plot_2D(results)
    else:
        advanced_plot()

def Plot_Subs_together():
    if choice_mode.get()=='basic':
        Substitution_plot(rezultat_substitution)
    else:
        advanced_plot()

def Plot_BIG_together():
    if choice_mode.get()=='basic':
        BIG_plot(BIG_data)
    else:
        advanced_plot()

def Plot_Comb_together():
    if choice_mode.get()=='basic':
        Plot_COMB(results_comb)
    else:
        advanced_plot()

plotBIG_btn=Button(frame, text='  Plot BIG Iteration  ', width=20, command= Plot_BIG_together); plotBIG_btn.grid(row=24, column=18)
plotSmall_btn=Button(frame, text=' Plot Small Iteration ', width=20, command= Plot_2D_together); plotSmall_btn.grid(row=25, column=18)
plotComb_btn=Button(frame, text='    Plot Combinations   ', width=20, command= Plot_Comb_together); plotComb_btn.grid(row=26, column=18)
plotSubs_btn=Button(frame, text='    Plot Substitution     ', width=20, command= Plot_Subs_together); plotSubs_btn.grid(row=27, column=18)

def parske_entalpije(inpt):
    global fruc
    x=normalize(inpt)
    fruc=pari.COMBINATIONS(x)
    pari_Hmix = {}
    for i in fruc[0]:
        aa= re.sub( r"([A-Z])", r" \1", i).split()
        h_mix = 4*parske_kombinacije[i]*x[aa[0]]*x[aa[1]]
        pari_Hmix[i] = round(h_mix, 2) 
    return pari_Hmix

Button(frame, text='    Show Details   ', command= get_details).grid(row=31, column=17, sticky=S, rowspan=2)
Button(frame, text='       Set Limits     ', command= settings).grid(row=31, column=18, sticky=S, rowspan=2)
Label(frame, textvariable=counter_display, bg=frame_bg, font=('Calibry Body', '16'), fg='#008060').grid(row=33, column=3, columnspan=14, rowspan=2, sticky=N)
Label(frame, textvariable=limit_mode, bg='#47476b', font=('Calibry Body', '10', 'bold'), fg='#ff9900').grid(row=33, column=18, rowspan=2)

def get_the_literature():
    from hea_data_extraction import data_from_literature as majka
    literatura=majka()
    return literatura

q=get_the_literature()

r=[]
for i in range(len(q)):
    try:
        if 'BMG' in q[i][1]:
            pack=0.64       
        elif 'BCC' in q[i][1] and 'FCC' not in q[i][1]:
            pack=0.68      
        elif 'FCC' in q[i][1] and 'BCC' not in q[i][1]:
            pack=0.74     
        elif 'FCC' in q[i][1] and 'BCC' in q[i][1]:
            pack=0.71
        else:
            pack=0.71
            
        r.append(splosni_izracun(q[i][0], pack)[:-1]+ (q[i][1],))
    except KeyError as error:
        #print(error)
        pass

def nearest_literature(array, e,f):
    d={}
    for j,i in enumerate(OPTIONS[1:]):
            d[i]=j
    dist=[]; gt=(1,2,3)
    x=[]; y=[]
    if choice_mode.get()=='modified':
        x += [n[0] for n in array]   
        y += [m[1] for m in array]  
    elif choice_mode.get()=='basic':
        x+= [n[d[x_varka.get()]] for n in array]
        y+= [m[d[y_varka.get()]] for m in array]
    tx=abs(max(x)- min(x))
    ty=abs(max(y)- min(y)) 
    r=ty/tx
    y=[i/r for i in y]
    z=[(i,j)for i,j in zip(x,y)]
    for i in z:
        d=math.sqrt((i[0]-e)**2 + (i[1]-f/r)**2)
        dist.append(d)
    n=dist.index(min(dist))
    return array[n]

def onclick_literature(event):
    global good_job_lit, pozicija_detajlov
    afna= (event.xdata, event.ydata)
    good_job_lit = nearest_literature(joint, afna[0], afna[1])
    try:
        pozicija_detajlov=details.winfo_geometry()
        if choice_mode.get()=='basic':
            #details.destroy()
            pass
        else:
            pass
        get_details()
    except  NameError:
        pass

#import numpy as np

def Literature_plot(e):
    #global fig, ax, sc, calc_x, calc_y, d
    switch.set('total_plot')
    d={}
    for j,i in enumerate(OPTIONS[1:]):
        d[i]=j 
    try:
        plt.close()
    except NameError:
        pass
    barve=['blue','red', 'magenta', '#00ff00', '#006600', '#00e6e6', '#ff9900', '#666699', '#003366', 'yellow', '#cc0000', '#006666', '#9900cc', '#ff9999']
    markerji=['s', 's', 's', '^', 'o', 'o', 'o', '.', 'v', 'v', '>', '>', '<', 's']
    plt.ion()
    fig, ax = plt.subplots()
    #plt.draw()
    iksi_fcc=[]; ipsiloni_fcc=[]; iksi_bcc=[]; ipsiloni_bcc=[]; iksi_fccbcc=[]; ipsiloni_fccbcc=[]; iksi_bmg=[]; ipsiloni_bmg=[] ; calc_x=[] ; calc_y=[]
    iksi_im=[]; ipsiloni_im=[]; BCC_plus_x=[]; BCC_plus_y=[]; FCC_plus_x=[]; FCC_plus_y=[]; FCCBCC_plus_x=[]; FCCBCC_plus_y=[]
    
    for i in range(len(e)):
        if e[i][-1] =='FCC':
            iksi_fcc.append(e[i][d[x_varka.get()]])
            ipsiloni_fcc.append(e[i][d[y_varka.get()]])
        elif e[i][-1] =='BCC':
            iksi_bcc.append(e[i][d[x_varka.get()]])
            ipsiloni_bcc.append(e[i][d[y_varka.get()]])
        elif e[i][-1] =='FCC+BCC' or e[i][-1] =='BCC+FCC':
            iksi_fccbcc.append(e[i][d[x_varka.get()]])
            ipsiloni_fccbcc.append(e[i][d[y_varka.get()]])
        elif e[i][-1] =='BMG':
            iksi_bmg.append(e[i][d[x_varka.get()]])
            ipsiloni_bmg.append(e[i][d[y_varka.get()]])
        elif e[i][-1] =='IM':
            iksi_im.append(e[i][d[x_varka.get()]])
            ipsiloni_im.append(e[i][d[y_varka.get()]])
        elif e[i][-1] !='BCC' and e[i][-1] != 'FCC+BCC' and 'BCC' in e[i][-1] and 'FCC' not in e[i][-1]:
            BCC_plus_x.append(e[i][d[x_varka.get()]])
            BCC_plus_y.append(e[i][d[y_varka.get()]])
        elif e[i][-1] !='FCC' and e[i][-1] != 'FCC+BCC' and 'FCC' in e[i][-1] and 'BCC' not in e[i][-1]:
            FCC_plus_x.append(e[i][d[x_varka.get()]])
            FCC_plus_y.append(e[i][d[y_varka.get()]])
        elif 'FCC' in e[i][-1] and 'BCC' in e[i][-1] and e[i][-1] != 'FCC+BCC' :
            FCCBCC_plus_x.append(e[i][d[x_varka.get()]])
            FCCBCC_plus_y.append(e[i][d[y_varka.get()]])
        else:
            calc_x.append(e[i][d[x_varka.get()]])
            calc_y.append(e[i][d[y_varka.get()]])
            
    ax.scatter(iksi_fcc, ipsiloni_fcc, c='aqua', marker='s', label='FCC')
    ax.scatter(iksi_bcc, ipsiloni_bcc, c='lime', marker='s', label='BCC')
    ax.scatter(iksi_fccbcc, ipsiloni_fccbcc, c='#b366ff', marker='.', label='FCC+BCC', s=100)
    ax.scatter(iksi_bmg, ipsiloni_bmg, c='#cc0000', marker='o', label='BMG')
    ax.scatter(iksi_im, ipsiloni_im, c='orange', marker='o', label='IM')
    ax.scatter(BCC_plus_x, BCC_plus_y, c='green', marker=',', label='BCC plus')
    ax.scatter(FCC_plus_x, FCC_plus_y, c='#003399', marker=',', label='FCC plus')
    ax.scatter(FCCBCC_plus_x, FCCBCC_plus_y, c='#6600cc', marker='.', label='FCC+BCC plus', s=33)
    ax.scatter(calc_x, calc_y, c='magenta', marker='^', label='calculated', s=50)
    #sc=ax.scatter(calc_x, calc_y, c='magenta', marker='^', label='calculated', s=50)
    #sc.set_offsets(np.c_[calc_x,calc_y])
    
    enote=['  %', '  %', '  -', '  kJ/mol', '  J/mol K', '  '+u"\u00b0"+'C', '  -', '  -']
    enota_x= enote[d[x_varka.get()]]; enota_y= enote[d[y_varka.get()]]
    plt.xlabel(x_varka.get()+'('+enota_x+')'); plt.ylabel(y_varka.get()+'('+enota_y+')'); plt.title(y_varka.get()+'('+enota_y+')  '+' vs. '+x_varka.get()+'('+enota_x+')', fontsize=18)
    plt.legend(loc='best')    
    fig.canvas.mpl_connect('button_press_event', onclick_literature)
    #fig.canvas.draw_idle()
    #plt.pause(5)
    #plt.waitforbuttonpress()
    plt.show()

def advanced_plot():
    global advanced
    advanced=[]
    bit=nosilec_podatkov.copy()[-1]
    
    if bit[1]=='results' or bit[1]=='rezultat_substitution' or bit[1]=='BIG_data':
        for i in range(len(bit[0])):
            for j in range(len(bit[0][i])):
                if len(bit[0][i][j])==17:
                    g=list(bit[0][i][j])
                    del g[0]; del g[0]
                    bit[0][i][j]=tuple(g)
        advanced=bit[0].copy()
        
        if bit[1]=='results':
            Plot_2D(advanced)
        elif bit[1]=='BIG_data':
            BIG_plot(advanced)
        else:
            Substitution_plot(advanced)
   
    elif bit[1]=='results_comb':
        ma=bit[0]
        advanced+=ma
        Plot_COMB(advanced)
 
    elif bit[1]=='single_case':
        ma=list(itertools.chain(*bit[0]))#one_shot
        advanced+=bit[0]

def total_plot():
    global joint
    joint=[]
    joint+=r
    for bit in nosilec_podatkov:
        try:
            if bit[1]=='BIG_data':
                if choice_mode.get()=='basic':
                    ma=list(itertools.chain(*bit[0]))
                    ma=[i for i in ma]
                else:
                    for i in range(len(bit[0])):
                        for j in range(len(bit[0][i])):
                            if len(bit[0][i][j])==17:
                                g=list(bit[0][i][j])
                                del g[0]; del g[0]
                                bit[0][i][j]=tuple(g)
                    ma=list(itertools.chain(*bit[0]))
                joint+=ma
        except NameError:
            pass
        try:
            if bit[1]=='results' or bit[1]=='rezultat_substitution':
                if choice_mode.get()=='basic':
                    ma=list(itertools.chain(*bit[0]))
                    ma=[i[2:] for i in ma]
                else:
                    for i in range(len(bit[0])):
                        for j in range(len(bit[0][i])):
                            if len(bit[0][i][j])==17:
                                g=list(bit[0][i][j])
                                del g[0]; del g[0]
                                bit[0][i][j]=tuple(g)
                    ma=list(itertools.chain(*bit[0]))
                joint+=ma
        except NameError:
            pass
        try:
            if bit[1]=='results_comb':
                ma=bit[0]
                joint+=ma         
        except NameError:
            pass
        try:
            if bit[1]=='single_case':
                ma=list(itertools.chain(*bit[0]))#one_shot
                joint+=ma
        except NameError:
            pass 
    if choice_mode.get()=='basic':
        Literature_plot(joint)
    elif choice_mode.get()=='modified':
        corr(joint)
    return joint

def single_shot_plot(z):
    """
    for single in nosilec_podatkov:
        try:
            if single[1]=='single_case':
                ma=list(itertools.chain(*single[0]))
                joint+=ma
        except NameError:
            pass
    """
    yas=list(itertools.chain(*z[0]))
    for i in yas:
        calc_x.append(i[d[x_varka.get()]])
        calc_y.append(i[d[y_varka.get()]])
        sc.set_offsets(np.c_[calc_x,calc_y])
        fig.canvas.draw_idle()
        plt.pause(5)
     
def brisi_joint():
    global joint
    del nosilec_podatkov[:]
    del dolzine[:]
    del one_shot[:]
    nps.config(fg='red')
    nosilec_podatkov_status.set('NO DATA'); print(len(nosilec_podatkov))
    counter_display.set('')
    joint=[]
    joint+=r
    return nosilec_podatkov, joint
    

Button(frame, text='Show All', width=12, command=total_plot).grid(row=24, column=13, columnspan=2, sticky=W) # prikaze vse na nosilcu podatkov
Button(frame, text='Clear Data', width=12, command=brisi_joint).grid(row=24, column=14, columnspan=2, sticky=E) # izbrise nosilec podatkov


def corr(e):
    switch.set('corr_plot')
    reload(CORRELATIONS)
    from CORRELATIONS import enka as fun_x
    from CORRELATIONS import dvojka as fun_y

    plt.ion()
    fig, ax = plt.subplots()
    iksi_fcc=[]; ipsiloni_fcc=[]; iksi_bcc=[]; ipsiloni_bcc=[]; iksi_fccbcc=[]; ipsiloni_fccbcc=[]; iksi_bmg=[]; ipsiloni_bmg=[] ; calc_x=[] ; calc_y=[]
    iksi_im=[]; ipsiloni_im=[]; BCC_plus_x=[]; BCC_plus_y=[]; FCC_plus_x=[]; FCC_plus_y=[];  FCCBCC_plus_x=[] ;  FCCBCC_plus_y=[]
    
    zaloga_x=[] ; zaloga_y=[]
    for i in range(len(e)):
        zaloga_x.append(fun_x(e,i)[0])
        zaloga_y.append(fun_y(e,i)[0])
        if e[i][-1] =='FCC':
            iksi_fcc.append(fun_x(e,i)[0])
            ipsiloni_fcc.append(fun_y(e,i)[0])
        elif e[i][-1] =='BCC':
            iksi_bcc.append(fun_x(e,i)[0])
            ipsiloni_bcc.append(fun_y(e,i)[0])
        elif e[i][-1] =='FCC+BCC' or e[i][-1] =='BCC+FCC':
            iksi_fccbcc.append(fun_x(e,i)[0])
            ipsiloni_fccbcc.append(fun_y(e,i)[0])
        elif e[i][-1] =='BMG':
            iksi_bmg.append(fun_x(e,i)[0])
            ipsiloni_bmg.append(fun_y(e,i)[0])
        elif e[i][-1] =='IM':
            iksi_im.append(fun_x(e,i)[0])
            ipsiloni_im.append(fun_y(e,i)[0])
        elif e[i][-1] !='BCC' and e[i][-1] != 'FCC+BCC' and 'BCC' in e[i][-1] and 'FCC' not in e[i][-1]:
            BCC_plus_x.append(fun_x(e,i)[0])
            BCC_plus_y.append(fun_y(e,i)[0])
        elif e[i][-1] !='FCC' and e[i][-1] != 'FCC+BCC' and 'FCC' in e[i][-1] and 'BCC' not in e[i][-1]:
            FCC_plus_x.append(fun_x(e,i)[0])
            FCC_plus_y.append(fun_y(e,i)[0])
        elif 'FCC' in e[i][-1] and 'BCC' in e[i][-1] and e[i][-1] != 'FCC+BCC' :
            FCCBCC_plus_x.append(fun_x(e,i)[0])
            FCCBCC_plus_y.append(fun_y(e,i)[0])
        else:
            calc_x.append(fun_x(e,i)[0])
            calc_y.append(fun_y(e,i)[0])
            
    ax.scatter(iksi_fcc, ipsiloni_fcc, c='aqua', marker='s', label='FCC')
    ax.scatter(iksi_bcc, ipsiloni_bcc, c='lime', marker='s', label='BCC')
    ax.scatter(iksi_fccbcc, ipsiloni_fccbcc, c='#b366ff', marker='.', label='FCC+BCC', s=100)
    ax.scatter(iksi_bmg, ipsiloni_bmg, c='#cc0000', marker='o', label='BMG')
    ax.scatter(iksi_im, ipsiloni_im, c='orange', marker='o', label='IM')
    ax.scatter(BCC_plus_x, BCC_plus_y, c='green', marker=',', label='BCC plus')
    ax.scatter(FCC_plus_x, FCC_plus_y, c='#003399', marker=',', label='FCC plus')
    ax.scatter(FCCBCC_plus_x, FCCBCC_plus_y, c='#6600cc', marker='.', label='FCC+BCC plus', s=33)
    ax.scatter(calc_x, calc_y, c='magenta', marker='^', label='calculated', s=50)

    nn=-1
    for i,j in zip(zaloga_x, zaloga_y):
        nn+=1
        e[nn] = (j,)+ e[nn]
        e[nn] = (i,)+ e[nn]
        
    plt.xlabel(fun_x(e,0)[1], fontsize=15); plt.ylabel(fun_y(e,0)[1], fontsize=15)  
    plt.legend(loc='best'); plt.title('Advanced Plot of All', fontsize=20)
    fig.canvas.mpl_connect('button_press_event', onclick_literature)
    plt.show()


choice_mode=StringVar(); choice_mode.set('modified')
basic=Radiobutton(frame, text='BASIC', command=lambda: choice_mode.set('basic'), variable=choice_mode, value='basic', bg=frame_bg)
modified=Radiobutton(frame, text='ADVANCED', command=lambda: choice_mode.set('modified'), variable=choice_mode, value='modified', bg=frame_bg)
basic.grid(row=33, column=1)
modified.grid(row=33, column=1, sticky=E)
    
#root.mainloop()

# final label at the bottom of Iterator Machine frame
Label(frame, text='', bg=frame_bg).grid(row=34, column=1)
# ****** THIS IS THE END ***************

