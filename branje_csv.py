raw_data = []

g=open('podatki_o_atomih.csv', 'r')
for vrstica in g.readlines()[2:]:
     narazen = vrstica.split(',')
     raw_data.append(narazen)

for j,i in enumerate(raw_data):
    raw_data[j] = i[0].split(";")
    raw_data[j][-1]= raw_data[j][-1].replace("\n","")

raw_data=raw_data[:78]

radii={}; kappa_data_Pauling={};  VEC_data={}; Tm_data={};e_a_data={};molske_data={}; density_data={}; kappa_data_Allen={}
for i in raw_data:
     radii[i[1]]=float(i[3])
     kappa_data_Pauling[i[1]]=float(i[4])
     VEC_data[i[1]]=float(i[5])
     Tm_data[i[1]]=int(i[6])
     e_a_data[i[1]]= int(i[7])
     molske_data[i[1]]= float(i[8])
     density_data[i[1]]=float(i[9])
     try:
          kappa_data_Allen[i[1]]=float(i[10])
     except IndexError:
          kappa_data_Allen[i[1]]=float(i[4])

""" PARSKE KOMBINACIJE """

baza_parskih = []

g=open('parske_entalpije.csv', 'r')
for vrstica in g.readlines():
     narazen = vrstica.split(',')
     baza_parskih.append(narazen)

for j,i in enumerate(baza_parskih):
    baza_parskih[j] = i[0].split(";")
    baza_parskih[j][-1]= baza_parskih[j][-1].replace("\n","")

parske_kombinacije = {}; entalpije_trokarevski={}

for i in baza_parskih:
     parske_kombinacije[i[0]+i[1]]=int(i[2])
     parske_kombinacije[i[1]+i[0]]=int(i[2])

for i in baza_parskih:
     try:
          entalpije_trokarevski[i[0]+i[1]]=int(i[4])
          entalpije_trokarevski[i[1]+i[0]]=int(i[4])
     except IndexError:
          pass
