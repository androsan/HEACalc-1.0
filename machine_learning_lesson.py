""" Program za učenje strojenga učenja. Vir:
     https://machinelearningmastery.com/machine-learning-in-python-step-by-step/
"""

# load libraries
import pandas
import seaborn
from pandas.plotting import scatter_matrix
#import matplotlib.pyplot as plt
from sklearn import model_selection
from sklearn.metrics import classification_report
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
import pari_machine_learning
import matplotlib.pyplot as plt

# TEST databases @ UCI Machine Learning Repository:

# https://archive.ics.uci.edu/ml/machine-learning-databases/ionosphere/ionosphere.data     IONOSPHERE
#"https://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data"   IRIS

# 2.2 Load dataset
#url = "https://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data"
#names = ['sepal-length', 'sepal-width', 'petal-length', 'petal-width', 'class']
#dataset = pandas.read_csv(url, names=names)

""" --------------------------------------------------------------------------------------------- """
def read_data(file):
    raw_data = []
    g=open(file, 'r')
    for vrstica in g.readlines()[:5]:
         narazen = vrstica.split(',')
         raw_data.append(narazen)

    for j,i in enumerate(raw_data):
        raw_data[j] = i[0].split(";")
        raw_data[j][-1]= raw_data[j][-1].replace("\n","")

    #raw_data=raw_data[:5]
    return raw_data
#q=read_data('hea.csv')

imena=['Elements/family' , 'Composition', 'Structure', 'Ref.', 'δ%', 'Hmix kJ/mol', 'Δχ%', 'Ω', 'μ', 'VEC', 'e/a', 'sm%', 'Km%']
file=r'hea3.csv'
dataset=pandas.read_csv(file, names=imena)[:322]

del dataset['μ']
del dataset['e/a']
del dataset['sm%']
del dataset['Km%']

imena=['Elements/family' , 'Composition', 'Structure', 'Ref.', 'δ%', 'Hmix kJ/mol', 'Δχ%', 'Ω', 'VEC']
"""
a=[i for i in range(20)]
b=[2*i for i in a]
c=[]
for i in a:
    if i<10:
        c.append('x')
    else:
        c.append('y')

d={'A':a, 'B':b, 'result':c}
df=pandas.DataFrame(data=d)
imena=['prva_kolona','druga_kolona','rezultat']
"""

""" --------------------------------------------------------------------------------------------- """

# 3.1 Dimensions of dataset (shape)
# print(dataset.shape)

# 3.2 Quick look at dataset (head(number_of_rows))
# print(dataset.head(20))# first 20 rows of dataset

# 3.3 STATISTICAL SUMMARY (count, mean, std, min, percentiles, max)
# print(dataset.describe())

# 3.4 Class distribution (število vrstic, ki opisujejo posamezno kategorijo, npr. HEA, BMG, IM)
# print(dataset.groupby('class').size())
"""-----------------------------------"""
# Filter Dataframe
fildf = pandas.DataFrame(columns=imena)
def FD(baza, ime_kolone, imena_atributov):        # "imena_atributov" je list, npr. ['BCC', 'FCC', 'FCC+BCC']
    internal=pandas.DataFrame(columns=imena)
    for j in imena_atributov:
        for i in range(len(baza)):
                if baza.iloc[i][ime_kolone]== j:
                        fildf.loc[i] = baza.iloc[i]
                        internal.loc[i] = baza.iloc[i]
    return fildf, internal

#FD(dataset, 'Structure', ['FCC'])#, 'FCC+BCC', 'BCC'] )

def izberi_kategorije(baza, ime_kolone, izberi_kategorije):
    frames=[]
    for i in izberi_kategorije:
        rezultat= FD(baza, ime_kolone, [i])
        rezultat=rezultat[1][pandas.notnull(rezultat[1][:])]
        frames.append(rezultat)
    kategorije= pandas.concat(frames)
    return kategorije

izbor= ['FCC','BCC','FCC+BCC','BMG'] #['FCC','BCC','FCC+BCC','BMG']
datas= izberi_kategorije(dataset, 'Structure', izbor)
#fildf=fildf[pandas.notnull(fildf[:])]       # Izbriše vse vrstice, ki na kateremkoli mestu vsebujejo NaN (i.e. "not available")
"""-----------------------------------"""
# 4.1 UNIVARIATE PLOTS
def SP(baze, x_ime, y_ime):         # "baze" je list, ki ga dobiš s filtriranjem atributov(npr. 'FCC', 'BCC', itd.)glavne baze 
    # single plot
    barve=['blue','red','green', 'magenta', 'orange', 'black']
    legenda=[]
    for barva, j in zip(barve,baze):
        x=[j.iloc[i][x_ime] for i in range(len(j))]
        y=[j.iloc[i][y_ime] for i in range(len(j))]
        t=plt.scatter(x,y,lw=0, alpha=0.5, facecolor=barva)
        legenda.append(t)
    plt.legend(tuple(legenda),tuple(['FCC','BCC','FCC+BCC','BMG']))
    plt.xlabel(x_ime)
    plt.ylabel(y_ime)
    plt.show()

#for i in pari_machine_learning.novo:
#	SP(datas, i[0], i[1])

def box_and_whisker_plots(podatki):
    # box and whisker plots
    podatki.plot(kind='box', subplots=True, layout=(2,2), sharex=False, sharey=False)
    plt.show()

def histogrami(podatki):
    # histograms
    podatki.hist()
    plt.show()

# 4.2 MULTIVARIATE PLOTS
def sipanje(podatki):
    # scatter plot matrix
    scatter_matrix(podatki)
    plt.show()

def sipanje_multi(podatki, ime_kolone):
    # scatter plot matrix with multiple input data
    seaborn.pairplot(podatki, hue=ime_kolone, size=3, diag_kind="kde")
    plt.show()

sipanje_multi(datas, 'Structure')

# 5.1 CREATE A VALIDATION DATASET

def CVD(podatki, x_min,x_max, y, validation_size, seed):     # CVD = Create Validation Dataset
    # Split-out validation dataset 
    array = podatki.values
    X = array[:,x_min:x_max]     # x_min = 0,  x_max=4   (x pomeni številka kolone v podatkih za vhodni parameter, npr. delta, dHmix,...)
    Y = array[:,y]                        # y = 4                              (y pomeni številka kolone v podatkih za razultat, npr. HEA, BMG, IM)
    #validation_size = 0.20
    #seed = 7
    X_train, X_validation, Y_train, Y_validation = model_selection.train_test_split(X, Y, test_size=validation_size, random_state=seed)
    return X_train, X_validation, Y_train, Y_validation

Validacija = CVD(datas, 4,9,2,0.2,7)
#Validacija = CVD(df, 0,2,2,0.2,7)


# 5.2 TEST HARNESS : Test options and evaluation metric
seed = 7
scoring = 'accuracy'

# 5.3 BUILD MODELS

def build_models(X_train, Y_train, splits_num, scoring):
    # Spot Check Algorithms
    models = []
    models.append(('LR', LogisticRegression()))
    models.append(('LDA', LinearDiscriminantAnalysis()))
    models.append(('KNN', KNeighborsClassifier()))
    models.append(('CART', DecisionTreeClassifier()))
    models.append(('NB', GaussianNB()))
    models.append(('SVM', SVC()))
    # evaluate each model in turn
    results = []
    names = []
    for name, model in models:
            kfold = model_selection.KFold(n_splits=splits_num, random_state=seed)  # splits=10 , seed=7
            cv_results = model_selection.cross_val_score(model, X_train, Y_train, cv=kfold, scoring=scoring)
            results.append(cv_results)
            names.append(name)
            msg = "%s: %f (%f)" % (name, cv_results.mean(), cv_results.std())
            print(msg)
    return results, names, models

Modeli = build_models(Validacija[0], Validacija[2], 15, scoring)

# 5.4 SELECT BEST MODEL
def SBM(results, names):
    # Compare Algorithms
    fig = plt.figure()
    fig.suptitle('Algorithm Comparison')
    ax = fig.add_subplot(111)
    plt.boxplot(results)
    ax.set_xticklabels(names)
    plt.show()

Izberi_najboljsi_model = SBM(Modeli[0], Modeli[1])

# 6. MAKE PREDICTIONS
def make_predictions(X_train, X_validation,Y_train,  Y_validation):
    # Make predictions on validation dataset
    best_model = DecisionTreeClassifier()#KNeighborsClassifier()
    best_model.fit(X_train, Y_train)
    predictions = best_model.predict(X_validation)
    print(accuracy_score(Y_validation, predictions))
    print(confusion_matrix(Y_validation, predictions))
    print(classification_report(Y_validation, predictions))
    return best_model

Predpostavke = make_predictions(Validacija[0],Validacija[1],Validacija[2],Validacija[3])

def napoved(vnesi):
    napoved=Predpostavke.predict(vnesi)
    return napoved


### NEODVISNI PODATKI """ - -  -    -        -
def independence(od, do):
    a=[i for i in range(od,do)]
    b=[2*i for i in a]
    test={'A':a, 'B':b}
    df=pandas.DataFrame(data=test)
    test = df.values
    return test

def insert_HEA_parameters(a1, b1, c1, d1, e1):
    d={'δ%':[a1], 'Hmix kJ/mol':[b1], 'Δχ%':[c1], 'Ω':[d1], 'VEC':[e1]}
    df=pandas.DataFrame(data=d, columns=d.keys())    # ukaz " columns=d.keys()" ohrani vrstni red pri pretvorbi iz slovarja v DataFrame
    d=df.values
    n=napoved(d)
    return n




