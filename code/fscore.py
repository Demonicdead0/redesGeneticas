import os
import numpy as np
from sklearn.metrics import f1_score
from matplotlib import pyplot as plt

def leer_tabla(ruta_archivo: str):
    with open(ruta_archivo, 'r', encoding='utf8') as f:
        lineas = f.readlines()[1:]
        tabla = [list(map(int, linea.split())) for linea in lineas]
    return np.array(tabla)

def fscore(tabla1, tabla2):
    etiquetas1 = tabla1.flatten()
    etiquetas2 = tabla2.flatten()

    return f1_score(etiquetas1, etiquetas2, average='weighted')

carpeta: str = "../data/"

datos_verdaderos: list = []

for archivo in os.listdir(carpeta):
    if archivo.startswith('deg2'):
        ruta_archivo = os.path.join(carpeta, archivo);
        if os.path.isfile(ruta_archivo):
            datos_verdaderos.append(ruta_archivo)

carpeta: str = "../res/bf/"

lista_bruteforce: list = []
for archivo in os.listdir('../res/bf/'):
    ruta_archivo = os.path.join(carpeta, archivo)
    lista_bruteforce.append(ruta_archivo)

carpeta: str = "../res/sfs/"
lista_sfs: list = []
for archivo in os.listdir('../res/sfs/'):
    ruta_archivo = os.path.join(carpeta, archivo)
    lista_sfs.append(ruta_archivo)

carpeta: str = "../res/bfb/"
lista_bfb: list = []
for archivo in os.listdir('../res/bfb/'):
    ruta_archivo = os.path.join(carpeta, archivo)
    lista_bfb.append(ruta_archivo)

carpeta: str = "../res/dp/"
lista_dp: list = []
for archivo in os.listdir('../res/dp/'):
    ruta_archivo = os.path.join(carpeta, archivo)
    lista_dp.append(ruta_archivo)


carpeta: str = "../res/bfb2/"
lista_bfb2: list = []
for archivo in os.listdir('../res/bfb2/'):
    ruta_archivo = os.path.join(carpeta, archivo)
    lista_bfb2.append(ruta_archivo)

print(len(datos_verdaderos))
print(len(lista_bruteforce))
print(len(lista_sfs))
print(len(lista_bfb))
print(len(lista_bfb2))
print(len(lista_dp))


x: list = []
y1: list = []
y2: list = []
y3: list = []
y4: list = []
y5: list = []

for i in range(len(datos_verdaderos)):
    x.append(i)
    tabla1 = leer_tabla(datos_verdaderos[i])
    #tabla2 = leer_tabla(lista_bruteforce[i])
    #tabla3 = leer_tabla(lista_sfs[i])
    tabla4 = leer_tabla(lista_bfb[i])
    tabla5 = leer_tabla(lista_dp[i])
    tabla6 = leer_tabla(lista_bfb2[i])
    #y1.append(fscore(tabla1, tabla2))
    #y2.append(fscore(tabla1, tabla3))
    y3.append(fscore(tabla1, tabla4))
    y4.append(fscore(tabla1, tabla5))
    y5.append(fscore(tabla1, tabla6))


print("Datos procesados")
x = np.array(x)

# Ancho de las barras
width = 0.35  

fig, ax = plt.subplots()

n_bins = 200
# Desplazar las barras y graficarlas una al lado de otra
plt.hist(y4, bins=n_bins, alpha=0.7, label='DP MI')
#plt.hist(y1, bins=n_bins, alpha=0.7, label='fuerza Bruta reducidad MI')
#plt.hist(y2, bins=n_bins, alpha=0.7, label='SFS MI')
plt.hist(y3, bins=n_bins, alpha=0.7, label='SBS MI')
plt.hist(y5, bins=n_bins, alpha=0.7, label='SBS MIP')


plt.title('Histograma para visualizar distribución')
plt.xlabel('Valores')
plt.ylabel('Frecuencia')
plt.legend()

plt.show()