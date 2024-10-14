import numpy 
import pandas as pd
import os

pathData = "/home/usuario/Documentos/GitHub/CaNS/run/data/"
df = pd.read_csv(pathData + 'prueba.txt')

column_counts = df.count(axis=1)
print(column_counts)