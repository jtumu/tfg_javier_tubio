#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Javier Tubio Muñoz
"""

# Randomizador
# Para usar el script:
    # 1. cambiar el os.chdir a donde esté tu archivo multifasta
    # 2. cambiar la seed (semilla)
    # 3. Cambiar el nombre a Organismo_random_seed

# El script toma las proteínas del proteoma, las aleatoriza cambiando de orden los 
# aminoácidos, manteniendo la longitud y los aminoácidos originales. Mantiene la
# Metionina en primer lugar siempre.
import os
os.chdir("/home/jtb/Desktop/TFG/")
print("Carpeta actual:", os.getcwd())
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


ARCHIVO_ENTRADA = "s_cerevisiae.fasta"
ARCHIVO_SALIDA = "S_Cerevisiae_random_555"
seed =555

random.seed(seed)

def randomizar_manteniendo_M(secuencia_str):
 
    if secuencia_str.startswith("M"):
        inicio = secuencia_str[0]
        resto = list(secuencia_str[1:])
        random.shuffle(resto)
        return inicio + "".join(resto)
    else:
        lista_completa = list(secuencia_str)
        random.shuffle(lista_completa)
        return "".join(lista_completa)

contador_original = 0
for record in SeqIO.parse(ARCHIVO_ENTRADA, "fasta"):
    contador_original += 1
    
print(f"Hay {contador_original} secuencias en '{ARCHIVO_ENTRADA}'.")

registros_nuevos = []
for record in SeqIO.parse(ARCHIVO_ENTRADA, "fasta"):
    
    seq_original = str(record.seq)
    
    seq_nueva = randomizar_manteniendo_M(seq_original)
    
    nuevo_record = SeqRecord(
        Seq(seq_nueva),
        id=record.id,
        description=record.description,
        name=record.name
    )
    
    registros_nuevos.append(nuevo_record)


count = SeqIO.write(registros_nuevos, ARCHIVO_SALIDA, "fasta")

print(f"Se han randomizado y guardado {count} secuencias en '{ARCHIVO_SALIDA}'.")





