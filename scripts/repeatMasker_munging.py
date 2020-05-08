import pandas as pd
import os, uuid, sys

def writeToBed(df, selected_columns, output):
    df[selected_columns].to_csv(output, sep =  "\t", index = False, header = False)

selected_columns = ['genoName', 'genoStart', 'chromEnd', 'name', 'score', 'strand']

col_names = ['#bin', 'swScore', 'milliDiv', 'milliDel', 'milliIns', 'genoName', 
             'genoStart', 'genoLeft', 'strand', 'repName', 'repClass','repFamily', 
             'repStart', 'repEnd', 'repLeft', 'id']

repeat_file = 'repeatMaskerhg38.gz'

df = pd.read_csv(repeat_file, sep="\t", header = 0, names = col_names)
