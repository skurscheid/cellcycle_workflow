import pandas as pd
import os, uuid, sys

idr_cutoff = snakemake.params.globalIDRCutoff
signal_value = snakemake.params.signalValue
idr_file = snakemake.input.idr_file
selected_columns = ['chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand']
idr_peaks = snakemake.output.idr_peaks
other_peaks = snakemake.output.other_peaks

def writeToBed(df, selected_columns, output):
    df[selected_columns].to_csv(output, sep =  "\t", index = False, header = False)

def run_script():
    try:
        df = pd.read_csv(idr_file, sep="\t", header = None, names = ["chrom", "chromStart", "chromEnd", "name", "score", "strand",
                                                             "signalValue", "p-value", "q-value", "summit", "localIDR", "globalIDR",
                                                             "rep1_chromStart", "rep1_chromEnd", "rep1_signalValue", "rep1_summit",
                                                             "rep2_chromStart", "rep2_chromEnd", "rep2_signalValue", "rep2_summit"])
        
        other = df[(df['signalValue'] >= signal_value) & (df['globalIDR'] <= idr_cutoff)]
        idr = df[df["globalIDR"] > idr_cutoff]
        writeToBed(other, selected_columns, other_peaks)
        writeToBed(idr, selected_columns, idr_peaks)
    except Exception as e:
        print(e)

# Main method.
if __name__ == '__main__':
    run_script()