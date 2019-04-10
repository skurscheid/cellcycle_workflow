export baseDir="/home/sebastian/Data/Tremethick/CellCycle"
export bwG1Dir="${baseDir}/ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/deepTools/bamCompare/log2/GRCh38_ensembl84/G1"
export bwMDir="${baseDir}/ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/deepTools/bamCompare/log2/GRCh38_ensembl84/M"
export idrG1Dir="${baseDir}/ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/idr/callpeak_input_only/BEDs/GRCh38_ensembl84/G1"
export idrMDir="${baseDir}/ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/idr/callpeak_input_only/BEDs/GRCh38_ensembl84/M"
export matrixOut="${baseDir}/ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/deepTools/computeMatrix/GRCh38_ensembl84"

export matrixOutM="${baseDir}/ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/deepTools/computeMatrix/GRCh38_ensembl84/M/combined"
export matrixOutG1="${baseDir}/ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/deepTools/computeMatrix/GRCh38_ensembl84/G1/combined"

export pdfOut="${baseDir}/ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/deepTools/plotHeatmap/"

export beforeRegionStartLength="1000"
export afterRegionStartLength="1000"
export threads="32"

# plotting coverage of different libraries across H2A.Z peaks
computeMatrix reference-point --scoreFileName $(ls $bwG1Dir/*.bw|tr "\n" " "; ls $bwMDir/*.bw|tr "\n" " ") \
                              --regionsFileName $idrG1Dir/H2AZG1/H2AZG1_idr.sorted.bed $idrMDir/H2AZM/H2AZM_idr.sorted.bed \
                              --outFileName $matrixOut/matrix_H2AZ_peaks_log2.gz\
                              --outFileSortedRegions $matrixOut/matrix_H2AZ_peaks_log2.bed\
                              --numberOfProcessors 32\
                              --referencePoint center\
                              --beforeRegionStartLength $beforeRegionStartLength\
                              --afterRegionStartLength $afterRegionStartLength\
                              

plotHeatmap --matrixFile $matrixOut/matrix_H2AZ_peaks_log2.gz\
            --outFileName $pdfOut/heatmap_H2AZ_peaks_log2.pdf\
            --averageTypeSummaryPlot mean\
            --refPointLabel "Peak center"\

# H2A.Z has extremely high number of peaks (> 100k)
# removing it for visualisation of other libraries
computeMatrix reference-point --scoreFileName $(ls $bwMDir/*.bw|tr "\n" " ") \
                              --regionsFileName $(ls $idrMDir/*/*idr.bed | grep -v "H2AZ" | tr "\n" " ")\
                              --outFileName $matrixOutM/matrix_wo_H2AZ.gz\
                              --outFileSortedRegions $matrixOutM/matrix_wo_H2AZ.bed\
                              --referencePoint center\
                              --beforeRegionStartLength $beforeRegionStartLength\
                              --afterRegionStartLength $afterRegionStartLength\
                              --numberOfProcessors $threads

plotHeatmap --matrixFile $matrixOutM/matrix_wo_H2AZ.gz\
            --outFileName heatmap_M_wo_H2AZ.pdf\
            --averageTypeSummaryPlot mean\
            --refPointLabel "Peak center"\


computeMatrix reference-point --scoreFileName $(ls $bwG1Dir/*.bw|tr "\n" " ") \
                              --regionsFileName $(ls $idrG1Dir/*/*idr.bed | grep -v "H2AZ" | tr "\n" " ")\
                              --outFileName $matrixOutG1/matrix_wo_H2AZ.gz\
                              --outFileSortedRegions $matrixOutG1/matrix_wo_H2AZ.bed\
                              --referencePoint center\
                              --beforeRegionStartLength $beforeRegionStartLength\
                              --afterRegionStartLength $afterRegionStartLength\
                              --numberOfProcessors $threads

plotHeatmap --matrixFile $matrixOutG1/matrix_wo_H2AZ.gz\
            --outFileName heatmap_G1_wo_H2AZ.pdf\
            --averageTypeSummaryPlot mean\
            --refPointLabel "Peak center"\

computeMatrix reference-point --scoreFileName $(ls $bwMDir/*.bw|tr "\n" " ") \
                              --regionsFileName $(ls $idrG1Dir/*/*idr.bed | grep -v "H2AZ" | tr "\n" " ")\
                              --outFileName $matrixOutG1/matrix_wo_H2AZ_M_peaks.gz\
                              --outFileSortedRegions $matrixOutG1/matrix_wo_H2AZ_M_peaks.bed\
                              --referencePoint center\
                              --beforeRegionStartLength $beforeRegionStartLength\
                              --afterRegionStartLength $afterRegionStartLength\
                              --numberOfProcessors $threads

plotHeatmap --matrixFile $matrixOutG1/matrix_wo_H2AZ_M_peaks.gz\
            --outFileName heatmap_G1_wo_H2AZ_M_peaks.pdf\
            --averageTypeSummaryPlot mean\
            --refPointLabel "Peak center"


computeMatrix reference-point --scoreFileName $(ls $bwG1Dir/*.bw|tr "\n" " ") \
                              --regionsFileName $(ls $idrMDir/*/*idr.bed | grep -v "H2AZ" | tr "\n" " ")\
                              --outFileName $matrixOutM/matrix_wo_H2AZ_G1_peaks.gz\
                              --outFileSortedRegions $matrixOutM/matrix_wo_H2AZ_G1_peaks.bed\
                              --referencePoint center\
                              --beforeRegionStartLength $beforeRegionStartLength\
                              --afterRegionStartLength $afterRegionStartLength\
                              --numberOfProcessors $threads

plotHeatmap --matrixFile $matrixOutM/matrix_wo_H2AZ_G1_peaks.gz\
            --outFileName heatmap_M_wo_H2AZ_G1_peaks.pdf\
            --averageTypeSummaryPlot mean\
            --refPointLabel "Peak center"
