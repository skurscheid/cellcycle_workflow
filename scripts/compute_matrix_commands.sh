export bwG1Dir="/data/experiment/ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/deepTools/bamCompare/GRCh38_ensembl84/G1"
export bwMDir="/data/experiment/ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/deepTools/bamCompare/GRCh38_ensembl84/M"
export idrG1Dir="/data/experiment/ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/idr/BEDs/GRCh38_ensembl84/G1"
export idrMDir="/data/experiment/ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/idr/BEDs/GRCh38_ensembl84/M"
export matrixOutM="/data/experiment/ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/deepTools/computeMatrix/GRCh38_ensembl84/M/combined"
export matrixOutG1="/data/experiment/ChIP-Seq/LR1807201/N08851_SK_LR1807201_SEQ/deepTools/computeMatrix/GRCh38_ensembl84/G1/combined"

export beforeRegionStartLength="1000"
export afterRegionStartLength="1000"
export threads="8"

computeMatrix reference-point --scoreFileName $(ls $bwMDir/*.bw|tr "\n" " ") \
                              --regionsFileName $(ls $idrMDir/*/*idr.bed|tr "\n" " ")\
                              --outFileName $matrixOutM/matrix.gz\
                              --outFileSortedRegions $matrixOutM/matrix.bed\
                              --referencePoint center\
                              --beforeRegionStartLength $beforeRegionStartLength\
                              --afterRegionStartLength $afterRegionStartLength\
                              --numberOfProcessors $threads

plotHeatmap --matrixFile $matrixOutM/matrix.gz\
            --outFileName heatmap.pdf\
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
