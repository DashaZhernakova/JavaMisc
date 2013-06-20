package misc;

import com.lowagie.text.DocumentException;
import java.io.IOException;
import umcg.genetica.graphics.Heatmap;
import umcg.genetica.math.matrix.DoubleMatrixDataset;

/**
 *
 * @author dashazhernakova
 */
public class DrawHeatMap {
    
    public static void main(String[] args) throws IOException, DocumentException{
        
        //String fName = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/deepSAGE_transcr/NORMALIZED/lincRNA/expresion_table_tmp.txt", 
        //        out = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/deepSAGE_transcr/NORMALIZED/lincRNA/heatmap.pdf";
        String fName = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/lincRNA_Yale+Argonne/expression_table_all_yale+argonne.txt.expressedInAllSamples.txt.200genes.sortedByName.txt.QuantileNormalized.Log2Transformed.txt.gz",
                out = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/lincRNA_Yale+Argonne/heatmap.QuantileNormalized.Log2Transformed_200_sortedByName.pdf";
        //String fName = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/lincRNA_Montgomery/expression_table_all.txt.expressedInAllSamples.txt.200genes.sortedByName.txt.QuantileNormalized.Log2Transformed.txt.gz",
        //        out = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/lincRNA_Montgomery/heatmap.QuantileNormalized.Log2Transformed_200_sortedByName.pdf";
        //String fName = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/deepSAGE_transcr/NORMALIZED/lincRNA/expresion_table.txt.expressedInAllSamples.txt.200genes.sortedByName.txt.QuantileNormalized.Log2Transformed.txt.gz",
          //      out = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/deepSAGE_transcr/NORMALIZED/lincRNA/heatmap.QuantileNormalized.Log2Transformed_200_sortedByName.pdf";
        DoubleMatrixDataset matrix = new DoubleMatrixDataset(fName);
        String[] rowNames = new String[matrix.nrRows], colNames = new String[matrix.nrCols];
        
        for (int i = 0; i < matrix.nrRows; i++){
            rowNames[i] = matrix.rowObjects.get(i).toString();
        }
        for (int i = 0; i < matrix.nrCols; i++){
            colNames[i] = matrix.colObjects.get(i).toString();
        }
        Heatmap.drawHeatmap(matrix.rawData, rowNames, colNames, 1000, 2700, out, Heatmap.Output.PDF);
        //Heatmap.drawHeatmap(matrix.rawData, rowNames, colNames, 800, 2200, out, Heatmap.Output.PDF);
        //Heatmap.drawHeatmap(matrix.rawData, rowNames, colNames, 1000, 1000, out, Heatmap.Output.PDF);
    }
}
