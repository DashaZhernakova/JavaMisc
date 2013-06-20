/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package processtmap;

import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeMap;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author dashazhernakova
 */
public class GetHalfExpressionTable {
    float threshold;
    String expr_file;
    public GetHalfExpressionTable(String f_name){
        expr_file = f_name;
    }

    public void calculateAverages() throws IOException{
        TextFile expr = new TextFile(expr_file, false);
        expr.open();
        String[] fields= expr.readLineElems(TextFile.tab);
        int len = fields.length;
        TreeMap<Float, ArrayList<String>> avgs = new TreeMap<Float,ArrayList<String>>();
        while ( (fields = expr.readLineElems(TextFile.tab)) != null ){
            float sum = 0;
            for (int i = 1; i < len; i++)
                sum += Float.valueOf(fields[i]);
            float avg = sum / (len - 1);
            
        }
        expr.close();
    }
    public static void main(String[] args) throws IOException {
        GetHalfExpressionTable ht = new GetHalfExpressionTable("");
    }
}
