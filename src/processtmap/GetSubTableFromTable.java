/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package processtmap;

import java.io.IOException;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author dashazhernakova
 */
public class GetSubTableFromTable {
    public void get(int[] posToLeave, String in_f, String out_f) throws IOException{
        TextFile in = new TextFile(in_f, false);
        TextFile out = new TextFile(out_f, true);
        in.open();
        out.open();
        String[] els;
        
        while( (els = in.readLineElems(TextFile.tab)) != null ){
            String new_str = "";
            for (int i = 0; i < posToLeave.length; i++)
                new_str += "\t" + els[posToLeave[i]];
            out.writeln(new_str.replaceFirst("\t", ""));
        }
        in.close();
        out.close();
    }
    
    
    public static void main(String[] args) throws IOException {
        GetSubTableFromTable g = new GetSubTableFromTable();
        int[] x = {2,3};
        g.get(x ,"/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/MontgomeryGenotypes/RNASEQ60_snps.full.txt" , "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/MontgomeryGenotypes/RNASEQ60_snps.full.map");
    }
}
