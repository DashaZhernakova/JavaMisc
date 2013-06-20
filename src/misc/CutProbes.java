/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author dashazhernakova
 */
public class CutProbes {
    public void cut(String f_name) throws IOException{
        TextFile in = new TextFile(f_name, false);
        TextFile out = new TextFile(f_name.replace(".txt", "_cutProbes.txt"), true);
        in.open();
        out.open();
        
        boolean notZero;
        String[] els = in.readLineElems(TextFile.tab);
        out.writelnTabDelimited(els);
        int len = els.length;
        while( (els = in.readLineElems(TextFile.tab)) != null ){
            notZero = false;

            for (int i = 1; i < len; i++){
                float cur = Float.valueOf(els[i]);
                if ( cur != 0){
                    notZero = true;
                    break;
                }
            }
            if (notZero)
                out.writelnTabDelimited(els);
        }
        
        in.close();
        out.close();
        
    }
    public static void main(String[] args) throws IOException {
        CutProbes cp = new CutProbes();
        cp.cut("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromTophat_coverageBed_exonwise_Pickrell_yale/expression_table.txt");
        
    }
}
