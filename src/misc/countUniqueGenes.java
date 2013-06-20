/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.IOException;
import java.util.HashSet;
import java.util.Set;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author dashazhernakova
 */
public class countUniqueGenes {
    public void count(String probes_fname) throws IOException{
        TextFile probes = new TextFile(probes_fname, false);
        probes.open();
        String[] fields = probes.readLineElems(TextFile.tab);
        int uniqueProbes = 0;
        Set<String> genes = new HashSet<String>();
        while ( (fields = probes.readLineElems(TextFile.tab)) != null ){
            uniqueProbes ++;
            genes.add(fields[16]);
                
        }
        System.out.println("unique probes: " + uniqueProbes);
        System.out.println("unique genes: " + genes.size());
        probes.close();
    }
    
    public static void main(String[] args) throws IOException{
        countUniqueGenes c = new countUniqueGenes();
        c.count("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromTophat_coverageBed_exonwise_Montgomery/Cis-0/eQTLProbesFDR0.05.txt");
        
        
    }
}
