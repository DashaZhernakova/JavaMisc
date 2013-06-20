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
public class PutExonInAnnotation {
 public void put(String f_n) throws IOException{
     TextFile in = new TextFile(f_n, false);
     TextFile out = new TextFile(f_n+"_new", true);
     in.open();
     out.open();
     String[] els = in.readLineElems(TextFile.tab);
     out.writeln("platform\tprobe\tgene\tchr\tstart\tend");
     while ( (els = in.readLineElems(TextFile.tab)) != null ){
         out.writeln("exonwise\t" + els[0] + ":exon:" + els[1] + "\t" + els[2] + "\t" + els[3]+ "\t" + els[4]+ "\t" + els[5]);
     }
     in.close();
     out.close();
 }   
 public static void main(String[] args) throws IOException {
     PutExonInAnnotation p = new PutExonInAnnotation();
     p.put("/Users/dashazhernakova/Documents/UMCG/hg19/exons_hg19_ranks.txt");        
 }
}
