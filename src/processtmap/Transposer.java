/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package processtmap;

import java.io.File;
import java.io.IOException;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author dashazhernakova
 */
public class Transposer {
    String[][] matrix;
    int numL, len;
    
    public Transposer(){
        
    }
    public void processDir(String dir_name, String ext) throws IOException{
        File dir = new File(dir_name);
        for (File ch : dir.listFiles()) {
            if (ch.getName().endsWith(ext)){
                System.out.println("processing " + ch.getPath());
                readFromFile(ch.getPath());
                writeTransposed(ch.getPath().replace(ext, "_transp" + ext), " ");
            }
        }

    }
    public void readFromFile(String fname) throws IOException{
        TextFile t = new TextFile(fname, false);
        t.open();
        
        numL = t.countLines();
        String[] els = t.readLineElems(TextFile.space);
        len = els.length;
        int lineNum = 0;
        matrix = new String[numL][len];
        while (lineNum < numL){
            for (int i = 0; i < len; i++){
                if (els.length < len)
                    System.out.println("panika");
                matrix[lineNum][i] = els[i];
            }
            els = t.readLineElems(TextFile.space);
            lineNum++;
        }
        t.close();
    }
    public void writeTransposed(String fname, String sep)throws IOException{
        TextFile t = new TextFile(fname, true);
        t.open();
        for (int col = 0; col < len; col++){
            t.write(matrix[0][col]);
            for (int row = 1; row < numL; row++){
                t.write(sep+ matrix[row][col]);
            }
            t.writeln();
        }
        t.close();
        
    }
    public static void main(String[] args) throws IOException{
        Transposer t = new Transposer();
        //t.readFromFile("/Users/dashazhernakova/Documents/UMCG/tmp_snps.txt");
        //t.writeTransposed("/Users/dashazhernakova/Documents/UMCG/tmp_snps_new.txt", "\t");
        t.processDir("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/HapMap3YRI", ".ped");
    }
}
