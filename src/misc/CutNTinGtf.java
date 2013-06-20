/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

/**
 *
 * @author dashazhernakova
 */
public class CutNTinGtf {
    public void cut() throws FileNotFoundException, IOException{
        String eol = System.getProperty( "line.separator" );
        BufferedReader gtf = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/UMCG/hg19/Homo_sapiens.GRCh37.65_cut.gtf"));
        BufferedWriter new_gtf = new BufferedWriter (new FileWriter("/Users/dashazhernakova/Documents/UMCG/hg19/Homo_sapiens.GRCh37.65_NTcut.gtf"));
        
        String line= gtf.readLine();
        while ((line) != null){
            if ((line.startsWith("X")) || (line.startsWith("Y")) || (line.startsWith("MT")) || (Character.isDigit(line.charAt(0)))){
                new_gtf.write(line);
                line = gtf.readLine();
                if (line != null) new_gtf.write(eol);
            }
            else
                line = gtf.readLine();
        }
        
        gtf.close();
        new_gtf.close();
        
    }
    public static void main(String[] args) throws IOException {
        CutNTinGtf c = new CutNTinGtf();
        c.cut();
    }
}
