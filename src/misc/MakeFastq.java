/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.TreeSet;

/**
 *
 * @author dashazhernakova
 */
public class MakeFastq {
    public void make(String fname) throws IOException{
        BufferedReader in = new BufferedReader(new FileReader(fname));
        BufferedWriter out = new BufferedWriter(new FileWriter(fname.replaceFirst(".txt", ".fastq")));
        TreeSet<String> ids = new TreeSet<String>();
        String line = "";
        while ( (line = in.readLine()) != null){
            String[] spl = line.split("\t");
            if (! ids.contains(spl[0])){
            out.write("@" + spl[0] + "\n");
            String seq = spl[1].replaceAll("N+$", "").replaceAll("^N+", "");
            out.write(seq + "\n");
            out.write("+" + spl[0] + "\n");
            //char[] qual = new char[spl[1].length()];
            //Arrays.fill(qual, 'I');
            for (int i= 0; i < seq.length(); i++)
                out.write("I");
            out.write("\n");
            ids.add(spl[0]);
            }
        }
        out.close();
        in.close();
    }
    public static void main(String[] args) throws IOException {
        MakeFastq mf = new MakeFastq();
        //System.out.print("nnaaaasnnnnssdffgnnnnnnnnnnn".replaceAll("n+$", "").replaceAll("^n+", ""));
        mf.make("/Users/dashazhernakova/Documents/UMCG/forLude/GPL6400_2005-03-16_HG17_WG_CGH.txt");
    }
}
