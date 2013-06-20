
package misc;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

/**
 * trims reads in fastq file to reads f the length toLength
 * @author dashazhernakova
 */
public class ReadsTrimmer {
    int toLength;
    public ReadsTrimmer(int len){
        toLength = len;
    }
    public void trimRead(String fname) throws IOException{
        BufferedReader in = new BufferedReader(new FileReader(fname));
        BufferedWriter out = new BufferedWriter(new FileWriter(fname.replaceFirst(".fq", "_trimmed.fq")));
        
        String line = "";
        
        while ( (line = in.readLine()) != null){
            if (line.startsWith("@")){
                out.write(line + "\n");
                out.write(in.readLine().substring(1, toLength - 1) + "\n");
                out.write(in.readLine() + "\n");
                out.write(in.readLine().substring(1, toLength - 1) + "\n");
            }
        }
        in.close();
        out.close();
    }
    
    public static void main(String[] args) throws IOException {
        ReadsTrimmer rt = new ReadsTrimmer(25);
        if (args.length > 0)
            rt.trimRead(args[0]);
        else
            System.out.println("provide an input file please");
        //rt.trimRead("/Users/dashazhernakova/Documents/UMCG/ribosomal/111208_SN163_0449_BD0HDVACXX_L8_GTAGAG.fq");
    }
}
