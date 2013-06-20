
package processtmap;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.TreeSet;
/**
 * gets unaligned reads from sam and fastq files
 * @author dashazhernakova
 */
public class GetUnalignedReads {
    /**
     * gets unaligned reads from filtered! sam and fastq files
     * @param fastq - fastq file name
     * @param sam - filtered! sam file name
     * @param outp - output file name
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public void getUnaligned(String fastq, String sam, String outp) throws FileNotFoundException, IOException{
        BufferedReader fq = new BufferedReader(new FileReader(fastq));
        BufferedReader aln = new BufferedReader(new FileReader(sam));
        BufferedWriter out = new BufferedWriter(new FileWriter(outp));
        
        TreeSet<String> aln_ids = new TreeSet<String>();
        String line ="";
        
        //all aligned reads from filtered sam file
        while ( (line = aln.readLine()) != null){
            if (! line.startsWith("@"))
                aln_ids.add(line.split("\t")[0]);
        }
        //writing to the output reads from fastq absent in filtered sam file
        String cur_read = "";
        while ( (line = fq.readLine()) != null){
            if (line.startsWith("@")){
                if (! cur_read.isEmpty())
                    out.write(cur_read);
            cur_read = "";
            if (! aln_ids.contains(line.subSequence(1, line.length() - 2)))
                cur_read += line + "\n";
            }
            else if (! cur_read.isEmpty())
                cur_read += line + "\n";
        }
        out.write(cur_read);
        out.close();
    }
    
    public static void main(String[] args) throws IOException {
        GetUnalignedReads gur = new GetUnalignedReads();
        if (args.length > 3)
            gur.getUnaligned(args[1], args[2], args[3]);
        else
            System.out.println("not enough arguments");
        //gur.getUnaligned("/Users/dashazhernakova/Documents/UMCG/ribosomal/111208_SN163_0449_BD0HDVACXX_L8_GTAGAG_trimmed.fq", "/Users/dashazhernakova/Documents/UMCG/ribosomal/accepted_hits.filtered.sam", "/Users/dashazhernakova/Documents/UMCG/ribosomal/unaligned_java.fq");
    }
}
