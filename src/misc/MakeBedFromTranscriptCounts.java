/*
 * To change this template+ "\t" + choose Tools | Templates
 * and open the template in the editor.
 */
package misc;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author dashazhernakova
 */
public class MakeBedFromTranscriptCounts{

    public void makeBedFromSimpleCounts(String inFname) throws FileNotFoundException, IOException{
        Map<String, String> transcr_counts = new HashMap<String, String>();

        TextFile counts = new TextFile(inFname, false);
        counts.open();
        //read transcript counts
        transcr_counts = counts.readAsHashMap(0, 1);
        
        TextFile annot = new TextFile("/Users/dashazhernakova/Documents/UMCG/hg19/annotation_transcr_hg19.txt", false);
        TextFile out = new TextFile(inFname.replace(".txt", ".bed"), true);
        
        //for every transcript from annotation from chr6 and non-zero count write a line with score = count
        annot.open();
        String[] elems;
        while ( (elems = annot.readLineElems(TextFile.tab)) != null){
            if (elems[3].equals("6")){
                String count = transcr_counts.get(elems[1]);
                if ( (count != null) && (! count.equals("0.0")) )
                    out.writelnTabDelimited(new String[] {elems[3], Integer.toString(Integer.valueOf(elems[4]) - 1), elems[5], elems[1], count});
            }
        }
        
        
        counts.close();
        annot.close();
        out.close();
    }
    public void makeGtfFromSimpleCounts(String inFname) throws FileNotFoundException, IOException{
        Map<String, String> transcr_counts = new HashMap<String, String>();

        TextFile counts = new TextFile(inFname, false);
        counts.open();
        //read transcript counts
        transcr_counts = counts.readAsHashMap(0, 1);
        
        TextFile annot = new TextFile("/Users/dashazhernakova/Documents/UMCG/hg19/chr6.gtf", false);
        TextFile out = new TextFile(inFname.replace(".txt", ".gtf"), true);
        
        //for every transcript from annotation from chr6 and non-zero count write a line with score = count
        annot.open();
        String line;
        while ( (line = annot.readLine()) != null){
            //if (elems[3].equals("6")){
                String count = transcr_counts.get(line.split("\t")[8].split("\"")[3]);
                if ( (count != null) && (! count.equals("0.0")) )
                    out.writeln(line.replace("\t*", " count " + count + ";\t*"));
                    //out.writelnTabDelimited(new String[] {elems[3], ".", "transcript", elems[4], elems[5], ".", ".", ".", "transcript_id \"" + elems[1] + "\"; " + "gene_id \"" + elems[2] + "\"; " + "count " + count + ";"});
            //}
        }
        
        
        counts.close();
        annot.close();
        out.close();
    }
    
    public void makeBedFromMISO(String inF) throws IOException{
        TextFile counts = new TextFile(inF, false);
                
        TextFile annot = new TextFile("/Users/dashazhernakova/Documents/UMCG/hg19/annotation_transcr_hg19.txt", false);
        TextFile out = new TextFile(inF + ".bed", true);
        
        //read transcript counts from MISO output
        Map<String, String> transcr_counts = new HashMap<String, String>();
        String[] count_els = counts.readLineElems(TextFile.tab);
        while ( (count_els = counts.readLineElems(TextFile.tab)) != null){
            List<String> ids = new ArrayList<String>();
            for (String id : count_els[4].split(",")){
                ids.add(id.split(":")[1]);
            }
            for (String count : count_els[6].split(",")){
                int num = Integer.valueOf(count.split(":")[0]);
                String c = count.split(":")[1];
                
                if (! c.equals("0"))
                    transcr_counts.put(ids.get(num), c);
            }            
        }
        
        //for every transcript from annotation from chr6 write a line with score = count
        annot.open();
        String[] elems;
        while ( (elems = annot.readLineElems(TextFile.tab)) != null){
            if (elems[3].equals("7")){
                String count = transcr_counts.get(elems[1]);
                if (count != null)
                    out.writelnTabDelimited(new String[] {elems[3], Integer.toString(Integer.valueOf(elems[4]) - 1), elems[5], elems[1], count});
            }
        }
        counts.close();
        annot.close();
        out.close();

    }
        public void makeBedFromIsoemFpkm() throws IOException{
        TextFile counts = new TextFile("/Users/dashazhernakova/Documents/UMCG/hapmapData/hg19/1382_1/6.filtered.sorted_by_name_isoforms_w_fpkm.gtf", false);
        counts.open();
        
        //TextFile annot = new TextFile("/Users/dashazhernakova/Documents/UMCG/hg19/annotation_transcr_hg19.txt", false);
        TextFile out = new TextFile("/Users/dashazhernakova/Documents/UMCG/hapmapData/hg19/1382_1/6.filtered.sorted_by_name_isoforms_w_fpkm.bed", true);
        
        //read transcript counts from MISO output
        Map<String, String> transcr_counts = new HashMap<String, String>();
        String[] els = counts.readLineElems(TextFile.tab);
        els = counts.readLineElems(TextFile.tab);
        Transcript cur_tr = new Transcript(els[0], els[3], els[4],els[5], els[6], els[8].split("\"")[1], els[8].split("\"")[3].replaceFirst(":.*", ""));
        //String cur_tr = "";
        while ( (els = counts.readLineElems(TextFile.tab)) != null){
            String tr_id = els[8].split("\"")[3].replaceFirst(":.*", "");
            if (! tr_id.equals(cur_tr.id)){
                out.writeln(cur_tr.toString());
                cur_tr = new Transcript(els[0], els[3],els[4],els[5], els[6], els[8].split("\"")[1], tr_id);
            
            }
            else
                cur_tr.end = els[4];
        }
       
        counts.close();
        out.close();

    }
    public void processInFolders(String dirName) throws IOException{
        File dir = new File(dirName);
        ArrayList<String> files = new ArrayList<String>();
        for (File ch : dir.listFiles())
            if (ch.isDirectory())
                for (File child : ch.listFiles()){
                    if (child.isDirectory())
                        if ( child.getName().contains("miso")){
                            files.add(child.getPath() + "/summary/miso.miso_summary");
                            makeBedFromMISO(child.getPath() + "/summary/miso.miso_summary");
                            break;
                        }
                }
        
    }
  
    public static void main(String[] args) throws FileNotFoundException, IOException{
        MakeBedFromTranscriptCounts m = new MakeBedFromTranscriptCounts();
        //m.makeGtfFromSimpleCounts("/Users/dashazhernakova/Documents/UMCG/pickrell_data/tophat_out/NA18522_yale/rsem_counts.txt");
        m.makeBedFromMISO("/Users/dashazhernakova/Documents/UMCG/flux_simulation/chr7/tophat_out/miso/summary/7.miso_summary");
        //m.makeBedFromIsoemFpkm();
        //m.processInFolders("/Users/dashazhernakova/Documents/UMCG/pickrell_data/tophat_out");
    }
}

class Transcript{
    String st;
    String end;
    String str;
    String chr;
    String fpkm;
    String gene;
    String id;
    Transcript(String c, String s, String e, String fp, String stra, String g, String i){
        st = s;
        end = e;
        str = stra;
        chr = c;
        fpkm = fp;
        gene = g;
        id = i;
    }
    @Override
    public String toString(){
        //return chr+ "\t" + "isoviz"+ "\t" + "transcript"+ "\t" + st+ "\t" + end+ "\t" + fpkm+ "\t" + str+ "\t" + "."+ "\t" + "gene_id \"" + gene + "\"; transcript_id \"" + id + "\";";
        return chr+ "\t" + st+ "\t" + end+ "\t" + id+ "\t" + fpkm;
    }
}
