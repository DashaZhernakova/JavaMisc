package misc;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.String;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
import umcg.genetica.io.text.TextFile;

/**
 * cuts some samples from expression table
 * @author dashazhernakova
 */
public class CutSamples {
    ArrayList<Integer> clmnToDel;
    ArrayList<Integer> clmnToLeave;
    ArrayList<String> mixups;
    public CutSamples(){
        clmnToDel = new ArrayList<Integer>();
        clmnToLeave = new ArrayList<Integer>();
    }
    /**
     * 
     * @param mix_ups - mixups ids (tab-separated) to delete
     */
    public CutSamples(String mix_ups){
        mixups = new ArrayList<String>();
        for (String mu : mix_ups.split("\t"))
            mixups.add(mu);
        clmnToDel = new ArrayList<Integer>();
        clmnToLeave = new ArrayList<Integer>();
    }
    /**
     * gets column numbers to delete, if column id matches pattern or id in mixups
     * @param line - table header
     * @param pattern - pattern of id to delete
     * @return - new header with deleted ids
     */
    public String getPos(String line, String pattern){
        String[] spl = line.split("\t");
        String new_line = spl[0]; 
        for (int i = 1; i < spl.length; i++){
            if ( mixups != null){
            if ((spl[i].matches(pattern)) || (mixups.contains(spl[i])))
                clmnToDel.add(i);
            else
                new_line += "\t" + spl[i];
            }
            else{
                if ((spl[i].matches(pattern)) )
                    clmnToDel.add(i);
            else
                new_line += "\t" + spl[i];
            }
                
        }
        return new_line;
    }
    
    public void setClmnToLeave(String header, String fname) throws IOException{
        TextFile f = new TextFile(fname, false);
        Set<String> idsToLeave = new HashSet<String>();
        String id = "";
        while ( (id = f.readLine()) != null ){
            idsToLeave.add(id);
        }
        
        String[] spl = header.split("\t");
        for (int i = 1; i < spl.length; i++){
            if (idsToLeave.contains(spl[i]))
                clmnToLeave.add(i);
        }
        
        
        f.close();
    }
    /**
     * gets repeated ids positions
     * @param header - table header
     * @return new header
     */
    public String cutRepeated(String header){
        String[] spl = header.split("\t");
        ArrayList<String> ids = new ArrayList<String>();
        String new_header = "";
        
        for (int i = 0; i < spl.length; i++){
            if (! ids.contains(spl[i])){
                ids.add(spl[i]);
                new_header += "\t" + spl[i];
            }
                
            else
                clmnToDel.add(i);
        }
        System.out.println(clmnToDel);
        //System.out.println(new_header.replaceFirst("\t", ""));
        
        return new_header.replaceFirst("\t", "");
    }
    
    /**
     * generates a list of random positions from [0:len]
     * @param numPosToLeave - number of positions to leave
     * @param len - overall number of initial positions
     */
    public void getRamdomPosToDel(int numPosToLeave, int len){
        ArrayList<Integer> pos = new ArrayList<Integer>();
        
        for (int i = 1; i < len; i++)
            pos.add(i);
        
        Collections.shuffle(pos);
        clmnToLeave.addAll(pos.subList(0, numPosToLeave));
    }
    
    /*
     * generates a list of random positions from ArrayList<Integer> pos
     */
    public void getRamdomPosToDel(int numPosToLeave, ArrayList<Integer> pos){

        Collections.shuffle(pos);
        System.out.println(pos.size());
        clmnToLeave.addAll(pos.subList(0, numPosToLeave));
     }
    
    /**
     * from an expression table header extracts random positions of ids that belong to the population specified in f_name file and are present in Individuals.txt
     * @param header - expression table header
     * @param population_f_name - name of a file containing ids that belong to a population in tab-separated format
     * @param numPosToLeave - number of positions to use in the result table 
     * @return
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public ArrayList<Integer> getGroupsPos(String header, String population_f_name, int numPosToLeave) throws FileNotFoundException, IOException{
        BufferedReader ids_f = new BufferedReader(new FileReader(population_f_name));
        BufferedReader individuals = new BufferedReader(new FileReader("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/HapMapPhase3GenotypeData/Individuals.txt"));
        ArrayList<Integer> pos = new ArrayList<Integer>();
        ArrayList<String> individs = new ArrayList<String>();
        String line = ids_f.readLine();
        ArrayList<String> ids = new ArrayList<String>();

        ids.addAll(Arrays.asList(line.split("\t")));
        System.out.println("ids.size()" + ids.size());
        //getting ids of individuals that are present in GenotypeMatrix
        while ( (line = individuals.readLine()) != null){
            individs.add(line);
        }
        String[] spl = header.replaceAll("GSM[0-9]+_","").split("\t");
        System.out.println("header len " + spl.length);
        //String[] spl = header.split("\t");
        for (int i = 0; i < spl.length; i++)
            if ( (ids.contains(spl[i])) && (individs.contains(spl[i])) )
                pos.add(i);
        
        //get random positions from chosen
        getRamdomPosToDel(numPosToLeave, pos);
        //System.out.println(pos.size());
        return pos;
    }
    
    /**
     * from expression table file generates 3 files for 3 populations with expression data for random individuals from that population
     * @param f_name - expression table file name
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public void cutRandomPosPopulations(String f_name) throws FileNotFoundException, IOException{
        BufferedReader in;
        BufferedWriter out;
        //files with population information (tab-separated ids)
        String[] files = {"/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/pedinfo2sample_ASIA.txt",
            "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/pedinfo2sample_CEU.txt", 
            "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/pedinfo2sample_YRI.txt"
        };
        //loop over populations
        for (int j = 0; j < files.length; j++){
            in = new BufferedReader(new FileReader(f_name));
            clmnToLeave = new ArrayList<Integer>();
            String line = in.readLine();
            String[] spl = line.split("\t");
            int len = spl.length;
            
            String f_n = files[j]; //current population file name
            out = new BufferedWriter(new FileWriter(f_name.replaceFirst(".txt", "_" + f_n.substring(f_n.length() - 7, f_n.length() - 4) + ".txt")));
            
            //gets random positions. last arg - number of columns to leave
            if (j == 0)
                getGroupsPos(line, f_n, 72);
            else
                getGroupsPos(line, f_n, 73);
            line = line.replaceAll("GSM[0-9]+_","");
            while ((line ) != null){
                spl = line.split("\t");
                String cur_line = spl[0]; 
                //for (int i = 1; i < 9; i++) //if first 9 columns about probe positions are present
                //    cur_line += "\t" + spl[i];

                for (int i = 1; i < len; i++)
                    if (clmnToLeave.contains(i))
                        
                        cur_line += "\t" + spl[i];
                out.write(cur_line + "\n");
               line = in.readLine();
            }
            out.close();
        }
    }
    
    /**
     * cuts columns from expression table specified in clmnToDel
     * @param f_name - expression table
     * @param pattern - pattern that should be matched to put the id in clmnToDel
     * @throws FileNotFoundException
     * @throws IOException 
     */
     
    public void cutPos(String f_name, String pattern) throws FileNotFoundException, IOException{
        BufferedReader in = new BufferedReader(new FileReader(f_name));
        BufferedWriter out = new BufferedWriter(new FileWriter(f_name.replaceFirst(".txt", "_cut.txt")));
        
        String line = in.readLine();
        String[] spl;
        
        //uncomment to put to clmnToDel positions of ids matching a pattern 
        out.write(getPos(line, pattern));
        
        //uncomment to put to clmnToDel positions of repeated ids
        //out.write(cutRepeated(line));
        
        while ((line = in.readLine()) != null){
            spl = line.split("\t");
            String cur_line = spl[0];
            for (int i = 1; i < spl.length; i++)
                if (! clmnToDel.contains(i))
                    cur_line += "\t" + spl[i];
            out.write("\n" + cur_line);
        }
        out.close();
    }
    /**
     * using the list of samples to delete (from file) deletes those samples
     */
    public void cutIdsFromFile(String original_file, String idsToDel_file) throws IOException{
        TextFile toDel = new TextFile(idsToDel_file, false);
        toDel.open();
        TextFile expression = new TextFile(original_file, false);
        expression.open();
        TextFile new_expr = new TextFile(original_file.replaceAll("WithDuplicates", "WithoutDuplicates"), true);
        new_expr.open();
        
        String line;
        Set<String> idsToDel = new HashSet<String>();
        //reading from file ids to delete
        while ( (line = toDel.readLine()) != null)
            idsToDel.add(line);
        
        String[] header = expression.readLineElems(TextFile.tab);
        //System.out.println(header.toString());
        String new_header = header[0];
        
        //to delete the special fields
        for (int i = 1; i < 9; i++)
            clmnToDel.add(i);
        //getting the positions of ids to delete listed in the file
        for (int i = 9; i < header.length; i++){
            if ( idsToDel.contains( header[i].substring(0, header[i].length() - 4) ) )
                clmnToDel.add(i);
            else //writing new header
                new_header += "\t" + header[i];
        }
        //System.out.println(clmnToDel);
        new_expr.writeln(new_header);
        //System.out.println(new_header);
        
        // rewriting the expression table skipping values at clmnToDel
        String[] spl;
        boolean first = true;
        while ((spl = expression.readLineElems(TextFile.tab)) != null){
            String cur_line = spl[0];
            for (int i = 1; i < spl.length; i++)
                if (! clmnToDel.contains(i))
                    cur_line += "\t" + spl[i];
            new_expr.write("\n" + cur_line);
            
        }
        
        toDel.close();
        expression.close();
        new_expr.close();
    }
    /**
     * using a list of samples to leave, takes only those samples from the table
     */
    public void leaveIdsFromFile(String fName, String samplesToLeaveFname) throws IOException{
        TextFile table = new TextFile(fName, false);
        TextFile new_table = new TextFile(fName.replace(".txt", "_new.txt"), true);
        
        setClmnToLeave(table.readLine(), samplesToLeaveFname);
        table.close();
        table.open();
        String[] els;
        while ( (els = table.readLineElems(TextFile.tab)) != null ){
            new_table.write(els[0]);
            for (int i = 1; i < els.length; i++){
                if (clmnToLeave.contains(i))
                    new_table.write("\t" + els[i]);
            }
            new_table.write("\n");
        }
        table.close();
        new_table.close();
    }
    
    /**
     * from 10th column averages expression values over 4 consecutive columns
     * @param f_name - expression table
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public void averagePos(String f_name) throws FileNotFoundException, IOException{
        BufferedReader in = new BufferedReader(new FileReader(f_name));
        BufferedWriter out = new BufferedWriter(new FileWriter(f_name.replaceFirst(".txt", "_avg.txt")));
        
        String line = in.readLine();
        String[] spl;
        spl = line.split("\t");
        out.write(spl[0]);
        //for (int i = 1; i < 9; i++)
        //        out.write("\t" + spl[i]);
        
        //write the header with the end of the id cut 
        String cur_id = "";
        for (int i = 1; i < spl.length - 3; i+= 4){
            if (i % 4 == 1){
               
                if (! (spl[i].substring(0, spl[i].length() - 4)).equals(spl[i+1].substring(0, spl[i+1].length() - 4)) &&
                        (spl[i+1].substring(0, spl[i+1].length() - 4)).equals(spl[i+2].substring(0, spl[i+2].length() - 4)) &&
                        (spl[i+2].substring(0, spl[i+2].length() - 4)).equals(spl[i+3].substring(0, spl[i+3].length() - 4))){
                    System.out.println("Problems with sample ids!!!");
                    return;
                }
            }
            out.write("\t" + spl[i].substring(0, spl[i].length() - 4));
            
        }
        while ((line = in.readLine()) != null){
            float avg = 0;
            spl = line.split("\t");
            String cur_line = spl[0]; 
            //for (int i = 1; i < 9; i++)
            //    cur_line += "\t" + spl[i];
            for (int i = 1; i < spl.length - 3; i+= 4){
                avg = (Float.valueOf(spl[i]) + Float.valueOf(spl[i+1]) + Float.valueOf(spl[i+2]) + Float.valueOf(spl[i+3]))/4;
                cur_line += "\t" + Float.toString(avg);
            }
            out.write("\n" + cur_line);
        }
        out.close();
    }
    public static void main(String[] args) throws FileNotFoundException, IOException{
        CutSamples cs = new CutSamples();
        //CutSamples cs = new CutSamples("NA18498");
        
        //cs.cutPos("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromTophat_coverageBed_Pickrell_yale/expression_table_full.txt", ".*_.+_yale");
        
        //cs.cutRandomPosPopulations("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Stranger/IlluminaExpressionDataOriginal_cut_cut.txt");
        
        //cs.cutRepeated("a	b	c	a	g	d	b	c	a	g	r	t	b"); 
        
        //cs.cutPos("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Stranger/IlluminaExpressionDataOriginal_cut.txt", "");
        //cs.cutPos("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Stranger/IlluminaExpressionDataOriginalWithDuplicates.txt", null);
        
        //System.out.println("GSM159536_NA12003	GSM159537_NA12004	GSM159538_NA12005	GSM159539_NA12006".replaceAll("GSM[0-9]+_","");
        
        //cs.averagePos("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Stranger/IlluminaExpressionDataOriginalWithDuplicates.txt");
        //if (args.length > 1)
          //  cs.cutPos(args[0],args[1]);
        /*TreeMap<String, ArrayList<String>> map = new TreeMap<String, ArrayList<String>>();
        ArrayList<String> tmp = new ArrayList<String>();
        tmp.add("aaa");
        map.put("1", tmp);
        
        tmp.add("bbb");
        map.get("1").add("ccc");
        */
        
        
        //cs.cutIdsFromFile("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Stranger/IlluminaExpressionDataOriginalWithDuplicates.txt", 
        //        "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Stranger/samplesToExculde.txt");
        
        //cs.averagePos("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Stranger/IlluminaExpressionDataOriginalWithoutDuplicates.txt");
        
        //cs.cutRandomPosPopulations("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Stranger/IlluminaExpressionDataOriginalWithoutDuplicates_avg.txt");
        
        //cs.cutPos("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromTophat_coverageBed_exonwise_Pickrell_argonne/expression_table.txt", "NA[0-9]*_[0-9].*");
        
        cs.leaveIdsFromFile("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Stranger/IlluminaExpressionDataOriginalWithoutDuplicates_avg_sampleIds.txt", 
                "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Stranger/CEU_samplesToInclude.txt");
    }
}
