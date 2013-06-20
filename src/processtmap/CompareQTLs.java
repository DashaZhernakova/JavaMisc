/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package processtmap;

//import eqtlmappingpipeline.pcaoptimum.eQTLFileCompare;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;
import umcg.genetica.io.text.TextFile;
import umcg.genetica.io.trityper.SNP;
import umcg.genetica.io.trityper.SNPLoader;
import umcg.genetica.io.trityper.TriTyperGenotypeData;
import umcg.genetica.io.trityper.util.DetermineLD;

/**
 *
 * @author dashazhernakova
 */

public class CompareQTLs {
    ArrayList<eQTL> different; //different pairs SNP-gene
    ArrayList<eQTL> shared; //same eQTLs shared by 2 lists
    ArrayList<eQTL> repeated; //repeated eQTLs
    ArrayList<eQTL> wrong_dir; //wrong direction eQTLs
    ArrayList<eQTL> eqtls1; //1st array to compare
    ArrayList<eQTL> eqtls2; //2nd array to compare
    public CompareQTLs(){
        different = new ArrayList<eQTL>();
        shared = new ArrayList<eQTL>();
        repeated = new ArrayList<eQTL>();
        wrong_dir = new ArrayList<eQTL>();
        eqtls1 = new ArrayList<eQTL>();
        eqtls2 = new ArrayList<eQTL>();        
    }
    public void cutShared() throws FileNotFoundException{
       for (eQTL eqtl : eqtls1)
            if (! repeated.contains(eqtl))
                different.add(eqtl);
    }
    /**
     * gets shared pairs SNP-gene, finds repeated and SNP-genes pairs with different directions
     */
    public void getShared(){
        System.out.println("2 " + eqtls2.size());
        System.out.println("1 " + eqtls1.size());
        for (eQTL eqtl2 : eqtls2){            
            for (eQTL eqtl1 : eqtls1)
                if ( (eqtl1.gene != null) && (eqtl2.gene != null)){
                if ((eqtl1.snp.equals(eqtl2.snp)) && (eqtl1.gene.equals(eqtl2.gene))){
                    if (checkDirection(eqtl1,eqtl2)){
                        if (shared.contains(eqtl1))
                            repeated.add(eqtl1);
                        else
                            shared.add(eqtl1);
                    }
                    else{
                        if (! wrong_dir.contains(eqtl1))
                            wrong_dir.add(eqtl1);
                        if (! wrong_dir.contains(eqtl2))
                            wrong_dir.add(eqtl2);
                    } }
                }}
    }
    public void getLinesFromExprTableByTrId(Set <String> absent_tr) throws IOException{
        TextFile e = new TextFile("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromTophat_coverageBed_Montgomery/expression_table.txt", false);
        Set<String> my_trs = new HashSet<String>();
        e.open();
        String line = e.readLine();
        while ( (line = e.readLine()) != null){
            String[] spl = line.split("\t");
            my_trs.add(spl[0]);
            if (absent_tr.contains(spl[0]))
                System.out.println(line);
        }
        for (String tr : absent_tr){
            if (! my_trs.contains(tr))
                System.out.println("absent in my expression table " + tr);
        }
        e.close();
    }
    public Set <String> getSharedGenes() throws IOException{
        Set <String> genes1 = new HashSet<String>();
        Set <String> shared_genes = new HashSet<String>();
        Set <String> absent_g = new HashSet<String>();
        GeneNameConverter conv = new GeneNameConverter();
        for (eQTL eqtl1 : eqtls1)
            genes1.add(eqtl1.gene);
        for (eQTL eqtl2 : eqtls2){
            if ( genes1.contains(eqtl2.gene) )
                shared_genes.add(eqtl2.gene);
            else{
                absent_g.add(eqtl2.whole_line.split("\t")[3]);
                //ArrayList<String> trs = conv.geneIdsToTrIds(eqtl2.gene);
                //if (trs != null )
                //absent_g.addAll(trs);
            }
           
        }
        
        //getLinesFromExprTableByTrId(absent_g);
        System.out.println("# shared genes: " + shared_genes.size());
        return shared_genes;
    }
    /**
     * gets shared pairs SNP-gene, finds repeated and SNP-genes pairs with different directions
     * adds shared pairs from BOTH files to shared
     */
        public void getBothSharedNoLD(){
            eQTLList l = new eQTLList();
            Set<String> genes = new HashSet<String>();
            for (eQTL eqtl2 : eqtls2)
            for (eQTL eqtl1 : eqtls1)
                if ( (eqtl1.gene != null) && (eqtl2.gene != null)){
                if ((eqtl1.snp.equals(eqtl2.snp)) && (eqtl1.gene.equals(eqtl2.gene))){
                    //if ( (checkDirection(eqtl1,eqtl2)) && (! l.containsSNPinLD(shared, eqtl1.gene, eqtl2.r_value)) && (! l.containsSNPinLD(shared, eqtl1.gene, eqtl2.r_value))){
                //if ( (eqtl1.gene.equals(eqtl2.gene)) && (eqtl1.gene.equals(eqtl2.gene)) ){    
                if ( (checkDirection(eqtl1,eqtl2)) ){
                        if (shared.contains(eqtl1))
                            repeated.add(eqtl1);
                        else{
                            shared.add(eqtl1);
                            genes.add(eqtl1.gene);
                            //System.out.println("1: " + eqtl1);
                            //shared.add(eqtl2);
                            //System.out.println("2: " + eqtl2);
                            break;
                        }
                        
                    }
                    else{
                        if (! wrong_dir.contains(eqtl1))
                            wrong_dir.add(eqtl1);
                        if (! wrong_dir.contains(eqtl2))
                            wrong_dir.add(eqtl2);
                    } 
                }}
                
            //ArrayList<eQTL> l = new eQTLList(shared).getWithoutLD();
            //for (eQTL e : l)
            //    System.out.println(e);   
        //String[] one = null;
        //String[] two = null;
        String one = "";
        String two = "";
            System.out.println("number of shared genes " + genes.size());
        /*for (int i = 0; i < shared.size()/2; i++){
            //one+= "," + Double.toString(java.lang.Math.pow(Float.valueOf(shared.get(2*i).whole_line.split("\t")[17]),2));
            //two+= "," + Double.toString(java.lang.Math.pow(Float.valueOf(shared.get(2*i + 1).whole_line.split("\t")[17]),2));
            one+= "," + Double.toString(java.lang.Math.pow(shared.get(2*i).r_value,2));
            two+= "," + Double.toString(java.lang.Math.pow(shared.get(2*i + 1).r_value,2));
            
        }
        System.out.println("my_mon<-c(" + one.replaceFirst(",", "") + ")");
        System.out.println("mon<-c(" + two.replaceFirst(",", "") + ")");
        */
        
            //getGenesNotPresentInThisList(0);
        //getNotPresentInThisList(2);
    }
    
    
    /**
     * gets gene names that are not present in eQTL list number num, but present in list 1 - num
     * @param num = 0 => gets gene names from eqtls2 that are not present in eqtls1; = 1 => gets gene names from eqtls1 that are not present in eqtls2;
     * @return list of gene names absent in the list
     */
    public Set<String> getGenesNotPresentInThisList(int num){
        ArrayList<eQTL> eqtls;
        ArrayList<eQTL> other;
        Set<String> absent_genes = new HashSet<String>();
        Set<String> genes = new HashSet<String>();
        if (num == 0){
            eqtls = eqtls1; 
            other = eqtls2;
        }
        else{
            eqtls = eqtls2;
            other = eqtls1;
        }
        for (eQTL ee : eqtls){
            genes.add(ee.gene);
        }
        int i = 0;
        for (eQTL e : other){
            if (! genes.contains(e.gene)){
                absent_genes.add(e.gene);
                System.out.println("absnt " + e.gene);
            }
            else
                i++;
        }
        System.out.println("shared genes " + i);
        return absent_genes;
    }
    /**
     * gets eQTLs that aren't present in list number num (that are unique for list number (1 - num) ) 
     * @param num = 0 => gets eQTLs from eqtls2 that are not present in eqtls1; = 1 => gets eQTLs from eqtls1 that are not present in eqtls2;
     * @return list of eQTLs absent in the list
     */
    public ArrayList<eQTL> getNotPresentInThisList(int num){
        ArrayList<eQTL> eqtls;
        ArrayList<eQTL> absent = new ArrayList<eQTL>();
        Set<String> genes = new HashSet<String>();
        
        if (num == 0){
            eqtls = eqtls1;
        for (eQTL ee : eqtls2){
            genes.add(ee.gene);
        }}
        else{
            eqtls = eqtls2;
            for (eQTL ee : eqtls1){
            genes.add(ee.gene);
        }
        }
        
        for (eQTL e : eqtls){
            //if ( (! shared.contains(e)) && (! genes.contains(e.gene))){
            if  (! shared.contains(e)){
                    absent.add(e);
                System.out.println("absent eQTL " + e);
            }
        }
        return absent;
    }
   
    /**
     * given 2 eQTLs checks whether the allelic direction is the same
     * @param eqtl1
     * @param eqtl2
     * @return true if the same direction
     */
    public boolean checkDirection(eQTL eqtl1, eQTL eqtl2){
        String[] alleles1 = eqtl1.snp_type.split("/");
        String[] alleles2 = eqtl2.snp_type.split("/");
        if ((eqtl1.direction.equals("*")) || (eqtl2.direction.equals("*")))
                return true;
        //if given alleles and directions are ok
        if ((eqtl1.direction.equals(eqtl2.direction)) && (eqtl1.allele.equals(eqtl2.allele)))
            return true;
        else
            //if given alleles are opposite and directions are opposite => same direction
            if ((! eqtl1.direction.equals(eqtl2.direction)) && (eqtl1.allele.equals(otherAllele(eqtl2.allele, alleles2))))
                return true;
            else{    
            //problems with directions
                for (String all : eqtl2.snp_type.split("/"))
                    if (! ((all.equals(alleles1[1])) || all.equals(alleles1[0])))
                        System.out.append("problems with snp types and directions");
                return false;
            }
    }
    
    /**
     * gets the other allele from list 
     * @param all E.g. A
     * @param alleles E.g.[A,G]
     * @return E.g. G
     */
    public String otherAllele(String all, String[] alleles){
        if (alleles[0].equals(all))
            return alleles[1];
        if (alleles[1].equals(all))
            return alleles[0];
        return null;
    }
    public String complement(String nucleotide){
        if (nucleotide.equals("C"))
            return "G";
        if (nucleotide.equals("G"))
            return "C";
        if (nucleotide.equals("T"))
            return "A";
        if (nucleotide.equals("A"))
            return "T";
        return null;
    }
    /**
     * prints list of eQTLs
     * @param list - list name
     * @param f_name - where to
     * @throws IOException 
     */
    public void printToFile(ArrayList<eQTL> list, String f_name) throws IOException{
        BufferedWriter out = new BufferedWriter(new FileWriter(f_name));
        for (eQTL e : list)
            out.write(e.whole_line + "\n");
        out.close();
    }
    /**
     * replaces repeated pairs SNP-gene
     * @param list - list in which to replace
     * @return - list with repeated pairs replaced
     */
    public ArrayList<eQTL> replaceRepeated(ArrayList<eQTL> list){
        ArrayList<eQTL> unique = new ArrayList<eQTL>();
        
        for (eQTL e : list){
            boolean rep = false;
            for (eQTL u_e : unique){
                if ((e.snp.equals(u_e.snp)) && (e.gene.equals(u_e.gene))){
                    rep = true;
                    break;
                }
            }
            if (! rep)
                unique.add(e);
        }
        return unique;
    }
    /**
     * reads eQTLs from 2 files into List<eQTL>
     * @param f_name1 - file 1
     * @param f_name2 - file 2
     * @throws FileNotFoundException
     * @throws IOException 
     */
    public void getFromFiles(String f_name1, String f_name2) throws FileNotFoundException, IOException{
        BufferedReader one = new BufferedReader(new FileReader(f_name1));
        BufferedReader two = null;
        GeneNameConverter converter = new GeneNameConverter();
        converter.transcrIdsToGeneNames();
        String line="";
        line = one.readLine();
        while ((line = one.readLine()) != null){
            String[] spl=line.split("\t");
            //for (String gene : spl[16].split("///"))
            //    eqtls1.add(new eQTL(spl[1],spl[8], spl[9], Float.toString(Math.signum(Float.valueOf(spl[10]).floatValue())), gene, line, Float.valueOf(spl[17])));
            //System.out.println("line " + line);
            //eqtls1.add(new eQTL(spl[1],spl[8], spl[9], Float.toString(Math.signum(Float.valueOf(spl[10]).floatValue())), converter.trIdsToGeneIds.get(spl[4]), line));
            eqtls1.add(new eQTL(spl[1],spl[8], spl[9], Float.toString(Math.signum(Float.valueOf(spl[10]).floatValue())), spl[16], line));
            //eqtls1.add(new eQTL(spl[1],spl[8], spl[9], Float.toString(Math.signum(Float.valueOf(spl[10]).floatValue())), converter.trIdsToGeneIds.get(spl[4]), line, Float.valueOf(spl[17])));
        }
        if (! f_name2.isEmpty()){
            two = new BufferedReader(new FileReader(f_name2));
            //if 2 different files, read 2nd file
            if (! f_name1.equals(f_name2)){
                line = two.readLine();
                while ((line = two.readLine()) != null){
                    String[] spl=line.split("\t");
                    for (String gene : spl[16].split("///"))
                        eqtls2.add(new eQTL(spl[1],spl[8], spl[9], Float.toString(Math.signum(Float.valueOf(spl[10]).floatValue())), gene, line, Float.valueOf(spl[17])));

                }
            }
            else
                eqtls2 = eqtls1;
        }
    }
    
    public ArrayList<String> getTranscriptsByGeneName(String annot_fname,Set<String> genes) throws IOException{
        //Map<String, ArrayList<String>> transcripts = new HashMap<String,ArrayList<String>>();
        ArrayList<String> transcripts = new ArrayList<String>();
        TextFile annot = new TextFile(annot_fname, false);
        annot.open();
        String[] fields;
        
        while ( (fields = annot.readLineElems(TextFile.tab)) != null){
            if ( genes.contains(fields[2]) )
                transcripts.add(fields[1]);
        }
                return transcripts;
    }
    
    public void getTranscriptExpression(String expr_fname, String out_fname, ArrayList<String> transcripts) throws IOException{
        TextFile expr = new TextFile(expr_fname, false);
        expr.open();
        TextFile out = new TextFile(out_fname, true);
        String line;
        while ( (line = expr.readLine()) != null){
            if (transcripts.contains(line.split("\t")[0]))
                out.write(line + "\n");
        }
        out.close();
    }
    public void getQTLsFromDerm(String fname) throws IOException{
        TextFile f = new TextFile(fname, false);
        f.open();
        //GeneNameConverter converter = new GeneNameConverter();
        //converter.makeIdToHUGOmap();
        String line= f.readLine();
        while ( (line = f.readLine()) != null ){
            String[] fields = line.split("\t");
            //System.out.println(fields[9]);
            if ( (fields[12].equals("0.001")) || (fields[12].equals("0.0001")) ) //only below perm thresh
                eqtls2.add(new eQTL(fields[0], "*", "*", "*", fields[2], line, Float.valueOf(fields[9])));
        }
        eQTLList l = new eQTLList(eqtls2);
        eqtls2 = l.filterMostSignificantPerGene();
        printToFile(eqtls2, "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromMontgomeryPaper_belowThr.txt");
        f.close();
    }
 
    
    public static void main(String[] args) throws IOException {
        
        //String path = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromTophat_coverageBed_Pickrell_yale/PCremove/Cis-10PCAsRemoved-GeneticVectorsNotRemoved";
        //String arg_path = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromTophat_coverageBed_Pickrell_argonne/Cis-10PCAsRemoved-GeneticVectorsNotRemoved";
        //String path = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/meta_out/hg19";
        //String arg_path = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Stranger-CEU/eQTLMapping-Uncorrected";
        //String path = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromTophat_coverageBed_Montgomery/Cis-5PCAsRemoved-GeneticVectorsNotRemoved";
        //String arg_path = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromTophat_intersectBed";
        
       
        CompareQTLs c = new CompareQTLs();
        
        
        //c.cutShared();
        //c.printToFile(c.replaceRepeated(c.eqtls1),path + "/probes_unique.txt");
        //c.printToFile(c.replaceRepeated(c.eqtls2),arg_path + "/probes_unique.txt");
        
        /*String path = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromTophat_coverageBed_exonwise_Montgomery_noSNPs/";
        c.getFromFiles("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromTophat_coverageBed_exonwise_Montgomery_noSNPs/Cis-5/eQTLProbesFDR0.05.txt",
                "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromTophat_coverageBed_exonwise_Montgomery/Cis-5PCremoved/eQTLsFDR0.05.txt");
                //"/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/meta_out/hg19/eQTLsFDR0.05_P+M_gene_ids.txt","");
                //"/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromTophat_coverageBed_Montgomery/meta/eQTLsFDR0.05.txt","");
        */
        //c.getQTLsFromDerm("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/RNASEQ60_CRG_INSERT_FULL-CIS1MB.observed.wperm.txt");
        //c.getQTLsFromDerm("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromMontgomeryPaper_mostSignificant_belowThr.txt");
        // "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Stranger/PCremove_ceu/eQTLProbesFDR0.05_allSNPs.txt"
        
        
        //c.getBothSharedNoLD();
        /*c.getShared();
        //c.printToFile(c.different, path + "/different_probes.txt");
        c.printToFile(c.wrong_dir, path + "diff_direction.txt");
        c.getSharedGenes();
        c.printToFile(c.shared, path + "shared.txt");
        System.out.println("number of shared " + c.shared.size());
        System.out.println("number of dif dir " + c.wrong_dir.size());
        */
        //System.out.println("number of significant in Montg " + c.eqtls2.size());
        //System.out.println("number of eqtls in my " + c.eqtls1.size());
        
       /*
       c.getTranscriptExpression("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromTophat_coverageBed_Montgomery/expression_table.txt"
                , "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromTophat_coverageBed_Montgomery/someTrExpr.txt", 
                c.getTranscriptsByGeneName("/Users/dashazhernakova/Documents/UMCG/hg19/annotation_transcr_hg19.txt", 
               c.getGenesNotPresentInThisList(0)) );
       */
        String mp = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Yale+Argonne/PCremove_newGenotypes/metaqtl_Montgomery_eqtls/eQTLProbesFDR0.05.txt",
               dp = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Yale+Argonne/PCremove_newGenotypes/metaqtl_deepSAGE_eqtls/eQTLProbesFDR0.05.txt",
                p = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Yale+Argonne/PCremove_newGenotypes/metaqtl/eQTLProbesFDR0.05.txt",
               pm = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromTophat_coverageBed_Montgomery/metaqtl_pickrell_eqtls/eQTLProbesFDR0.05.txt",
               dm = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromTophat_coverageBed_Montgomery/metaqtl_deepSAGE_eqtls/eQTLProbesFDR0.05.txt",
                m = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromTophat_coverageBed_Montgomery/meta/eQTLProbesFDR0.05.txt",
               pd = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/deepSAGE_transcr/meta/metaqtl_pickrell_eqtls/eQTLProbesFDR0.05.txt",
               md = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/deepSAGE_transcr/meta/metaqtl_Montgomery_eqtls/eQTLProbesFDR0.05.txt",
                d = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/deepSAGE_transcr/meta/eQTLProbesFDR0.05.txt";
        
        
    
        eQTLFileCompare co = new eQTLFileCompare();
        //co.compareOverlapAndZScoreDirectionTwoEQTLFiles(dm,d,
        //        "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/comparisons/dS2Mont-dS.txt");
        co.compareOverlapAndZScoreDirectionTwoEQTLFiles(
                "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/deepSAGE_transcr/meta/tags2transcr/xTagProbes_eqtls.txt",
                "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/deepSAGE_tag/metaqtl/eQTLsFDR0.05.txt",
                "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/deepSAGE_tag/metaqtl/comp_with_xTagProbes_pr.txt"
                );
        
        /*String t1 = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/tmp.txt", 
                t2 = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/tmp2.txt";
        co.compareOverlapAndZScoreDirectionTwoEQTLFiles(t1,t1, "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/tmp3.txt");
         * 
         */
         }
    /*
    public void testLD(String snp1, String snp2) throws IOException {
        TriTyperGenotypeData tgds = new TriTyperGenotypeData();
        tgds.load("/dir/where/genotypematrix.dat is");
        Integer id1 = tgds.getSnpToSNPId().get(snp1);
        Integer id2 = tgds.getSnpToSNPId().get(snp2);
        
        SNPLoader loader = tgds.createSNPLoader();
        
        DetermineLD ldcalc = new DetermineLD();
        
        if(id1 != null && id2 != null){
            SNP snpobj1 = tgds.getSNPObject(id1);
            SNP snpobj2 = tgds.getSNPObject(id2);
            
            loader.loadGenotypes(snpobj1);
            loader.loadGenotypes(snpobj2);
            
            double r2 = ldcalc.getRSquared(snpobj1, snpobj2, ldcalc.RETURN_R_SQUARED, ldcalc.INCLUDE_CASES_AND_CONTROLS, false);
        }
    }*/
}
