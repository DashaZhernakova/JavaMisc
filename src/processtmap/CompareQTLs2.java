package processtmap;

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

/**
 *
 * @author dashazhernakova
 */
public class CompareQTLs2 {
    
    ArrayList<eQTL> different; //different pairs SNP-gene
    ArrayList<eQTL> shared; //same eQTLs shared by 2 lists
    ArrayList<eQTL> repeated; //repeated eQTLs
    ArrayList<eQTL> wrong_dir; //wrong direction eQTLs
    ArrayList<eQTL> eqtls1; //1st array to compare
    ArrayList<eQTL> eqtls2; //2nd array to compare
    public CompareQTLs2(){
        different = new ArrayList<eQTL>();
        shared = new ArrayList<eQTL>();
        repeated = new ArrayList<eQTL>();
        wrong_dir = new ArrayList<eQTL>();
        eqtls1 = new ArrayList<eQTL>();
        eqtls2 = new ArrayList<eQTL>();        
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
            eqtls1.add(new eQTL(spl[1], spl[4], spl[8], spl[9], Float.toString(Math.signum(Float.valueOf(spl[10]).floatValue())), spl[16], line));
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
                        eqtls2.add(new eQTL(spl[1],spl[4],spl[8], spl[9], Float.toString(Math.signum(Float.valueOf(spl[10]).floatValue())), gene, line));

                }
            }
            else
                eqtls2 = eqtls1;
        }
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
                            else{
                                shared.add(eqtl1);
                                //shared.add(eqtl2);
                            }
                        }
                        else{
                            if (! wrong_dir.contains(eqtl1))
                                wrong_dir.add(eqtl1);
                            if (! wrong_dir.contains(eqtl2))
                                wrong_dir.add(eqtl2);
                        } 
                    }
                }
        }
    }
    public ArrayList<eQTL> getNotShared(int num)
    {
        
        ArrayList<eQTL> eqtls = new ArrayList<eQTL>();
        ArrayList<eQTL> notShared = new ArrayList<eQTL>();
        if (num == 0){
            eqtls = eqtls1; 
        }
        else{
            eqtls = eqtls2;
        }
        System.out.println("eqtls size " + eqtls.size());
        for (eQTL e : eqtls){
            if (! shared.contains(e)){
                notShared.add(e);
            }
        }
        return notShared;
    }
    
    public ArrayList<eQTL> getAbsentInNumFromAllEQTLS(int num, String eqtls_filen) throws IOException{
        TextFile all_eqtls_file = new TextFile(eqtls_filen, false);
        all_eqtls_file.open();
        ArrayList<eQTL> notShared = new ArrayList<eQTL>();
        ArrayList<eQTL> out = new ArrayList<eQTL>();
        Set<eQTL> all_eqtls = new HashSet<eQTL>();
        String line = all_eqtls_file.readLine();
        String[] spl;
        while ( (line = all_eqtls_file.readLine()) != null ){
            spl = line.split("\t");
            all_eqtls.add(new eQTL(spl[1],spl[8], spl[9], Float.toString(Math.signum(Float.valueOf(spl[10]).floatValue())), spl[16], line));
        }
        all_eqtls_file.close();
        
        notShared = getNotShared(Math.abs(num - 1));
        System.out.println("size of notShared " + notShared.size());
        System.out.println("size of unique notShared " + confineToUniqueSNPGenePairs(notShared).size());
        notShared = confineToUniqueSNPGenePairs(notShared);
        for (eQTL e : notShared){
            //System.out.println("e=" + e);
            for (eQTL ee : all_eqtls)
                if ( ( (e.gene.equals(ee.gene)) && (e.snp.equals(ee.snp)) ) ){
                    out.add(e);
                    break;
                }
            
                
                
        }
        System.out.println("out of those foung in eQTLs.txt " + out.size());
        return out;
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
    
    public ArrayList<eQTL> confineToUniqueSNPGenePairs(ArrayList<eQTL> list){
        ArrayList<eQTL> newList = new ArrayList<eQTL>();
        boolean contained;
        for (eQTL e : list){
            contained = false;
            for (eQTL ee : newList){
                if ( ( (e.gene.equals(ee.gene)) && (e.snp.equals(ee.snp)) ) ){
                    contained = true;
                    break;
                }
            }
            if (! contained)
                newList.add(e);
        }
        
        return newList;
    }
    public ArrayList<eQTL> confineToUniqueGenes(ArrayList<eQTL> list){
        ArrayList<eQTL> newList = new ArrayList<eQTL>();
        boolean contained;
        for (eQTL e : list){
            contained = false;
            for (eQTL ee : newList){
                if (e.gene.equals(ee.gene)){
                    contained = true;
                    break;
                }
            }
            if (! contained)
                newList.add(e);
        }
        
        return newList;
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
    
    public static void main(String[] args) throws IOException {
        
          CompareQTLs2 c = new CompareQTLs2();
        
        
        
        String path = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Peter-Bram_array_data/eQTLmapping/PCremove/Cis-5PCAsRemoved-GeneticVectorsNotRemoved/";
        c.getFromFiles("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Peter-Bram_array_data/eQTLmapping/PCremove/Cis-5PCAsRemoved-GeneticVectorsNotRemoved/eQTLsFDR0.05.txt",
                "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromTophat_coverageBed_exonwise_Montgomery/Cis-5PCremoved/eQTLsFDR0.05.txt");
           c.getShared();
       
           //c.printToFile(c.wrong_dir, path + "diff_direction.txt");
           
         
        //c.printToFile(c.confineToUniqueSNPGenePairs(c.shared), path + "shared_unique.txt");
        //c.printToFile(c.confineToUniqueSNPGenePairs(c.confineToUniqueSNPGenePairs(c.getNotShared(0))), path + "not_shared_unique.txt");
        
        
        //c.printToFile(c.getAbsentInNumFromAllEQTLS(0, "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/fromTophat_coverageBed_exonwise_Montgomery_noSNPs/Cis-5/eQTLs.txt"),
        //        path + "fromAllEQTLS.txt");
        
        System.out.println("number of unique SNP-gene pairs eqtls1 " + c.confineToUniqueSNPGenePairs(c.eqtls1).size());
        System.out.println("number of unique SNP-gene pairs eqtls2 " + c.confineToUniqueSNPGenePairs(c.eqtls2).size());
        System.out.println("number of unique genes eqtls1 " + c.confineToUniqueGenes(c.eqtls1).size());
        System.out.println("number of unique genes eqtls2 " + c.confineToUniqueGenes(c.eqtls2).size());
  
        System.out.println("number of shared " + c.shared.size());
        System.out.println("number of dif dir " + c.wrong_dir.size());
        System.out.println("number of not shared " + c.confineToUniqueSNPGenePairs(c.getNotShared(0)).size());
        System.out.println("number of unique shared " + c.confineToUniqueSNPGenePairs(c.shared).size());
        System.out.println("number of unique dif dir " + c.confineToUniqueSNPGenePairs(c.wrong_dir).size());
    }
}
