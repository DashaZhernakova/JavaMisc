package misc;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author dashazhernakova
 */
public class eQTLs {
    Set<eQTL> eqtl_set;
    Map<String,Set<eQTL>> gene_eqlts_map;
    Map<String,Set<eQTL>> probe_eqlts_map;
    int size;
    int numUniqueGenes;
    int numUniqueProbes;
    int numUniqueSNPGenes;
    int numUniqueSNPProbes;
    String mode;
    
    public eQTLs(Set<eQTL> eqtls){
        eqtl_set = new HashSet<eQTL>();
        gene_eqlts_map = new HashMap<String, Set<eQTL>>();
        probe_eqlts_map = new HashMap<String, Set<eQTL>>();
        String gene = "", probe = "", snp = "";
        for (eQTL eqtl : eqtls){
            eqtl_set.add(eqtl);
            gene = eqtl.gene;
            probe = eqtl.probe;
            snp = eqtl.snp;
            Set<eQTL> set;
            boolean alreadyIn = false;
            //filling gene_eqlts_map
            if (gene_eqlts_map.containsKey(gene)){
                set = gene_eqlts_map.get(gene);
                for (eQTL e : set)
                    if (e.snp.equals(snp)){
                        alreadyIn = true;
                        break;
                    }
                if (! alreadyIn)
                    numUniqueSNPGenes++;
                set.add(eqtl);
            }
            else{
                set = new HashSet<eQTL>();
                set.add(eqtl);
                gene_eqlts_map.put(gene, set);
                numUniqueSNPGenes++;
            }
            
            alreadyIn = false;
            //filling probe_eqlts_map
            if (probe_eqlts_map.containsKey(probe)){
                set = probe_eqlts_map.get(probe);
                for (eQTL e : set)
                    if (e.snp.equals(snp)){
                        alreadyIn = true;
                        break;
                    }
                if (! alreadyIn)
                    numUniqueSNPProbes++;
                set.add(eqtl);
            }
            else{
                set = new HashSet<eQTL>();
                set.add(eqtl);
                probe_eqlts_map.put(probe, set);
                numUniqueSNPProbes++;
            }
            
        }
        size = eqtl_set.size();
        numUniqueGenes = gene_eqlts_map.keySet().size();
        numUniqueProbes = probe_eqlts_map.keySet().size();
        
    }
    public eQTLs(String f_name, boolean withHeader) throws IOException{
        eqtl_set = new HashSet<eQTL>();
        gene_eqlts_map = new HashMap<String, Set<eQTL>>();
        probe_eqlts_map = new HashMap<String, Set<eQTL>>();
         
        TextFile in = new TextFile(f_name, false);
        in.open();
        String line = "";
        if (withHeader)
            line = in.readLine();
        String gene = "", probe = "", snp = "";
        String[] spl;
        while ((line = in.readLine()) != null){
            spl = line.split("\t");
            gene = spl[16];
            probe = spl[4];
            snp = spl[1];
            
            
            eQTL eqtl = new eQTL(snp, probe, spl[8], spl[9], Float.toString(Math.signum(Float.valueOf(spl[10]).floatValue())), gene, line);
            eqtl_set.add(eqtl);
            
            
            Set<eQTL> set;
            boolean alreadyIn = false;
            //filling gene_eqlts_map
            if (gene_eqlts_map.containsKey(gene)){
                set = gene_eqlts_map.get(gene);
                for (eQTL e : set)
                    if (e.snp.equals(snp)){
                        alreadyIn = true;
                        break;
                    }
                if (! alreadyIn)
                    numUniqueSNPGenes++;
                set.add(eqtl);
            }
            else{
                set = new HashSet<eQTL>();
                set.add(eqtl);
                gene_eqlts_map.put(gene, set);
                numUniqueSNPGenes++;
            }
            
            alreadyIn = false;
            //filling probe_eqlts_map
            if (probe_eqlts_map.containsKey(probe)){
                set = probe_eqlts_map.get(probe);
                for (eQTL e : set)
                    if (e.snp.equals(snp)){
                        alreadyIn = true;
                        break;
                    }
                if (! alreadyIn)
                    numUniqueSNPProbes++;
                set.add(eqtl);
            }
            else{
                set = new HashSet<eQTL>();
                set.add(eqtl);
                probe_eqlts_map.put(probe, set);
                numUniqueSNPProbes++;
            }
            
        }
        size = eqtl_set.size();
        numUniqueGenes = gene_eqlts_map.keySet().size();
        numUniqueProbes = probe_eqlts_map.keySet().size();
        in.close();
    }
    
    @Override
    public String toString(){
        String out = "";
        for (eQTL e : eqtl_set){
            out += e + "; ";  
        }
        return out;
    }
    
    public String toStringAll(){
        String out = "";
        for (String g : gene_eqlts_map.keySet()){
            for (eQTL e : gene_eqlts_map.get(g))
                out += e.snp + " : " + e.probe + " : " + e.gene + " : " + e.snp_type + " : " + e.allele + " : " + e.direction + "\n";
        }
            
        return out;
    }
    public boolean containsGeneSNPPair(String gene, String snp){
        if (gene_eqlts_map.containsKey(gene)){
            for (eQTL e : gene_eqlts_map.get(gene)){
                if (e.snp.equals(snp))
                        return true;
            }
        }
            
        return false;
    }
    public boolean containsGene(String gene){
        if (gene_eqlts_map.containsKey(gene))
            return true;
        return false;
    }
    
    public boolean contains(eQTL eqtl){
        if (gene_eqlts_map.containsKey(eqtl.gene)){
            for (eQTL e : gene_eqlts_map.get(eqtl.gene)){
                if (e.whole_line.equals(eqtl.whole_line))
                        return true;
            }
        }
            
        return false;
    }
    
    public Set<eQTL> get(String gene, String snp){
        Set<eQTL> eqtls = new HashSet<eQTL>();
        if (gene_eqlts_map.containsKey(gene)){
            for (eQTL e : gene_eqlts_map.get(gene)){
                if (e.snp.equals(snp))
                    eqtls.add(e);
            }
        }
        
        return eqtls;
    }
    
    public Set<String> getSharedGenes(eQTLs eqtls2){
        Set<String> sharedGenes = new HashSet<String>();
        
        for (String gene : gene_eqlts_map.keySet())
            if (eqtls2.gene_eqlts_map.containsKey(gene))
                sharedGenes.add(gene);
        
        return sharedGenes;
    }
    public Set<String> getUniqueGenes(){
               
        return gene_eqlts_map.keySet();
    }
    public Set<String> getSharedProbes(eQTLs eqtls2){
        Set<String> sharedProbes = new HashSet<String>();
        
        for (String probe : probe_eqlts_map.keySet())
            if (eqtls2.probe_eqlts_map.containsKey(probe))
                sharedProbes.add(probe);
        
        return sharedProbes;
    }
    //old
    /*public eQTLcomparison getShared(eQTLs eqtls, boolean probeLevel, String mode_out){
        //Set<eQTL> shared = new HashSet<eQTL>();
        //Set<eQTL> wrong_dir = new HashSet<eQTL>();
        mode = mode_out;
        eQTLcomparison compare = new eQTLcomparison();
        if (probeLevel){
            for (Entry<String, Set<eQTL>> other_entry : eqtls.probe_eqlts_map.entrySet()){
                String other_name = other_entry.getKey();
                Set<eQTL> other_eqtls = other_entry.getValue();
                //check if this probe is present at all. If not - skip
                if (probe_eqlts_map.containsKey(other_name)){
                    Set<eQTL> this_eqtls = probe_eqlts_map.get(other_name);
                    //iterating only on eQTLs with this probe
                    compare.getSameSNPeQTLs(this_eqtls, other_eqtls, mode); 
                }
            }
        }
        else{ // if gene level:
            for (Entry<String, Set<eQTL>> other_entry : eqtls.gene_eqlts_map.entrySet()){
                String other_name = other_entry.getKey();
                Set<eQTL> other_eqtls = other_entry.getValue();
                //check if this probe is present at all. If not - skip
                if (gene_eqlts_map.containsKey(other_name)){
                    Set<eQTL> this_eqtls = gene_eqlts_map.get(other_name);
                    //iterating only on eQTLs with this probe
                    compare.getSameSNPeQTLs(this_eqtls, other_eqtls, mode);  
                }
            }
        }
        return compare;
    }*/
    public void printToFile(String f_name) throws IOException{
        TextFile out = new TextFile(f_name, true);
        out.open();
        for (eQTL e : eqtl_set)
            out.writeln(e.whole_line);
        out.close();
    }
    public void printToFileSorted(String f_name) throws IOException{
        TextFile out = new TextFile(f_name, true);
        out.open();
        ArrayList<eQTL> array;
        for (Entry < String, Set<eQTL>> e : gene_eqlts_map.entrySet()){
            array = new ArrayList<eQTL>();
            array.addAll(e.getValue());
            Collections.sort(array, new SNPComparator());
            for (eQTL eqtl : array){
                out.writeln(eqtl.whole_line);
            }
        }
        out.close();
    }
    
    public void outputNumbers(){
        System.out.println("Number of eQTLs overall: " + size);
        System.out.println("Number of unique genes: " + numUniqueGenes);
        System.out.println("Number of unique probes: " + numUniqueProbes);
        System.out.println("Number of unique pairs SNP-gene: " + numUniqueSNPGenes);
        System.out.println("Number of unique pairs SNP-probe: " + numUniqueSNPProbes);
    }
    public void outputComparisonNumbers(eQTLcomparison compare){
        System.out.println("\nCOMPARISON NUMBERS");
        System.out.println("Shared eQTLs: ");
        if (mode.equals("both"))
            System.out.println("Overall number of pairs of shared eQTLs: " + compare.numShared);
        else{
        System.out.println("Overall number of shared eQTLs: " + compare.numShared + " (out of " + size);
        
        System.out.println("eQTLs with wrong directions: ");
        }
    }
    public class SNPComparator implements Comparator<eQTL>{
        @Override
        public int compare(eQTL a, eQTL b) {
            return a.snp.compareTo(b.snp);
        }
    } 
    
    public static void main(String[] args) throws IOException {
        eQTLs eqtls1 = new eQTLs("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/tmp.txt", false);
        eQTLs eqtls2 = new eQTLs("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/tmp2.txt", false);
        System.out.println("eqtls1: ");
        eqtls1.outputNumbers();
        System.out.println("eqtls2: ");
        eqtls2.outputNumbers();
        /*
        eqtls1.getShared(eqtls2, true, "this").printSharedToFile("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/tmp1_out.txt");
         eqtls1.getShared(eqtls2, true, "this").outputNumbers();
        eqtls1.getShared(eqtls2, true, "other").printSharedToFile("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/tmp2_out.txt");
        eqtls1.getShared(eqtls2, true, "other").outputNumbers();
        eqtls1.getShared(eqtls2, true, "both").printSharedToFile("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/tmp_both_out.txt");
        eqtls1.getShared(eqtls2, true, "both").outputNumbers();
         * */
         
    }
}
