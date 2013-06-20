package processtmap;

//import eqtlmappingpipeline.pcaoptimum.eQTLFileCompare;
import java.io.IOException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.Set;
import umcg.genetica.io.text.TextFile;

/**
 *
 * @author dashazhernakova
 */
public class eQTLcomparison {
    //Set<eQTL> shared;
    //Set<eQTL> wrong_dir;
    eQTLs shared;
    eQTLs shared_probeLevel;
    eQTLs shared_geneLevel;
    eQTLs notShared_probeLevel;
    eQTLs notShared_geneLevel;
    eQTLs wrong_dir;
    eQTLs eqtls1;
    eQTLs eqtls2;
    int numShared;
    int numWrongDir;
    int numGeneSNPShared;
    int numProbeSNPShared;
    int numGeneSNPWrongDir;
    int numProbeSNPWrongDir; 
    int numGenesShared;
    int numProbesShared;
    String mode;
    boolean probeLevel;
    int overlap;
    public eQTLcomparison(String f_name1, String f_name2) throws IOException{
        //shared = new HashSet<eQTL>();
        //wrong_dir = new HashSet<eQTL>();
        eqtls1 = new eQTLs(f_name1, true);
        eqtls2 = new eQTLs(f_name2, true);
        
        System.out.println("first list ( " + f_name1 + "): ");
        eqtls1.outputNumbers();
        System.out.println("\nsecond list( " + f_name2 + "): ");
        eqtls2.outputNumbers();
        overlap = 0;
        
    }
    public eQTLcomparison(){}
     /**
     * from 2 sets of eqtls gets the ones that have the same snp name
     * @param this_eqtls
     * @param other_eqtls
     * @param mode - this, other or both - write eQTLs from this_eqtls, from other_eqtls or both
     */
    public Set<eQTL> getSameSNPeQTLs(Set<eQTL> eqtl_set1, Set<eQTL> eqtl_set2){
        Set<eQTL> shared = new HashSet<eQTL>();
        Set<eQTL> wrong_d = new HashSet<eQTL>();
        for (eQTL eqtl2 : eqtl_set2){
        for (eQTL eqtl1 : eqtl_set1){
            //if the same SNP
            if (eqtl1.snp.equals(eqtl2.snp)){
                //if the same allelic direction
                if (eqtl1.checkDirection(eqtl2)){
                   //writing one eQTL or the pair
                   if (mode.equals("this"))
                      shared.add(eqtl1);
                   else if (mode.equals("other"))
                      shared.add(eqtl2);
                   else if (mode.equals("both")){
                      shared.add(eqtl1);
                      shared.add(eqtl2);
                   }
                   else
                       System.out.println("wrong mode: " + mode);
                }
                else{ //different allelic direction
                      if (mode.equals("this"))
                         wrong_d.add(eqtl1);
                      else if (mode.equals("other"))
                         wrong_d.add(eqtl2);
                      else if (mode.equals("both")){
                          wrong_d.add(eqtl1);
                          wrong_d.add(eqtl2);
                      }
                      
                } 
            }
        }
        }
        return shared;
        
    }
    /**
     * deprecated
     * @param probe
     * @param mode_out 
     */
    public void getShared(boolean probe, String mode_out){
        //Set<eQTL> shared = new HashSet<eQTL>();
        //Set<eQTL> wrong_dir = new HashSet<eQTL>();
        mode = mode_out;
        //eQTLcomparison compare = new eQTLcomparison();
        probeLevel = probe;
        if (probeLevel){
            for (Entry<String, Set<eQTL>> other_entry : eqtls2.probe_eqlts_map.entrySet()){
                String other_name = other_entry.getKey();
                Set<eQTL> other_eqtls = other_entry.getValue();
                //check if this probe is present at all. If not - skip
                if (eqtls1.probe_eqlts_map.containsKey(other_name)){
                    Set<eQTL> this_eqtls = eqtls1.probe_eqlts_map.get(other_name);
                    //iterating only on eQTLs with this probe
                    shared_probeLevel = new eQTLs(getSameSNPeQTLs(this_eqtls, other_eqtls)); 
                }
            }
        }
        else{ // if gene level:
            for (Entry<String, Set<eQTL>> other_entry : eqtls2.gene_eqlts_map.entrySet()){
                String other_name = other_entry.getKey();
                Set<eQTL> other_eqtls = other_entry.getValue();
                //check if this probe is present at all. If not - skip
                if (eqtls1.gene_eqlts_map.containsKey(other_name)){
                    Set<eQTL> this_eqtls = eqtls1.gene_eqlts_map.get(other_name);
                    //iterating only on eQTLs with this probe
                    shared_geneLevel = new eQTLs(getSameSNPeQTLs(this_eqtls, other_eqtls));   
                }
            }
        }
        
    }
    public Set<eQTL> getNotShared(String mode_out){
        Set<eQTL> not_shared = new HashSet<eQTL>();
        
        if (mode.equals("this")){
            for (eQTL e : eqtls1.eqtl_set){
                if (! shared.contains(e))
                    not_shared.add(e);
            }
        }
        
        else if (mode.equals("other")){
            for (eQTL e : eqtls2.eqtl_set){
                if (! shared.contains(e))
                    not_shared.add(e);
            }
        }
        return not_shared;
    }
    
    public Set<eQTL> getNotSharedByGene(String mode_out){
        Set<eQTL> not_shared = new HashSet<eQTL>();
        
        if (mode.equals("this")){
            for (eQTL e : eqtls1.eqtl_set){
                if (! eqtls2.containsGene(e.gene))
                    not_shared.add(e);
            }
        }
        
        else if (mode.equals("other")){
            for (eQTL e : eqtls2.eqtl_set){
                if (! eqtls1.containsGene(e.gene))
                    not_shared.add(e);
            }
        }
        notShared_geneLevel = new eQTLs(not_shared);
        return not_shared;
    }
    
    public void getShared2(boolean probe, String mode_out){
        System.out.println("Getting shared");
        Set<eQTL> shared_eqtls = new HashSet<eQTL>();
        Set<eQTL> wrong_d = new HashSet<eQTL>();
        Set<String> g = new HashSet<String>();
        Set<String> s = new HashSet<String>();
        Set<String> w = new HashSet<String>();
        
        int i=0, j=0;
        mode = mode_out;
        probeLevel = probe;
        String geneOrProbe1, geneOrProbe2;
        for (eQTL eqtl1 : eqtls1.eqtl_set){
        for (eQTL eqtl2 : eqtls2.eqtl_set){
            if (probeLevel){
                geneOrProbe1 = eqtl1.probe;
                geneOrProbe2 = eqtl2.probe;
            }
            else{
                geneOrProbe1 = eqtl1.gene;
                geneOrProbe2 = eqtl2.gene;
            }
            
            //TEMP!!! TO REMOVE!!!
            
            //if (! eqtl1.probe.equals(eqtl2.probe)){
            //if the same SNP
            if ( (eqtl1.snp.equals(eqtl2.snp)) && (geneOrProbe1.equals(geneOrProbe2)) ){
                g.add(geneOrProbe1);
                
                overlap++;
                //if the same allelic direction
                if (eqtl1.checkDirection(eqtl2)){
                   //writing one eQTL or the pair
                    i++;
                    s.add(geneOrProbe1);
                   if (mode.equals("this"))
                      shared_eqtls.add(eqtl1);
                   else if (mode.equals("other"))
                      shared_eqtls.add(eqtl2);
                   else if (mode.equals("both")){
                      shared_eqtls.add(eqtl1);
                      shared_eqtls.add(eqtl2);
                   }
                   else
                       System.out.println("wrong mode: " + mode);
                }
                else{ //different allelic direction
                    j++;
                    w.add(geneOrProbe1);
                      if (mode.equals("this"))
                         wrong_d.add(eqtl1);
                      else if (mode.equals("other"))
                         wrong_d.add(eqtl2);
                      else if (mode.equals("both")){
                          wrong_d.add(eqtl1);
                          wrong_d.add(eqtl2);
                      }
                      
                //} 
            }
        }
        }
        }
        shared = new eQTLs(shared_eqtls);
        if (probeLevel)
            shared_probeLevel = new eQTLs(shared_eqtls);
        else
            shared_geneLevel = new eQTLs(shared_eqtls);
        wrong_dir = new eQTLs(wrong_d);
        
        System.out.println("Same dirs: " + i + "\nWrong dirs: " + j);
        System.out.println("Overall " + g.size() + "\nShared " + s.size() + "\nwrong " + w.size());
    }
    
    public void makeTableForScatterPlot(String fname) throws IOException{
        eQTLs eqtls = null;
        TextFile plot = new TextFile(fname, true);
        //choosing the eQTLs set that is NOT in shared
        if (mode.equals("this"))
            eqtls = eqtls2;
        else if (mode.equals("other"))
            eqtls = eqtls1;
        else
            System.out.println("wrong mode");
        String gene, snp;
        Set<eQTL> corresp_eqtls = null;
        for (Entry<String, Set<eQTL>> entry : shared.gene_eqlts_map.entrySet()){
            gene = entry.getKey();
            
            for (eQTL e : entry.getValue()){
                snp = e.snp;
                corresp_eqtls = eqtls.get(gene, snp);
                for (eQTL eq : corresp_eqtls){
                    plot.writeln(gene + "\t" + snp + "\t" + e.whole_line.split("\t")[10] + "\t" + eq.whole_line.split("\t")[10]);
                }
                
            }
        }
        plot.close();
    }
    public void printSharedToFile(String level, String f_name) throws IOException{
        TextFile out = new TextFile(f_name, true);
        out.open();
        if (level.equals("probe"))
        for (eQTL e : shared_probeLevel.eqtl_set)
            out.writeln(e.whole_line);
        else if (level.equals("gene"))
        for (eQTL e : shared_geneLevel.eqtl_set)
            out.writeln(e.whole_line);
        else
            System.out.println("wrong level value");
        out.close();
    }
    public void printToFile(Set<eQTL> eqtls, String out_fname) throws IOException{
        TextFile out = new TextFile(out_fname, true);
        for (eQTL e : eqtls)
            out.writeln(e.whole_line);
        out.close();
    }
    
    
    
    public eQTLs getWrongDir(eQTLs eqtls, String mode, boolean probeLevel){
        
        Set<eQTL> wrong_d = new HashSet<eQTL>();
        //mode = mode_out;
        //probeLevel = probe;
        Set<eQTL> val;
        String geneOrProbe1, geneOrProbe2;
        if (probeLevel){
            for (Entry<String, Set<eQTL>> e : eqtls.probe_eqlts_map.entrySet()){
                val = e.getValue();
                for (eQTL eqtl1 : val){
                    for (eQTL eqtl2 : val){
                        geneOrProbe1 = eqtl1.probe;
                        geneOrProbe2 = eqtl2.probe;
                    
                        //if the same SNP
                        if ( (eqtl1.snp.equals(eqtl2.snp)) && (geneOrProbe1.equals(geneOrProbe2)) ){
                            //if the same allelic direction
                            if (! eqtl1.checkDirection(eqtl2)){

                                //different allelic direction
                                  if (mode.equals("this"))
                                     wrong_d.add(eqtl1);
                                  else if (mode.equals("other"))
                                     wrong_d.add(eqtl2);
                                  else if (mode.equals("both")){
                                      wrong_d.add(eqtl1);
                                      wrong_d.add(eqtl2);
                                  }

                            } 
                        }
                    }
                }
            }
        }
        else{
            for (Iterator<Entry<String, Set<eQTL>>> it = eqtls.gene_eqlts_map.entrySet().iterator(); it.hasNext();) {
                Entry<String, Set<eQTL>> e = it.next();
                val = e.getValue();
                for (eQTL eqtl1 : val){
                    for (eQTL eqtl2 : val){
                    
                        geneOrProbe1 = eqtl1.gene;
                        geneOrProbe2 = eqtl2.gene;
                    
                        //if the same SNP
                        if ( (eqtl1.snp.equals(eqtl2.snp)) && (geneOrProbe1.equals(geneOrProbe2)) ){
                            //if the same allelic direction
                            if (! eqtl1.checkDirection(eqtl2)){

                                //different allelic direction
                                  if (mode.equals("this"))
                                     wrong_d.add(eqtl1);
                                  else if (mode.equals("other"))
                                     wrong_d.add(eqtl2);
                                  else if (mode.equals("both")){
                                      wrong_d.add(eqtl1);
                                      wrong_d.add(eqtl2);
                                  }

                            } 
                        }
                    }
                }
            }
        }
        
        
        return new eQTLs(wrong_d);
    }
    
    public eQTLs getWrongDir(){
        
        Set<eQTL> wrong_d = new HashSet<eQTL>();
        //mode = mode_out;
        //probeLevel = probe;
        Set<eQTL> val;
        String geneOrProbe1, geneOrProbe2;
        
        for (Iterator<Entry<String, Set<eQTL>>> it = eqtls1.gene_eqlts_map.entrySet().iterator(); it.hasNext();) {
                Entry<String, Set<eQTL>> e = it.next();
                val = e.getValue();
                if (e.getKey().equals("GSTP1")){
                    System.out.println("GSTP1");
                }
                if (eqtls2.gene_eqlts_map.containsKey(e.getKey())){
                    for (eQTL eqtl1 : val){
                        for (eQTL eqtl2 : eqtls2.gene_eqlts_map.get(e.getKey())){

                            geneOrProbe1 = eqtl1.gene;
                            geneOrProbe2 = eqtl2.gene;

                            //if the same SNP
                            if ( (eqtl1.snp.equals(eqtl2.snp)) && (geneOrProbe1.equals(geneOrProbe2)) ){
                                //if the same allelic direction
                                if (! eqtl1.checkDirection(eqtl2)){

                                    //different allelic direction
                                      
                                          wrong_d.add(eqtl1);
                                          wrong_d.add(eqtl2);
                                      

                                } }
                            }
                        }
                }
            }
       
        
        
        return new eQTLs(wrong_d);
    }
    
    public void outputNumbers(){
        System.out.println("\n\nWRONG ALLELIC DIRECTIONS");
        System.out.println("\neQTLs with wrong allelic direction in the first list: ");
        eQTLs wr1 = getWrongDir(eqtls1, "both", false);
        System.out.println(wr1.toStringAll());
        wr1.outputNumbers();
        
        System.out.println("\neQTLs with wrong allelic direction in the second list: ");
        eQTLs wr2 = getWrongDir(eqtls2, "both", false);
        System.out.println(wr2.toStringAll());
        wr2.outputNumbers();
        
        System.out.println("\n\n-------------------------------------------\n"
                + "COMPARISON NUMBERS");
        /*eQTLs shared_eqtls = new eQTLs(shared);
        numGeneSNPShared = shared_eqtls.numUniqueSNPGenes;
        numProbeSNPShared = shared_eqtls.numUniqueSNPProbes;
        numGenesShared = shared_eqtls.numUniqueGenes;
        numProbesShared = shared_eqtls.numUniqueProbes;
        eQTLs wrong_dir_eqtls = new eQTLs(wrong_dir);
        numGeneSNPWrongDir = wrong_dir_eqtls.numUniqueSNPGenes;
        numProbeSNPWrongDir = wrong_dir_eqtls.numUniqueSNPProbes;*/
        System.out.println("\nSHARED eQTLs: ");
        
        //System.out.println("Number of unique shared probes: " + eqtls1.getSharedProbes(eqtls2).size());
        //System.out.println("Number of unique shared genes: " + eqtls1.getSharedGenes(eqtls2).size() + "\n\n");
        
        if (probeLevel) System.out.println("Number of unique shared probes: " + shared.probe_eqlts_map.keySet().size());
        else System.out.println("Number of unique shared genes: " + shared.gene_eqlts_map.keySet().size() + "\n\n");
        
        
        if (probeLevel){
            System.out.println("-------------PROBE level comparison-------");
        System.out.println("Overall number of shared eQTLs:   " + shared_probeLevel.size + "   out of   " + 
                eqtls1.size +"   in this eQTL list and   " + eqtls2.size + "   in the other");
        
        System.out.println("Number of shared unique pairs SNP - probe: " + 
                shared_probeLevel.numUniqueSNPProbes + " out of " + eqtls1.numUniqueSNPProbes + " (" + 
                String.format("%.1f",100*(float)shared_probeLevel.numUniqueSNPProbes/eqtls1.numUniqueSNPProbes) + "%) in this eQTL list and " + 
                eqtls2.numUniqueSNPProbes +  " (" + 
                String.format("%.1f",100*(float)shared_probeLevel.numUniqueSNPProbes/eqtls2.numUniqueSNPProbes) + "%) in the other");
        
       // System.out.println("Number of shared unique pairs SNP - gene: " + shared_probeLevel.numUniqueSNPGenes + " out of " + eqtls1.numUniqueSNPGenes + " in this eQTL list and " + eqtls2.numUniqueSNPGenes + " in the other");
        //    System.out.println("shared: " + shared_probeLevel);
        
        }
        if (! probeLevel){
            System.out.println("-------------GENE level comparison-------");
        System.out.println("Overall number of shared eQTLs: " + shared_geneLevel.size + " out of " + eqtls1.size + " (" + shared_geneLevel.size/eqtls1.size + ") in this eQTL list and " + eqtls2.size + " in the other");
        
        //System.out.println("Number of shared unique pairs SNP - probe: " + shared_geneLevel.numUniqueSNPProbes + " out of " + eqtls1.numUniqueSNPProbes + " in this eQTL list and " + eqtls2.numUniqueSNPProbes + " in the other");
        
        System.out.println("Number of shared unique pairs SNP - gene:   " + 
                shared_geneLevel.numUniqueSNPGenes + "  out of  " + eqtls1.numUniqueSNPGenes + " ("+
                String.format("%.1f",100*(float)shared_geneLevel.numUniqueSNPGenes/eqtls1.numUniqueSNPGenes) + "%)     in this eQTL list and  " + 
                eqtls2.numUniqueSNPGenes + " (" +
                String.format("%.1f",100*(float)shared_geneLevel.numUniqueSNPGenes/eqtls2.numUniqueSNPGenes) + "%)    in the other");
        //System.out.println("shared: " + shared_geneLevel);
        
        }
       // eQTLs wrong_dir_eqtls = new eQTLs(wrong_dir);
        System.out.println("\neQTLs with wrong directions: ");
        wrong_dir.outputNumbers();
        /*for (eQTL eq : wrong_dir.eqtl_set){
            System.out.println(eq);
        }*/
        Set<eQTL> diff = new HashSet<eQTL>();
        int different = 0;
        for (eQTL e : wrong_dir.eqtl_set){
            if (! ((wr1.containsGeneSNPPair(e.gene, e.snp)) || (wr2.containsGeneSNPPair(e.gene, e.snp))) ){
                different++;
                diff.add(e);
                //System.out.println(e);
            }
            //else
            //    System.out.println("absent!!" + e.whole_line);
        }
        if (different > 0){
            System.out.println("wrong directions between lists: \n" + diff + " \n" + different + " out of " + Integer.toString(shared.probe_eqlts_map.keySet().size()+diff.size()) + "(" + Float.toString(100 - 100*(float)diff.size()/(shared.probe_eqlts_map.keySet().size()+diff.size())) +"%)");
            System.out.println("wrong directions between lists: \n" + diff + " \n" + different + " out of " + overlap + "(" + Float.toString(100 - 100*(float)diff.size()/overlap) +"%)");
        }else
            System.out.println("!All wrong directions are within the lists, not between them!");
        
        //System.out.println("wrong direction: " + wrong_dir);
       
        
    }

    public static void main(String[] args) throws IOException {
        
        String p_pr = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Y+A_noNorm/metaqtl_15/eQTLsFDR0.05_PvalThresholdGeneLevel.txt",
                m_pr = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/Montgomery_noNorm/metaqtl_10/eQTLProbesFDR0.05-GeneLevel.txt",
                d="/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/deepSAGE_transcr_noNorm/metaqtl_20/eQTLProbesFDR0.05-ProbeLevel.txt";
        String path = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/new/";
        String dS = path + "deepSAGE_transcr/noNorm/metaqtl_15/eQTLProbesFDR0.05-GeneLevel.txt",
p=path + "Pickrell/noNorm/metaqtl_15/eQTLProbesFDR0.05-GeneLevel.txt",
m=path + "Montgomery/noNorm/metaqtl_10/eQTLProbesFDR0.05-GeneLevel.txt",
p_dS=path + "deepSAGE_transcr/noNorm/metaqtl_15/Pickrell2deepSAGE/eQTLProbesFDR0.05-GeneLevel.txt",
m_dS=path + "deepSAGE_transcr/noNorm/metaqtl_15/Montgomery2deepSAGE/eQTLProbesFDR0.05-GeneLevel.txt",
m_p=path + "Pickrell/noNorm/metaqtl_15/Montgomery2Pickrell/eQTLProbesFDR0.05-GeneLevel.txt",
dS_p=path + "Pickrell/noNorm/metaqtl_15/deepSAGE2Pickrell/eQTLProbesFDR0.05-GeneLevel.txt",
dS_m=path + "Montgomery/noNorm/metaqtl_10/deepSAGE2Montgomery/eQTLProbesFDR0.05-GeneLevel.txt",
p_m=path + "Montgomery/noNorm/metaqtl_10/Pickrell2Montgomery/eQTLProbesFDR0.05-GeneLevel.txt",
d_gen = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/new/deepSAGE_transcr/genes/metaqtl_15/eQTLProbesFDR0.05-GeneLevel.txt";
        
        String pa = path + "Pickrell/array/metaqtl_10/eQTLProbesFDR0.05-GeneLevel.txt",
                ma = path + "Montgomery/array/metaqtl_10/eQTLProbesFDR0.05-GeneLevel.txt",
                da = path + "deepSAGE_transcr/array/metaqtl_10/eQTLProbesFDR0.05-GeneLevel.txt",
                ppa = path + "Pickrell/array/metaqtl_10/rnaseq_subset/eQTLProbesFDR0.05-GeneLevel.txt",
                mma = path + "Montgomery/array/metaqtl_10/rnaseq_subset/eQTLProbesFDR0.05-GeneLevel.txt",
                dda = path + "deepSAGE_transcr/array/metaqtl_10/rnaseq_subset/eQTLProbesFDR0.05-GeneLevel.txt";
        String subs = "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/new/deepSAGE_tag/randomSubsets/55samples/metaqtl_10_1/eQTLsFDR0.05_PvalThreshold-ProbeLevel.txt";
        
        eQTLcomparison comp = new eQTLcomparison("/Users/dashazhernakova/Documents/UMCG/data/geuvadis/expression_table/FIN/eQTLmapping/dna-seq/metaqtl_10/rna-seq_subset_impute_cr0.5_pr0.8_info0.3/eQTLProbesFDR0.05-ProbeLevel.txt",
                "/Users/dashazhernakova/Documents/UMCG/data/geuvadis/expression_table/FIN/eQTLmapping/rna-seq/imputed/metaqtl_10_CR0.5_pr0.8_info0.3/eQTLProbesFDR0.05-ProbeLevel.txt");
        
        
        
        
        //eQTLcomparison comp = new eQTLcomparison(dda, path + "deepSAGE_transcr/array/metaqtl_10/rnaseq_subset/eQTLs.txt");
        /*eQTLcomparison comp = new eQTLcomparison("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/new/deepSAGE_transcr/genes/metaqtl_15/array_subset/notShared.txt",
                "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/new/deepSAGE_transcr/array/metaqtl_10/eQTLProbesFDR0.05-GeneLevel.txt");
        */
        //eQTLcomparison comp = new eQTLcomparison("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/new/deepSAGE_transcr/genes/metaqtl_15/array_subset/notReplicated_in_tagwise/metaqtl_15/eQTLProbesFDR0.05-GeneLevel.txt",
        //        "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/new/deepSAGE_transcr/array/metaqtl_10/eQTLProbesFDR0.05-GeneLevel.txt");
        comp.getShared2(false, "this");
        comp.outputNumbers();
        /*eQTLs nsh = new eQTLs(comp.getNotShared("other"));
        System.out.println(nsh.size);
       //comp.shared.printToFile("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/tmp.txt");
        nsh.printToFile("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/tmp.txt");
        */
        /*
        eQTLcomparison comp = new eQTLcomparison(
                "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/new/deepSAGE_tag/metaqtl_15/eQTLsFDR0.05_PvalThreshold-ProbeLevel.txt", 
                "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/new/deepSAGE_tag/metaqtl_0/eQTLProbesFDR0.05-ProbeLevel.txt");
        comp.getShared2(true, "both");
        */
        //System.out.println("\n" + comp.overlap);
        //comp.shared.printToFile("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/new/deepSAGE_tag/metaqtl_15/tmp_all_eqtls/shared.txt");
        //comp.wrong_dir.printToFile("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/new/deepSAGE_tag/metaqtl_15/PermutedEQTLsPermutationRound1_wrong_dir.txt");
        //comp.outputNumbers();
        //comp.getWrongDir().printToFileSorted("/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/new/deepSAGE_tag/metaqtl_15/tmp_all_eqtls/oppositeDirections.txt");
        
         /*
        eQTLFileCompare co = new eQTLFileCompare();
        co.compareOverlapAndZScoreDirectionTwoEQTLFiles(
                p_pr,m_pr,
              "/Users/dashazhernakova/Documents/UMCG/GeneticalGenomicsDatasets/tmp.txt");
      
        */
        
        
    }
}
