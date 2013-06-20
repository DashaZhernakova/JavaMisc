/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package processtmap;

/**
 *
 * @author dashazhernakova
 */
public class eQTL {
    String snp;
    String gene;
    String probe;
    String snp_type;
    String allele;
    String direction;
    String whole_line;
    float r_value;
    
    eQTL(String new_snp, String new_type, String new_allele, String new_dir, String new_gene, String line){
        snp=new_snp;
        snp_type=new_type;
        allele=new_allele;
        direction = new_dir;
        gene=new_gene;
        whole_line=line;
    }
    eQTL(String new_snp, String new_type, String new_allele, String new_dir, String new_gene, String line, float r){
        snp=new_snp;
        snp_type=new_type;
        allele=new_allele;
        direction = new_dir;
        gene=new_gene;
        whole_line=line;
        r_value = r;
    }
    eQTL(String new_snp, String pr, String new_type, String new_allele, String new_dir, String new_gene, String line){
        snp=new_snp;
        probe = pr;
        snp_type=new_type;
        allele=new_allele;
        direction = new_dir;
        gene=new_gene;
        whole_line=line;
        //r_value = r;
    }
    @Override
    public String toString(){
        return snp + " : " + probe + " : " + gene;
    }
    
    /**
     * given 2 eQTLs checks whether the allelic direction is the same
     * @param eqtl1
     * @param eqtl2
     * @return true if the same direction
     */
    public boolean checkDirection(eQTL eqtl2){
        String[] alleles1 = snp_type.split("/");
        String[] alleles2 = eqtl2.snp_type.split("/");
        if ((direction.equals("*")) || (eqtl2.direction.equals("*")))
                return true;
        //if given alleles and directions are ok
        if ((direction.equals(eqtl2.direction)) && (allele.equals(eqtl2.allele)))
            return true;
        else
            //if given alleles are opposite and directions are opposite => same direction
            if ((! direction.equals(eqtl2.direction)) && (allele.equals(otherAllele(eqtl2.allele, alleles2))))
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
    
}
