/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package makedecoyfasta;

import java.util.HashMap;

/**
 *
 * @author dfermin
 */
public class Peptide {
    String seq;
    String decoySeq;
    double MHplus;
    int pepLen;
    int pepStart;
    int pepEnd;
    
    HashMap<Character, Double> AAmap = new HashMap<Character, Double>();
    
    
    Peptide(String s, int p) {
        seq = s;
        pepLen = s.length();
        pepStart = p;
        pepEnd = p + pepLen - 1;
        MHplus = 0;
        decoySeq = "";
        
        AAmap.put('A', 71.03711);
        AAmap.put('R', 156.10111);
        AAmap.put('N', 114.04293);
        AAmap.put('D', 115.02694);
        AAmap.put('C', 103.00919);
        AAmap.put('E', 129.04259);
        AAmap.put('Q', 128.05858);
        AAmap.put('G', 57.02146);
        AAmap.put('H', 137.05891);
        AAmap.put('I', 113.08406);
        AAmap.put('L', 113.08406);
        AAmap.put('K', 128.09496);
        AAmap.put('M', 131.04049);
        AAmap.put('F', 147.06841);
        AAmap.put('P', 97.05276);
        AAmap.put('S', 87.03203);
        AAmap.put('T', 101.04768);
        AAmap.put('W', 186.07931);
        AAmap.put('Y', 163.06333);
        AAmap.put('V', 99.06841);
        
        calcMass();
        makeReverseSeq();
    }
    
    
    private void calcMass() {
        double H = 1.007825; // proton
        MHplus = 18.010565 + H; // mass of water + proton
        String approvedAA = "ACDEFGHIKLMNPQRSTVWY";
        
        for(int i = 0; i < pepLen; i++) {
            char c = seq.charAt(i);
            
            // skip over non-standard amino acid characters
            if( !approvedAA.contains( Character.toString(c)) ) continue;
            
            MHplus += AAmap.get(c);
        }
        
        // Round the mass to 5 decimal places
        double tmp = round_dbl(MHplus, 5);
        MHplus = tmp;
    }
    
    private void makeReverseSeq() {
        String x = seq;
        
        if(seq.endsWith("K") || seq.endsWith("R"))  
            x = seq.substring(0, (seq.length() - 1) ); 
        
        StringBuilder sb = new StringBuilder(x);
        
        decoySeq = sb.reverse().toString();
        
        if(seq.endsWith("K") || seq.endsWith("R"))
            decoySeq += seq.substring( (seq.length()-1) );
    }
    
    
    // Function returns a double rounded to the given number of integers
    public double round_dbl(double value, int numPlaces) {
        double ret = 0;
        double N = Math.pow(10, numPlaces);

        ret = (double)Math.round( (value * N) ) / N;

        return ret;
    }
    
    
    public void printPeptide() {
        String line = pepStart + "-" + pepEnd + "\t" + 
                      seq + " -" + decoySeq + "\t" + MHplus + " Da"; 
        
        System.out.println(line);
    }
}
