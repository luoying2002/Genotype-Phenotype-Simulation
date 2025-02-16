package edu.sysu.pmglab.ly.simulation;

import edu.sysu.pmglab.ccf.viewer.CCFViewer;
import edu.sysu.pmglab.ccf.viewer.CCFViewerReader;
import edu.sysu.pmglab.ly.simulation.genotype.GenotypeSim;
import edu.sysu.pmglab.ly.simulation.phenotype.PhenotypeSim;

import java.io.IOException;

public class Test {
    public static void main(String[] args) throws IOException {

        String basePath = "/Users/yiguoshabi/Desktop/这是真研究/QuantTraitPredictor/resources/simulation/";

        String CCF_MAP = basePath + "map.ccf";
        String LD_chr1 = basePath + "LD_chr1.ccf";

        int CHR = 1;
        int N = 100;
        int snpsNum = 1000;
        double modifyRate = 0.01;
        double Hg2 = 0.6;
        double p = 0.1;

        GenotypeSim genotypeSim = new GenotypeSim(CCF_MAP, LD_chr1, CHR, N, snpsNum, modifyRate);

        PhenotypeSim phenotypeSim  = new PhenotypeSim(genotypeSim, Hg2, p);

        phenotypeSim.exportAll(basePath + "simulation");

        phenotypeSim.gwasSummary2csv(basePath + "simulation_gwas_summary.csv");

        new CCFViewer(new CCFViewerReader(phenotypeSim.genFiles.get("causal_effect")), 100, true);

    }
}
