package edu.sysu.pmglab.ly.simulation.phenotype;

import edu.sysu.pmglab.bytecode.ImmutableByteCode;
import edu.sysu.pmglab.ccf.CCFMetaItem;
import edu.sysu.pmglab.ccf.CCFWriter;
import edu.sysu.pmglab.ccf.record.IRecord;
import edu.sysu.pmglab.ccf.type.IFieldType;
import edu.sysu.pmglab.ly.simulation.genotype.GenotypeSim;
import edu.sysu.pmglab.ly.simulation.gwas.GwasSummary;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

public class PhenotypeSim {

    private static final Logger logger = LoggerFactory.getLogger(PhenotypeSim.class);
    private static final Random rand = new Random();

    public int snpsNum;        // snp数量
    public int N;              // 样本数量
    public double Hg2;         // 遗传率（遗传力）
    public double p;           // 只有约 p 的SNP标记会被赋予非零效应（即对表型有影响）
    public GenotypeSim genotypeSim;    // 基因型数据对象
    public List<double[]> standardizeGenotypes;     // 标准化基因型数据
    public double[] causalEffect;      // 每个 snp 的因果效应值
    public double[] phenotypeScore;    // 样本表型值
    public double[][] gwasSummary;    // GWAS分析结果
    public Map<String, String> genFiles;    // 模拟生成文件路径

    public PhenotypeSim(GenotypeSim genotypeSim, double Hg2, double p) {
        this.snpsNum = genotypeSim.standardizeGenotypes.get(0).length;
        this.N = genotypeSim.standardizeGenotypes.size();
        this.Hg2 = Hg2;
        this.p = p;
        this.genotypeSim = genotypeSim;
        this.standardizeGenotypes = genotypeSim.standardizeGenotypes;
        this.causalEffect = new double[snpsNum];
        this.phenotypeScore = new double[N];
        this.genFiles = new HashMap<>(4);  // 依次为 _phenotype_standardize.ccf ; _phenotype.ccf ; _gwas_summary.ccf") ; _causal_effect.ccf ;

        // 计算每一个 snp 的因果效应值
        generateCausalEffects();

        // 计算每个样本的表型值
        generatePhenotypeScores();

        // GWAS分析
        this.gwasSummary = GwasSummary.performGWAS(N, snpsNum, phenotypeScore, standardizeGenotypes);  // BETA, SE, Zscore, P;
    }

    // 生成因果效应
    private void generateCausalEffects() {
        logger.info("开始生成因果效应值 ......");
        // 计算正态分布的标准差
        double variance = Hg2 / (snpsNum * p);

        // 遍历每个SNP，生成因果效应
        for (int i = 0; i < snpsNum; i++) {
            // p比例的SNP有效
            if (rand.nextDouble() < p) {
                // 有效的SNP，从正态分布中生成因果效应
                causalEffect[i] = rand.nextGaussian() * Math.sqrt(variance);
            } else {
                // 无效的SNP，因果效应为0
                causalEffect[i] = 0.0;
            }
        }
        logger.info("因果效应生成完毕！");
    }

    // 生成表型值
    private void generatePhenotypeScores() {
        logger.info("开始生成表型值 ......");
        // 计算噪声项的标准差
        double noiseVariance = 1 - Hg2;

        // 遍历每个样本，计算表型值
        for (int i = 0; i < N; i++) {
            double phenotype = 0.0;

            // 计算每个样本的表型值，累加所有SNP的效应
            for (int j = 0; j < snpsNum; j++) {
                phenotype += standardizeGenotypes.get(i)[j] * causalEffect[j];
            }

            // 加入噪声项
            phenotype += rand.nextGaussian() * Math.sqrt(noiseVariance);

            // 存储表型值
            phenotypeScore[i] = phenotype;
        }
        logger.info("表型值生成完毕！");
    }

    public void phenotype2ccf(String ccfPath) throws IOException {
        logger.info("正在导出表型模拟结果至 ccf 文件 ......");
        CCFWriter writer = new CCFWriter(new File(ccfPath));
        writer.addField("quantity", IFieldType.get("float32"));
        for(int i=0; i<this.snpsNum; i++){
            writer.addField("snp" + i, IFieldType.get("int8"));
        }
        IRecord record = writer.getRecord();
        int count;
        for(int i=0; i<this.N; i++){
            count = 0;
            record.getBox("quantity").char2Object(new ImmutableByteCode(String.valueOf(this.phenotypeScore[i])));
            for(int snp: this.genotypeSim.haplotypes.get(i)){
                record.getBox("snp" + count).char2Object(new ImmutableByteCode(String.valueOf(snp)));
                count += 1;
            }
            writer.write(record);
        }
        writer.close();
        logger.info("ccf 文件已保存至：{}", ccfPath);
    }

    public void phenotypeStandardize2ccf(String ccfPath) throws IOException {
        logger.info("正在导出表型模拟结果（标准化后）至 ccf 文件 ......");
        CCFWriter writer = new CCFWriter(new File(ccfPath));
        writer.addField("quantity", IFieldType.get("float32"));
        for(int i=0; i<this.snpsNum; i++){
            writer.addField("snp" + i, IFieldType.get("float32"));
        }
        IRecord record = writer.getRecord();
        int count;
        for(int i=0; i<this.N; i++){
            count = 0;
            record.getBox("quantity").char2Object(new ImmutableByteCode(String.valueOf(this.phenotypeScore[i])));
            for(double snp: this.standardizeGenotypes.get(i)){
                record.getBox("snp" + count).char2Object(new ImmutableByteCode(String.valueOf(snp)));
                count += 1;
            }
            writer.write(record);
        }
        writer.close();
        logger.info("ccf 文件（标准化后）已保存至：{}", ccfPath);
    }

    public void gwasSummary2ccf(String ccfPath) throws IOException {
        logger.info("正在导出 GWAS 分析结果至 ccf 文件 ......");
        CCFWriter writer = new CCFWriter(new File(ccfPath))
                .addField("VALUE::SNP", IFieldType.get("string"))
                .addField("VALUE::BETA", IFieldType.get("float32"))
                .addField("VALUE::SE", IFieldType.get("float32"))
                .addField("VALUE::Zscore", IFieldType.get("float32"))
                .addField("VALUE::P", IFieldType.get("float32"));

        IRecord record = writer.getRecord();
        int count = 0;
        for(double[] i: gwasSummary){
            record.getBox("SNP").char2Object(new ImmutableByteCode("snp" + count));
            record.getBox("BETA").char2Object(new ImmutableByteCode(String.valueOf(i[0])));
            record.getBox("SE").char2Object(new ImmutableByteCode(String.valueOf(i[1])));
            record.getBox("Zscore").char2Object(new ImmutableByteCode(String.valueOf(i[2])));
            record.getBox("P").char2Object(new ImmutableByteCode(String.valueOf(i[3])));
            writer.write(record);
            count++;
        }
        writer.close();
        logger.info("GWAS 分析的 ccf 文件已保存至：{}", ccfPath);
    }

    public void gwasSummary2csv(String csvPath) throws IOException {
        logger.info("正在导出 GWAS 分析结果至 csv 文件 ......");

        try (BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(Files.newOutputStream(Paths.get(csvPath)), StandardCharsets.UTF_8))) {
            // 写入文件头部
            writer.write("SNP,BETA,SE,Zscore,P\n");

            // 写入每一行数据
            for (int i = 0; i < snpsNum; i++) {
                // SNP 名称为 snp1, snp2, ...
                writer.write("snp" + (i + 1) + "," + gwasSummary[i][0] + "," + gwasSummary[i][1] + "," + gwasSummary[i][2] + "," + gwasSummary[i][3] + "\n");
            }
        } catch (IOException e) {
            logger.error(e.getMessage());
        }

        logger.info("GWAS 分析的 csv 文件已保存至：{}", csvPath);
    }

    public void causalEffect2ccf(String ccfPath) throws IOException {
        logger.info("正在导出 因果效应值 至 ccf 文件 ......");

        CCFMetaItem meta1 = CCFMetaItem.of("Hg2", String.valueOf(this.Hg2));
        CCFMetaItem meta2 = CCFMetaItem.of("p", String.valueOf(this.p));

        CCFWriter writer = new CCFWriter(new File(ccfPath))
                .addField("VALUE::SNP", IFieldType.get("string"))
                .addField("VALUE::causalEffect", IFieldType.get("float32"))
                .addMeta(meta1)
                .addMeta(meta2);

        IRecord record = writer.getRecord();
        int count = 0;
        for(double i: causalEffect){
            record.getBox("SNP").char2Object(new ImmutableByteCode("snp" + count));
            record.getBox("causalEffect").char2Object(new ImmutableByteCode(String.valueOf(i)));
            writer.write(record);
            count++;
        }
        writer.close();
        logger.info("因果效应值的 ccf 文件已保存至：{}", ccfPath);
    }

    public void exportAll(String pathAndPreixName) throws IOException {
        logger.info("正在导出至 ccf 文件 ......");

        genFiles.put("phenotype_standardize", pathAndPreixName + "_phenotype_standardize.ccf");
        genFiles.put("phenotype", pathAndPreixName + "_phenotype.ccf");
        genFiles.put("gwas_summary", pathAndPreixName + "_gwas_summary.ccf");
        genFiles.put("causal_effect", pathAndPreixName + "_causal_effect.ccf");

        phenotypeStandardize2ccf(genFiles.get("phenotype_standardize"));
        phenotype2ccf(genFiles.get("phenotype"));
        gwasSummary2ccf(genFiles.get("gwas_summary"));
        causalEffect2ccf(genFiles.get("causal_effect"));

    }

}
