package edu.sysu.pmglab.ly.simulation.gwas;

import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;

public class GwasSummary {

    private static final Logger logger = LoggerFactory.getLogger(GwasSummary.class);

    public static double[][] performGWAS(int N, int snpsNum, double[] phenotypeScore, List<double[]> genotypes) {
        logger.info("开始进行 GWAS 分析 ......");
        double[][] gwasSummary = new double[snpsNum][4];
        for (int i = 0; i < snpsNum; i++) {

            // 构造自变量x，x需要是一个二维数组，形状为[N][1]
            double[][] x = new double[N][1];  // 每个SNP只有一个特征（基因型数据）
            for (int j = 0; j < N; j++) {
                x[j][0] = genotypes.get(j)[i];
            }

            // 线性回归
            OLSMultipleLinearRegression regression = new OLSMultipleLinearRegression();
            regression.newSampleData(phenotypeScore, x);

            double[] coefficients = regression.estimateRegressionParameters();
            double[] stdErrors = regression.estimateRegressionParametersStandardErrors();

            // 填充GWAS结果
            gwasSummary[i][0] = coefficients[1];  // BETA
            gwasSummary[i][1] = stdErrors[1];    // SE
            gwasSummary[i][2] = coefficients[1] / stdErrors[1];  // Zscore
            gwasSummary[i][3] = 2 * (1 - normalCDF(Math.abs(gwasSummary[i][2])));  // P-value
        }

        return gwasSummary;
    }

    // 标准正态分布的累积分布函数
    private static double normalCDF(double x) {
        return 0.5 * (1 + Erf.erf(x / Math.sqrt(2)));
    }
}
