package edu.sysu.pmglab.ly.simulation.genotype;

import edu.sysu.pmglab.bytecode.ImmutableByteCode;
import edu.sysu.pmglab.ccf.CCFMeta;
import edu.sysu.pmglab.ccf.CCFReader;
import edu.sysu.pmglab.ccf.CCFWriter;
import edu.sysu.pmglab.ccf.filter.CCFFilter;
import edu.sysu.pmglab.ccf.record.IRecord;
import edu.sysu.pmglab.ccf.type.IFieldType;
import org.apache.commons.math3.linear.*;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.stream.IntStream;

public class GenotypeSim {

    private static final Logger logger = LoggerFactory.getLogger(GenotypeSim.class);
    private static final Random rand = new Random();

    public int N; // 样本量
    public int snpsNum; // SNP数量
    public double[][] C; // LD子矩阵
    public float[] P; // 等位基因频率
    public double[][] COV; // 协方差矩阵
    public List<double[]> standardizeGenotypes; // 模拟生成的标准化基因型数据
    public List<int[]> haplotypes; // 模拟生成的基因型数据（0，1，2）
    public int startSim;        // 模拟起始位置

    public GenotypeSim(String mapLdref, String ldref, int chr, int n, int snpsNum, double modifyRate) throws IOException {
        this.N = n;
        this.snpsNum = snpsNum;
        this.haplotypes = new ArrayList<>();
        this.standardizeGenotypes = new ArrayList<>();

        // 初始化等位基因频率 P
        this.initializeP(mapLdref, chr);

        // 初始化LD矩阵 C
        this.initializeC(ldref);

        // 计算协方差矩阵 COV
        this.calculateCovarianceMatrix();

        // 生成基因型数据
        this.generateAndStandardizeGenotypes(modifyRate);
        this.validateGenotypeFrequencies();
    }

    // 读取等位基因频率 P
    private void initializeP(String mapLdref, int chr) throws IOException {
        CCFReader reader = new CCFReader(mapLdref);
        CCFFilter ccfFilter = new CCFFilter(reader.getFile());
        ccfFilter.addRecordFilter(record -> ((Byte) record.getBox("chr").get()) == chr);

        List<IRecord> records = new ArrayList<>();
        long pointer;
        while ((pointer = ccfFilter.filter()) != -1L) {
            reader.seek(pointer);
            records.add(reader.read());
        }

        this.startSim = rand.nextInt(records.size() - snpsNum);

        logger.info("正在读取等位基因频率 P ......");
        this.P = new float[this.snpsNum];
        for (int i = this.startSim; i < this.startSim + this.snpsNum; i++) {
            this.P[i - this.startSim] = (float) records.get(i).getBox("af_UKBB").get();
        }
    }

    // 初始化LD矩阵 C
    private void initializeC(String ldref) throws IOException {
        CCFReader reader = new CCFReader(ldref);

        CCFMeta ccfMeta = reader.getMetas();

        // 提取元信息
        String scale = ccfMeta.get("scale").get(0).toString().split("=")[1];
        String[] scaleValues = scale.split(",");
        int i_length = Integer.parseInt(scaleValues[0]);
        int p_length = Integer.parseInt(scaleValues[1]);
        int x_length = Integer.parseInt(scaleValues[2]);

        int[] rowIndex = new int[i_length];
        int[] colPtr = new int[p_length];
        double[] values = new double[x_length];
        int count = 0;

        for (IRecord record : reader) {
            if (count < i_length) rowIndex[count] = (int) record.getBox("i").get();
            if (count < p_length) colPtr[count] = (int) record.getBox("p").get();
            if (count < x_length) values[count] = (float) record.getBox("x").get();
            count++;
        }

        reader.close();

        this.C = new double[snpsNum][snpsNum];

        // 填充 C 矩阵
        logger.info("正在初始化LD矩阵 C ......");
        int endSim = this.startSim + this.snpsNum - 1;
        for (int col = this.startSim; col <= endSim; col++) {
            int startIdx = colPtr[col] - 1;   // 列开始位置
            int endIdx = colPtr[col + 1] - 1; // 列结束位置

            // 遍历该列的非零元素
            for (int idx = startIdx; idx < endIdx; idx++) {
                int row = rowIndex[idx] - 1;  // 获取行索引（转换为 0-based）
                double value = values[idx];   // 获取非零值

                // 仅填充子矩阵内的元素
                if (row >= this.startSim && row <= endSim && col >= this.startSim && col <= endSim) {
                    int rowInSubMatrix = row - startSim;    // 计算在子矩阵中的行索引
                    int colInSubMatrix = col - startSim;    // 计算在子矩阵中的列索引
                    this.C[rowInSubMatrix][colInSubMatrix] = value;  // 填充值
                }
            }
        }
    }

    // 计算协方差矩阵 COV
    public void calculateCovarianceMatrix() {
        Double[] matrix = new Double[snpsNum * snpsNum]; // 线性化矩阵

        logger.info("正在计算协方差矩阵 COV ......");

        // 使用 parallelStream 并行计算协方差矩阵
        IntStream.range(0, snpsNum).parallel().forEach(i -> {
            IntStream.range(i, snpsNum).forEach(j -> {
                double cov = C[i][j] * Math.sqrt(P[i]*(1-P[i])*P[j]*(1-P[j]));
                matrix[i * snpsNum + j] = cov;  // 上三角
                if (i != j) {
                    matrix[j * snpsNum + i] = cov;  // 下三角
                }
            });
        });

        // 将一维数组转换为二维数组
        double[][] matrix2D = new double[snpsNum][snpsNum];
        for (int i = 0; i < snpsNum; i++) {
            for (int j = 0; j < snpsNum; j++) {
                matrix2D[i][j] = matrix[i * snpsNum + j];
            }
        }

        RealMatrix realMatrix = new Array2DRowRealMatrix(matrix2D);

        logger.info("正在检查协方差矩阵是否正定 ......");
        // 检查矩阵是否正定
        if (!checkPositiveDefinite(realMatrix)) {
            logger.warn("协方差矩阵不正定，修正中 ......");
            realMatrix = makePositiveDefinite(realMatrix);
        } else {
            logger.info("协方差矩阵正定.");
        }

        this.COV = realMatrix.getData();
    }

    // 检查矩阵是否正定
    private boolean checkPositiveDefinite(RealMatrix matrix) {
        try {
            // 尝试 Cholesky 分解
            new CholeskyDecomposition(matrix);
            return true;
        } catch (Exception e) {
            return false;
        }
    }

    // 如果矩阵不是正定的，进行正定处理（计算量太大，不推荐）
    private RealMatrix makePositiveDefinite(RealMatrix matrix) {
        EigenDecomposition eigenDecomposition = new EigenDecomposition(matrix);
        double[] realEigenvalues = eigenDecomposition.getRealEigenvalues();
        RealMatrix eigenvectors = eigenDecomposition.getV();

        // 确保特征值为正
        for (int i = 0; i < realEigenvalues.length; i++) {
            if (realEigenvalues[i] <= 1e-6) {
                realEigenvalues[i] = 1e-4;  // 小正数代替负特征值
            }
        }

        return eigenvectors.multiply(
                new DiagonalMatrix(realEigenvalues)
        ).multiply(eigenvectors.transpose());
    }

    // 检查协方差矩阵是否对称
    private boolean isSymmetric(RealMatrix matrix) {
        return matrix.subtract(matrix.transpose()).getNorm() < 1e-20;    // 小于某个阈值即认为是对称
    }

    // 生成基因型数据
    public void generateAndStandardizeGenotypes(double modifyRate) {

        // 使用协方差矩阵和正态分布来生成基因型数据
        RealMatrix covarianceMatrix = new Array2DRowRealMatrix(this.COV);

        // 对协方差矩阵进行对称化处理
        if (!isSymmetric(covarianceMatrix)){
            logger.warn("协方差矩阵不对称，修正中 ......");
            covarianceMatrix = makeSymmetric(covarianceMatrix);
        } else {
            logger.info("协方差矩阵对称.");
        }

        logger.info("开始生成基因型数据 ......");

        // 使用Cholesky分解得到下三角矩阵 L
        CholeskyDecomposition choleskyDecomposition = new CholeskyDecomposition(covarianceMatrix);
        RealMatrix lowerTriangularMatrix = choleskyDecomposition.getL();

        // 创建一个用于生成标准正态分布的随机数生成器
        RandomDataGenerator randomDataGenerator = new RandomDataGenerator();
        randomDataGenerator.reSeed(rand.nextInt());

        // 保存生成的标准化基因型数据
        for (int i = 0; i < this.N; i++) {  // 对每个样本生成数据
            // 生成标准正态分布的随机数 z_i
            double[] z = new double[this.snpsNum];
            for (int j = 0; j < this.snpsNum; j++) {
                z[j] = randomDataGenerator.nextGaussian(0, 1); // 均值为0，标准差为1
            }

            // 使用Cholesky矩阵L进行线性变换，得到标准化后的基因型数据
            double[] standardizedGenotype = new double[this.snpsNum];
            for (int j = 0; j < this.snpsNum; j++) {
                standardizedGenotype[j] = 0;
                for (int k = 0; k <= j; k++) {  // 只考虑L的下三角
                    standardizedGenotype[j] += lowerTriangularMatrix.getEntry(j, k) * z[k];
                }
            }

            // 基于标准化基因型将其离散化为 0, 1, 2
            int[] genotype = new int[this.snpsNum];
            for (int j = 0; j < this.snpsNum; j++) {
                // 获取当前 SNP 的等位基因频率
                double p = this.P[j];  // SNP j 的等位基因频率
                double q = 1 - p; // B 等位基因的频率

                // 排序标准化基因型
                double[] sortedGenotype = Arrays.copyOf(standardizedGenotype, standardizedGenotype.length);
                Arrays.sort(sortedGenotype);

                // 计算前 p^2, 2pq的阈值位置
                int p2Threshold = (int) (sortedGenotype.length * p * p);  // p^2 的位置
                int pqThreshold = (int) (sortedGenotype.length * (2 * p * q));  // 2pq 的位置

                // 将标准化基因型离散化
                if (standardizedGenotype[j] <= sortedGenotype[p2Threshold]) {
                    genotype[j] = 0;  // 0 表示 AA
                } else if (standardizedGenotype[j] <= sortedGenotype[p2Threshold + pqThreshold]) {
                    genotype[j] = 1;  // 1 表示 AB
                } else {
                    genotype[j] = 2;  // 2 表示 BB
                }
            }

            // 将生成的基因型保存
            this.haplotypes.add(genotype);
            this.standardizeGenotypes.add(standardizedGenotype);
        }

        logger.info("基因型数据已生成！");

        // 这里可以根据需要修改基因型数据，例如引入一些噪音或修正规则
        logger.info("随机交换一些样本，修改基因型 ......");
        modifyGenotypes(modifyRate);
    }

    // 验证生成的基因型数据与真实等位基因频率的匹配
    public void validateGenotypeFrequencies() {
        logger.info("开始验证等位基因频率 ......");

        // 1. 统计模拟生成的基因型数据中的等位基因频率
        double[] observedAlleleCounts = new double[2]; // 统计等位基因 A 和 B 的总计数
        double totalAlleles = this.haplotypes.size() * this.snpsNum * 2; // 总等位基因数，2个等位基因 / SNP

        // 统计每个样本的等位基因（AA, AB, BB）
        for (int[] haplotype : this.haplotypes) {
            for (int genotype : haplotype) {
                // 根据基因型（0, 1, 2）更新等位基因计数
                if (genotype == 0) { // AA 基因型
                    observedAlleleCounts[0] += 2; // 2 个 A 等位基因
                } else if (genotype == 1) { // AB 基因型
                    observedAlleleCounts[0] += 1; // 1 个 A 等位基因
                    observedAlleleCounts[1] += 1; // 1 个 B 等位基因
                } else if (genotype == 2) { // BB 基因型
                    observedAlleleCounts[1] += 2; // 2 个 B 等位基因
                }
            }
        }

        // 2. 计算模拟生成的等位基因频率
        double observedAlleleFrequencyA = observedAlleleCounts[0] / totalAlleles;
        double observedAlleleFrequencyB = observedAlleleCounts[1] / totalAlleles;

        // 3. 计算预期的等位基因频率
        // 等位基因频率 P[i] 已经在 P 数组中存储
        double[] expectedAlleleFrequency = new double[2];
        expectedAlleleFrequency[0] = 0; // A 等位基因频率
        expectedAlleleFrequency[1] = 0; // B 等位基因频率

        for (int i = 0; i < this.snpsNum; i++) {
            double p = this.P[i]; // SNP i 的等位基因频率
            expectedAlleleFrequency[0] += p; // A 等位基因的总频率
            expectedAlleleFrequency[1] += (1 - p); // B 等位基因的总频率
        }

        expectedAlleleFrequency[0] /= this.snpsNum; // 计算所有 SNP 的平均 A 等位基因频率
        expectedAlleleFrequency[1] = 1 - expectedAlleleFrequency[0]; // B 等位基因的频率 = 1 - A 等位基因频率

        // 4. 输出观察到的等位基因频率与期望的频率
        logger.info("观察到的等位基因频率：");
        logger.info("A 等位基因频率: {}", observedAlleleFrequencyA);
        logger.info("B 等位基因频率: {}", observedAlleleFrequencyB);

        logger.info("期望的等位基因频率：");
        logger.info("A 等位基因频率: {}", expectedAlleleFrequency[0]);
        logger.info("B 等位基因频率: {}", expectedAlleleFrequency[1]);

        // 5. 比较观察到的等位基因频率与期望的频率
        double threshold = 0.05; // 允许的差异阈值，例如 5%

        // 比较 A 等位基因频率
        double alleleDifferenceA = Math.abs(observedAlleleFrequencyA - expectedAlleleFrequency[0]);
        if (alleleDifferenceA > threshold) {
            logger.warn("警告：A 等位基因观察频率与期望频率差异较大: {} vs {}", observedAlleleFrequencyA, expectedAlleleFrequency[0]);
        } else {
            logger.info("A 等位基因观察频率与期望频率相近: {} vs {}", observedAlleleFrequencyA, expectedAlleleFrequency[0]);
        }

        // 比较 B 等位基因频率
        double alleleDifferenceB = Math.abs(observedAlleleFrequencyB - expectedAlleleFrequency[1]);
        if (alleleDifferenceB > threshold) {
            logger.warn("警告：B 等位基因观察频率与期望频率差异较大: {} vs {}", observedAlleleFrequencyB, expectedAlleleFrequency[1]);
        } else {
            logger.info("B 等位基因观察频率与期望频率相近: {} vs {}", observedAlleleFrequencyB, expectedAlleleFrequency[1]);
        }

        logger.info("等位基因频率验证完成！");
    }


    // 对协方差矩阵进行对称化处理
    private RealMatrix makeSymmetric(RealMatrix matrix) {
        int rows = matrix.getRowDimension();
        for (int i = 0; i < rows; i++) {
            for (int j = i + 1; j < rows; j++) {
                // 获取对称位置的值并计算平均值
                double avg = (matrix.getEntry(i, j) + matrix.getEntry(j, i)) / 2.0;

                // 直接修正对称位置
                matrix.setEntry(i, j, avg);
                matrix.setEntry(j, i, avg);
            }
        }
        return matrix;
    }

    // 修改基因型数据的逻辑（如果需要）
    private void modifyGenotypes(double modifyRate) {
        // 这里可以加入对生成的基因型数据的修改逻辑
        // 比如错误率，或者模拟不同的突变等

        // 例如对某些基因型数据引入一定的修改（例如随机交换）
        for (int i = 0; i < this.haplotypes.size(); i++) {
            int[] haplotype = this.haplotypes.get(i);
            for (int j = 0; j < haplotype.length; j++) {
                if (rand.nextDouble() < modifyRate) {  // 假设修改率为1%
                    haplotype[j] = rand.nextInt(3); // 随机修改为0, 1, 或2
                }
            }
        }
    }

    public void genotype2ccf(String path) throws IOException {
        logger.info("正在将 基因型数据 导出为ccf文件 ......");
        CCFWriter writer = new CCFWriter(new File(path));
        for(int i=0; i<this.snpsNum; i++){
            writer.addField("snp" + i, IFieldType.get("int8"));
        }
        IRecord record = writer.getRecord();
        int count;
        for(int[] sample: this.haplotypes){
            count = 0;
            for(int snp: sample){
                record.getBox("snp" + count).char2Object(new ImmutableByteCode(String.valueOf(snp)));
                count += 1;
            }
            writer.write(record);
        }
        writer.close();
        logger.info("基因型数据的 ccf文件 已保存至：{}", path);
    }
}
