package edu.sysu.pmglab.ly.simulation.utils;

import edu.sysu.pmglab.bytecode.ImmutableByteCode;
import edu.sysu.pmglab.ccf.CCFMetaItem;
import edu.sysu.pmglab.ccf.CCFWriter;
import edu.sysu.pmglab.ccf.record.IRecord;
import edu.sysu.pmglab.ccf.type.IFieldType;
import edu.sysu.pmglab.ccf.viewer.CCFViewer;
import edu.sysu.pmglab.ccf.viewer.CCFViewerReader;
import edu.sysu.pmglab.ly.simulation.genotype.GenotypeSim;
import org.rosuda.JRI.REXP;
import org.rosuda.JRI.Rengine;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.File;
import java.io.IOException;

public class RDSReader {

    private static final Logger logger = LoggerFactory.getLogger(GenotypeSim.class);

    private final Rengine re;

    public RDSReader() {
        this(Contants.R_HOME, Contants.JAVA_LIBRARY_PATH);
    }

    public RDSReader(String R_HOME, String JAVA_LIBRARY_PATH) {
        // 设置 R 环境路径
        System.setProperty("R_HOME", R_HOME);
        System.setProperty("java.library.path", JAVA_LIBRARY_PATH);

        // 启动 R 引擎
        re = new Rengine(new String[] {"--vanilla"}, false, null);

        // 检查 R 引擎是否启动成功
        if (!re.waitForR()) {
            logger.info("Rengine initialization failed");
        }
    }

    /**
     * 读取 RDS 文件（数据框）并返回该文件的 REXP 对象。
     *
     * @param filePath RDS 文件的路径
     * @return REXP 返回的 R 对象
     */
    public REXP readDataFrameRDS(String filePath) {
        // 在 R 中读取 RDS 文件
        re.eval("data <- readRDS('" + filePath + "')");

        // 获取 R 中的对象
        return re.eval("data");
    }

    /**
     * 将 RDS 文件（这里主要是数据框）转为 CCF 文件。
     *
     * @param rdsFilePath RDS 文件的路径
     * @param ccfFilePath ccf 文件存储路径
     */
    public void map2ccf(String rdsFilePath, String ccfFilePath) throws IOException {
        // [chr, pos, a0, a1, rsid, af_UKBB, ld, pos_hg17, pos_hg18, pos_hg38]
        REXP data = readDataFrameRDS(rdsFilePath);

        // 编码并写入文件
        CCFWriter writer = new CCFWriter(new File(ccfFilePath))
                .addField("VALUE::chr", IFieldType.get("int8"))
                .addField("VALUE::pos", IFieldType.get("varInt32"))
                .addField("VALUE::a0", IFieldType.get("string"))
                .addField("VALUE::a1", IFieldType.get("string"))
                .addField("VALUE::rsid", IFieldType.get("string"))
                .addField("VALUE::af_UKBB", IFieldType.get("float32"))
                .addField("VALUE::ld", IFieldType.get("float32"))
                .addField("VALUE::pos_hg17", IFieldType.get("varInt32"))
                .addField("VALUE::pos_hg18", IFieldType.get("varInt32"))
                .addField("VALUE::pos_hg38", IFieldType.get("varInt32"));

        IRecord record = writer.getRecord();

        int rows = data.asList().at("chr").asIntArray().length;

        for (int i=0; i<rows; i++){
            record.getBox("chr").char2Object(new ImmutableByteCode(String.valueOf(data.asList().at("chr").asIntArray()[i])));
            record.getBox("pos").char2Object(new ImmutableByteCode(String.valueOf(data.asList().at("pos").asIntArray()[i])));
            record.getBox("a0").char2Object(new ImmutableByteCode(data.asList().at("a0").asStringArray()[i]));
            record.getBox("a1").char2Object(new ImmutableByteCode(data.asList().at("a1").asStringArray()[i]));
            record.getBox("rsid").char2Object(new ImmutableByteCode(data.asList().at("rsid").asStringArray()[i]));
            record.getBox("af_UKBB").char2Object(new ImmutableByteCode(String.valueOf(data.asList().at("af_UKBB").asDoubleArray()[i])));
            record.getBox("ld").char2Object(new ImmutableByteCode(String.valueOf(data.asList().at("ld").asDoubleArray()[i])));
            record.getBox("pos_hg17").char2Object(new ImmutableByteCode(String.valueOf(data.asList().at("pos_hg17").asIntArray()[i])));
            record.getBox("pos_hg18").char2Object(new ImmutableByteCode(String.valueOf(data.asList().at("pos_hg18").asIntArray()[i])));
            record.getBox("pos_hg38").char2Object(new ImmutableByteCode(String.valueOf(data.asList().at("pos_hg38").asIntArray()[i])));

            writer.write(record);
        }
        writer.close();
        logger.info("map_LD映射文件转换完毕。");
    }

    /**
     * 读取 RDS 文件（稀疏矩阵）并提取稀疏矩阵的基本元素存储为ccf文件。
     *
     * @param rdsPath  输入的 RDS 文件路径
     * @param ccfPath  输出的 CCF 文件路径
     */
    public void sparseMatrixRDS2ccf(String rdsPath, String ccfPath) throws IOException {
        // 读取 RDS 文件
        re.eval("LD_ref <- readRDS('" + rdsPath + "')");

        // 提取稀疏矩阵的基本元素
        int[] i = re.eval("i <- LD_ref@i").asIntArray();
        int[] p = re.eval("p <- LD_ref@p").asIntArray();
        double[] x = re.eval("x <- LD_ref@x").asDoubleArray();
        int[] Dim = re.eval("Dim <- LD_ref@Dim").asIntArray();

        CCFMetaItem<String[]> meta1 = CCFMetaItem.of("Dim", new String[]{
                String.valueOf(Dim[0]),
                String.valueOf(Dim[1])
        });

        CCFMetaItem<String[]> meta2 = CCFMetaItem.of("scale", new String[]{
                String.valueOf(i.length),
                String.valueOf(p.length),
                String.valueOf(x.length)
        });

        // 存储为ccf文件
        CCFWriter writer = new CCFWriter(new File(ccfPath))
                .addField("VALUE::i", IFieldType.get("varInt32"))
                .addField("VALUE::p", IFieldType.get("varInt32"))
                .addField("VALUE::x", IFieldType.get("float32"))
                .addMeta(meta1)
                .addMeta(meta2);

        IRecord record = writer.getRecord();

        int maxLength = Math.max(Math.max(i.length, p.length), x.length);

        // 写进 CCF 文件
        for (int j = 0; j < maxLength; j++) {
            record.getBox("i").char2Object(new ImmutableByteCode(j < i.length ? String.valueOf(i[j]) : ""));
            record.getBox("p").char2Object(new ImmutableByteCode(j < p.length ? String.valueOf(p[j]) : ""));
            record.getBox("x").char2Object(new ImmutableByteCode(j < x.length ? String.valueOf(x[j]) : ""));
            writer.write(record);
        }

        writer.close();
        logger.info("稀疏矩阵转ccf文件完毕。");
    }

    /**
     * 结束 R 引擎。
     */
    public void close() {
        if (re != null) {
            re.end();
        }
    }
}
