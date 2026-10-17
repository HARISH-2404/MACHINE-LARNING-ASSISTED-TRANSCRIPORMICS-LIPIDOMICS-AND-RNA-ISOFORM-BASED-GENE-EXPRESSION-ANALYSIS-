package org.multiomics;

import java.io.IOException;
import java.nio.file.*;
import java.util.*;

/** Cross-platform Java 17+ launcher for the Python/R multi-omics workflow.
 * Java orchestrates tools; statistical inference remains in DESeq2 and the specialist packages.
 */
public class MultiOmicsPipeline {
    private static final Path ROOT = Paths.get("").toAbsolutePath().normalize();
    private static void run(List<String> command) throws IOException, InterruptedException {
        System.out.println("\n>>> " + String.join(" ", command));
        ProcessBuilder pb = new ProcessBuilder(command).directory(ROOT.toFile()).inheritIO();
        Map<String,String> env = pb.environment();
        Process p = pb.start(); int code = p.waitFor();
        if (code != 0) throw new RuntimeException("Step failed (exit " + code + "): " + command);
    }
    private static String python() { return System.getProperty("python", System.getenv().getOrDefault("PYTHON", "python")); }
    private static String rscript() { return System.getProperty("rscript", System.getenv().getOrDefault("RSCRIPT", "Rscript")); }
    private static void usage() {
        System.out.println("""
        Usage from repository root:
          java -cp java/build org.multiomics.MultiOmicsPipeline demo
          java -cp java/build org.multiomics.MultiOmicsPipeline validate <data_dir>
          java -cp java/build org.multiomics.MultiOmicsPipeline count-qc <counts.csv> <metadata.csv> <outdir>
          java -cp java/build org.multiomics.MultiOmicsPipeline lipidomics <lipid.csv> <metadata.csv> <outdir>
          java -cp java/build org.multiomics.MultiOmicsPipeline integrate <metadata.csv> <gene.csv> <transcript.csv> <lipid.csv> <outdir>
          java -cp java/build org.multiomics.MultiOmicsPipeline deseq2 <counts.csv> <metadata.csv> <outdir>
        Set -Dpython=python3 or -Drscript=Rscript if executable names differ.
        FASTQ QC/Salmon are run as documented in scripts/ and docs; configure reference files first.
        """);
    }
    public static void main(String[] args) throws Exception {
        if (args.length == 0) { usage(); System.exit(2); }
        String py = python();
        switch (args[0]) {
            case "demo" -> {
                run(List.of(py,"scripts/generate_example_data.py"));
                run(List.of(py,"scripts/validate_inputs.py","--data-dir","data/example"));
                run(List.of(py,"python/02_count_matrix_qc.py","--counts","data/example/gene_counts.csv","--metadata","data/example/metadata.csv","--outdir","results/gene_count_qc"));
                run(List.of(py,"scripts/lipidomics_analysis.py","--data","data/example/lipidomics.csv","--metadata","data/example/metadata.csv","--outdir","results/lipidomics"));
                run(List.of(py,"scripts/cross_omics_associations.py","--metadata","data/example/metadata.csv","--gene","data/example/gene_counts.csv","--transcript","data/example/transcript_counts.csv","--lipid","data/example/lipidomics.csv","--outdir","results/cross_omics"));
                run(List.of(py,"scripts/integrate_ml.py","--metadata","data/example/metadata.csv","--gene","data/example/gene_counts.csv","--transcript","data/example/transcript_counts.csv","--lipid","data/example/lipidomics.csv","--outdir","results/integration"));
                System.out.println("Demo complete. Synthetic results are software tests, not biological evidence.");
            }
            case "validate" -> { require(args,2); run(List.of(py,"scripts/validate_inputs.py","--data-dir",args[1])); }
            case "count-qc" -> { require(args,4); run(List.of(py,"python/02_count_matrix_qc.py","--counts",args[1],"--metadata",args[2],"--outdir",args[3])); }
            case "lipidomics" -> { require(args,4); run(List.of(py,"scripts/lipidomics_analysis.py","--data",args[1],"--metadata",args[2],"--outdir",args[3])); }
            case "integrate" -> { require(args,6); run(List.of(py,"scripts/integrate_ml.py","--metadata",args[1],"--gene",args[2],"--transcript",args[3],"--lipid",args[4],"--outdir",args[5])); }
            case "deseq2" -> { require(args,4); run(List.of(rscript(),"R/differential_gene_expression.R",args[1],args[2],args[3])); }
            default -> { usage(); System.exit(2); }
        }
    }
    private static void require(String[] a,int n) { if (a.length != n) { usage(); throw new IllegalArgumentException("Wrong number of arguments"); } }
}
