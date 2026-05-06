import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

public class CrossValidationRunner {
    

    public static void main(String[] args) {
        String foldPrefix = null;
        Integer kFolds = null;
        String trainJar = null;
        String predictJar = null;
        String validateJar = null;
        String outDir = null;
        String mafSource = null;
        String format = "txt";
        String method = null;


        for (int i = 0; i < args.length; i++) {
            switch (args[i]) {
                case "--foldPrefix":
                    if (i + 1 >= args.length) {
                        printHelp();
                        return;
                    }
                    foldPrefix = args[++i];
                    break;

                case "--cv":
                    if (i + 1 >= args.length) {
                        printHelp();
                        return;
                    }
                    kFolds = Integer.parseInt(args[++i]);
                    break;

                case "--trainJar":
                    if (i + 1 >= args.length) {
                        printHelp();
                        return;
                    }
                    trainJar = args[++i];
                    break;

                case "--predictJar":
                    if (i + 1 >= args.length) {
                        printHelp();
                        return;
                    }
                    predictJar = args[++i];
                    break;

                case "--method":
                    if (i + 1 >= args.length) {
                        printHelp();
                        return;
                    }
                    method = args[++i];
                    break;


                case "--validateJar":
                    if (i + 1 >= args.length) {
                        printHelp();
                        return;
                    }
                    validateJar = args[++i];
                    break;

                case "--out":
                    if (i + 1 >= args.length) {
                        printHelp();
                        return;
                    }
                    outDir = args[++i];
                    break;

                case "--mafSource":
                    if (i + 1 >= args.length) {
                        printHelp();
                        return;
                    }
                    mafSource = args[++i];
                    break;

                case "--format":
                    if (i + 1 >= args.length) {
                        printHelp();
                        return;
                    }
                    format = args[++i];
                    break;

                default:
                    printHelp();
                    return;
            }
        }
        if (foldPrefix == null || kFolds == null || trainJar == null || predictJar == null || validateJar == null || outDir == null|| method == null ) {
            printHelp();
            return;
        }
        if (isGor5(method) && mafSource == null) {
            throw new IllegalArgumentException("GOR5 methods require --mafSource <alignment-folder>");
        }



        if (kFolds < 2) {
            throw new IllegalArgumentException("--cv must be at least 2");
        }

        java.io.File out = new java.io.File(outDir);
        if (!out.exists() && !out.mkdirs()) {
            throw new RuntimeException("Could not create output directory: " + outDir);
        }

        List<Double> foldQ3Means = new ArrayList<>();
        List<Double> foldSOVMeans = new ArrayList<>();

        for (int fold = 1; fold <= kFolds; fold++) {

            String trainFile = foldPrefix + "_fold" + fold + "_train.db";
            String testFile = foldPrefix + "_fold" + fold + "_test.db";

            System.out.println("Running fold " + fold);
            String base = baseMethod(method);

            String modelFile = outDir + "/fold" + fold + "_" + base + "_model.txt";
            String testFasta = outDir + "/fold" + fold + "_test.fasta";
            String mafDir = outDir + "/fold" + fold + "_maf";
            String predFile = outDir + "/fold" + fold + "_" + method + "_pred.txt";
            String summaryFile = outDir + "/fold" + fold + "_" + method + "_summary." + format;
            String detailFile = outDir + "/fold" + fold + "_" + method + "_detail." + format;

            // convert test db -> fasta for normal prediction
            if (!isGor5(method)) {
                writeFastaFromDb(testFile, testFasta);
            } else {
                prepareFoldMafDir(testFile, mafSource, mafDir);
            }

            // train model on train fold
            ProcessBuilder trainPb = new ProcessBuilder(
                    "java", "-jar", trainJar,
                    "--db", trainFile,
                    "--method", base,
                    "--model", modelFile
            );
            trainPb.inheritIO();
            printCommand("TRAIN", trainPb);
            runCommand(trainPb);

            // predict on test fasta, redirect stdout to prediction file
            ProcessBuilder predPb;
            if (isGor5(method)) {
                predPb = new ProcessBuilder(
                        "java", "-jar", predictJar,
                        "--model", modelFile,
                        "--format", "txt",
                        "--maf", mafDir
                );
            } else {
                predPb = new ProcessBuilder(
                        "java", "-jar", predictJar,
                        "--model", modelFile,
                        "--format", "txt",
                        "--seq", testFasta
                );
            }
            predPb.redirectOutput(new java.io.File(predFile));
            predPb.redirectError(ProcessBuilder.Redirect.INHERIT);
            printCommand("PREDICT", predPb);
            runCommand(predPb);

            java.io.File predOut = new java.io.File(predFile);
            if (!predOut.exists() || predOut.length() == 0) {
                throw new RuntimeException("Prediction file is empty: " + predFile);
            }

            // validate predictions against the test db
            ProcessBuilder valPb = new ProcessBuilder(
                    "java", "-jar", validateJar,
                    "-p", predFile,
                    "-r", testFile,
                    "-s", summaryFile,
                    "-d", detailFile,
                    "-f", format
            );
            valPb.inheritIO();
            printCommand("VALIDATE", valPb);
            runCommand(valPb);

            double[] foldMetrics = readFoldMeansFromDetail(detailFile);
            foldQ3Means.add(foldMetrics[0]);
            foldSOVMeans.add(foldMetrics[1]);

            System.out.println("done fold " + fold);
            System.out.println();
        }
        String cvSummaryFile = outDir + "/cv_summary.txt";
        try (PrintWriter w = new PrintWriter(new FileWriter(cvSummaryFile))) {
            double q3Mean = mean(foldQ3Means);
            double q3Std = stdDev(foldQ3Means);
            double sovMean = mean(foldSOVMeans);
            double sovStd = stdDev(foldSOVMeans);

            w.println("fold_metric\tmean\tstddev");
            w.println(String.format(Locale.US, "Q3\t%.4f\t%.4f", q3Mean, q3Std));
            w.println(String.format(Locale.US, "SOV\t%.4f\t%.4f", sovMean, sovStd));
        } catch (Exception e) {
            throw new RuntimeException(e);
        }

        System.out.println("wrote CV summary: " + cvSummaryFile);
    }


    //remove this later, only for debuggin
    public static void printCommand(String label, ProcessBuilder pb) {
//        System.out.println(label + ": " + String.join(" ", pb.command()));
    }

    //HELPERS
    public static void runCommand(ProcessBuilder pb) {
        try {
            Process p = pb.start();
            int exit = p.waitFor();
            if (exit != 0) {
                throw new RuntimeException("Command failed with exit code " + exit);
            }
        } catch (Exception e) {
            throw new RuntimeException(e);
        }


    }

    public static void writeFastaFromDb(String dbPath, String fastaPath) {
        List<DatasetEntry> entries = DatasetReader.readSeclib(dbPath);
        try (java.io.PrintWriter w = new java.io.PrintWriter(new java.io.FileWriter(fastaPath))) {
            for (DatasetEntry e : entries) {
                w.println(">" + e.id);
                if (e.as != null) {
                    w.println(e.as);
                }
            }
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    public static void prepareFoldMafDir(String testDbPath, String mafSource, String mafDirPath) {
        java.io.File mafDir = new java.io.File(mafDirPath);
        if (!mafDir.exists() && !mafDir.mkdirs()) {
            throw new RuntimeException("Could not create fold maf directory: " + mafDirPath);
        }

        List<DatasetEntry> testEntries = DatasetReader.readSeclib(testDbPath);
        java.io.File sourceDir = new java.io.File(mafSource);
        java.io.File[] sourceFiles = sourceDir.listFiles();
        if (sourceFiles == null) {
            throw new RuntimeException("Could not list files in maf source: " + mafSource);
        }

        for (DatasetEntry e : testEntries) {
            java.io.File src = findAlignmentFile(sourceFiles, e.id);
            if (src == null) {
                throw new RuntimeException("No alignment file found for test id: " + e.id);
            }

            java.io.File dst = new java.io.File(mafDir, src.getName());
            try {
                java.nio.file.Files.copy(
                        src.toPath(),
                        dst.toPath(),
                        java.nio.file.StandardCopyOption.REPLACE_EXISTING
                );
            } catch (Exception ex) {
                throw new RuntimeException(ex);
            }
        }
    }

    public static java.io.File findAlignmentFile(java.io.File[] files, String id) {
        //  try exact basename match with common extensions
        String[] exts = {".aln", ".ma", ".maf"};
        for (String ext : exts) {
            for (java.io.File f : files) {
                if (!f.isFile()) continue;
                if (f.getName().equals(id + ext)) return f;
            }
        }

        // try prefix match
        java.io.File hit = null;
        for (java.io.File f : files) {
            if (!f.isFile()) continue;
            String n = f.getName();
            if (!(n.endsWith(".aln") || n.endsWith(".ma") || n.endsWith(".maf"))) continue;
            if (n.startsWith(id)) {
                if (hit != null) {
                    throw new RuntimeException("Multiple alignment files match id: " + id);
                }
                hit = f;
            }
        }
        return hit;
    }

    public static double[] readFoldMeansFromDetail(String detailPath) {
        List<Double> q3s = new ArrayList<>();
        List<Double> sovs = new ArrayList<>();

        try (BufferedReader br = new BufferedReader(new FileReader(detailPath))) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.trim();
                if (!line.startsWith(">")) continue;

                String[] parts = line.substring(1).trim().split("\\s+");
                if (parts.length < 3) continue;

                // format: >id Q3 SOV QH QE QC SOV_H SOV_E SOV_C
                String q3Str = parts[1];
                String sovStr = parts[2];

                if (!q3Str.equals("-")) q3s.add(Double.parseDouble(q3Str));
                if (!sovStr.equals("-")) sovs.add(Double.parseDouble(sovStr));
            }
        } catch (Exception e) {
            throw new RuntimeException(e);
        }

        return new double[]{mean(q3s), mean(sovs)};
    }

    public static double mean(List<Double> values) {
        if (values.isEmpty()) return Double.NaN;
        double sum = 0.0;
        for (double v : values) sum += v;
        return sum / values.size();
    }

    public static double stdDev(List<Double> values) {
        if (values.isEmpty()) return Double.NaN;
        double m = mean(values);
        double sumSq = 0.0;
        for (double v : values) {
            double d = v - m;
            sumSq += d * d;
        }
        return Math.sqrt(sumSq / values.size());
    }

    public static void printHelp() {
        System.out.println("Usage: java CrossValidationRunner --foldPrefix <prefix> --cv <k> --trainJar <train.jar> --predictJar <predict.jar> --validateJar <validateGor.jar> --out <out-dir> --method <gor1|gor3|gor4> [--format <txt|html>]");
        System.out.println("--foldPrefix prefix used by TrainSetCreator output files");
        System.out.println("--cv number of folds to run");
        System.out.println("--trainJar path to training jar");
        System.out.println("--predictJar path to prediction jar");
        System.out.println("--validateJar path to validateGor.jar");
        System.out.println("--out output directory for models, predictions, and validation results");
        System.out.println("--mafSource folder containing per-protein alignment files (required for gor5)");
        System.out.println("--method gor method to use (gor1, gor3, gor4, gor5gor4, gor5gor1, gor5gor3)");
        System.out.println("--format output format for validation (default is txt)");
    }

    public static boolean isGor5(String method) {
        return method.startsWith("gor5");
    }

    public static String baseMethod(String method) {
        if (method.equals("gor1") || method.equals("gor3") || method.equals("gor4")) {
            return method;
        }
        if (method.equals("gor5gor1")) return "gor1";
        if (method.equals("gor5gor3")) return "gor3";
        if (method.equals("gor5gor4")) return "gor4";

        throw new IllegalArgumentException("Unknown method: " + method);
    }
}
