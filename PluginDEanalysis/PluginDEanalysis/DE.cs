using BaseLibS.Graph;
using BaseLibS.Param;
using PerseusApi.Document;
using PerseusApi.Generic;
using PerseusApi.Matrix;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using PluginInterop;
using System.Text;
using PluginDEanalysis.Properties;

namespace PluginDEanalysis
{
    public class DE : PluginInterop.R.MatrixProcessing
    {
        public override string Heading => "DE analysis";
        public override string Name => "Run DE analysis";
        public override string Description => "Run differential analysis";
        public override bool IsActive => true;
        public override string Url => "http://coxdocs.org/doku.php?id=perseus:user:activities:matrixprocessing:DEanalysis:DEcomputation";

        protected override bool TryGetCodeFile(Parameters param, out string codeFile)
        {
            byte[] code = (byte[])Resources.ResourceManager.GetObject("DEcompute");
            codeFile = Path.GetTempFileName();
            File.WriteAllText(codeFile, Encoding.UTF8.GetString(code));
            return true;
        }

        protected override string GetCommandLineArguments(Parameters param)
        {
            var tempFile = Path.GetTempFileName();
            param.ToFile(tempFile);
            return tempFile;
        }

        protected override Parameter[] SpecificParameters(IMatrixData mdata, ref string errString)
        {
            if (mdata.ColumnCount < 3)
            {
                errString = "Please add at least 3 main data columns to the matrix.";
                return null;
            }
            Parameters[] subParams = new Parameters[mdata.CategoryRowCount];
            for (int i = 0; i < mdata.CategoryRowCount; i++)
            {
                string[] values = mdata.GetCategoryRowValuesAt(i);
                int[] sel = values.Length == 1 ? new[] { 0 } : new int[0];
                subParams[i] =
                    new Parameters(new Parameter[]{
                        new SingleChoiceParam("Reference sample"){
                            Values = values,
                            Help = "The group that should be present to compute differential expression analysis.",
                        },new SingleChoiceParam("Target sample"){
                            Values = values,
                            Help = "The group that should be present to compute differential expression analysis."
                        }
                    });
            }
            return new Parameter[]
            {
                new SingleChoiceWithSubParams("Program")
                {
                    Values = new[] { "EdgeR", "DESeq2", "Limma and DEqMS", "SAM", "ROTS" },
                    SubParams = new[]{
                        new Parameters(new Parameter[]{ 
//                            new DoubleParam("Dispersion (without replicates)", 0.04){
//                            Help = "This dispersion value of EdgeR for the dataset without replicates." },
                            new SingleChoiceParam("Test method"){
                                Values = new[] { "Likelihood ratio test", "Quasi-likelihood F-test", "Exact Test" },
                                Help = "This test method for running EdgeR." },
                            new SingleChoiceParam("Normalization method (EdgeR)"){
                                Values = new[] { "TMM", "TMMwsp", "RLE", "upperquartile", "none" },
                                Help = "This method for calculating normalization factor." },
                            new BoolParam("Data in Log2 scale", true) {
                            Help = "Is the data in log2 scale?"}}),
                        new Parameters(new Parameter[] { new SingleChoiceParam("Fit type"){
                            Values = new[] { "parametric", "local", "mean" },
                            Help = "This fit type for running DESeq2. If your dataset without replicates, please use local or mean."},
                            new SingleChoiceWithSubParams("Data scale"){
                            Values = new[] { "Log2 scale", "Raw intensity or read counts" },
                            Help = "The scale of the input data.",
                            SubParams = new[]{new Parameters(new Parameter[] { new DoubleParam("Multiplication for Log2 values", 1000) {
                            Help = "Multiply this value to the main matrix. It will performa after log2."}}),
                            new Parameters(new Parameter[] { new SingleChoiceWithSubParams("Transformation"){
                                    Values = new[] { "Perform Log2", "Divided by constant", "None" },
                                    Help = "Which transformation that you want to apply, If None was selected, " +
                                    "error may occurs because DESeq2 may not be able to handle proteomics data at large scales " +
                                    "or has large values",
                                    SubParams = new[]{ new Parameters(new Parameter[] {new DoubleParam("Multiplication", 1000) {
                                        Help = "Multiply this value the main matrix. It will performa after log2."}}),
                                        new Parameters(new Parameter[] { new DoubleParam("Divided by", 1000) {
                                        Help = "The main matrix will be divided by this value."}}),
                                        new Parameters(new Parameter[] { new BoolParam("Force to run", false) {
                                        Help = "If the matrix contains the values which are too high, DESeq2 may produce error. " +
                                        "Please transform the matrix first. If you really want to try such raw matrix, you can turn this on."}})
                                    }}})
                                    }},
                            new BoolParam("Rounding", true) {
                            Help = "Round up the main matrix. It will perform after data transformation."}
                        }),
                        new Parameters(new Parameter[] { new BoolWithSubParams("Voom") {
                            Help = "Voom is a function in the limma package that modifies RNA-Seq data for use with limma.",
                                Value = false,
                                SubParamsTrue = new Parameters(new DoubleParam("Span", 0.5){
                                Help = "width of the lowess smoothing window as a proportion."})},
                            new BoolParam("Trend", true) {
                            Help = "An intensity-dependent trend is fitted to the prior variances."},
                            new BoolParam("Robust", false){
                            Help = "This is frequently useful to protect the empirical Bayes procedure against outliers." },
                            new SingleChoiceParam("Normalization method"){
                                Values = new[] { "scale", "quantile", "none" },
                                Help = "The normalization methods." },
                            new BoolWithSubParams("DEqMS") {
                            Help = "DEqMS package estimate prior variance for proteins quantified by different number of PSMs.",
                                Value = true,
                                SubParamsTrue = new Parameters(new SingleChoiceParam("Fit method"){
                                    Values = new[] { "loess", "nls", "spline" },
                                    Help = "The method for fitting variance against the number of peptides/PSM counts." },
                                    new SingleChoiceParam("Peptide/PSM counts"){
                                    Values = new[] { "Razor + unique peptides", "Unique peptides" },
                                    Help = "The columns of peptides/PSM counts for each main column."})}
                        }),
                        new Parameters(new Parameter[]{ new SingleChoiceWithSubParams("Method selection"){
                            Values = new[] { "SAM", "SAMseq" },
                            Help = "The type of input omics data.",
                            SubParams = new[]{ new Parameters(new Parameter[] { new SingleChoiceParam("Test statistic"){
                                Values = new[] { "Standard (t - statistic)", "Wilcoxon (Two-sample wilcoxon or Mann-Whitney test)" },
                                Help = "Test statistic to use in two class unpaired case." },
                                new SingleChoiceParam("Time summary type"){
                                Values = new[] { "Slope", "Signdeled area" },
                                Help = "Summary measure for each time course." },
                                new SingleChoiceParam("Regression method"){
                                Values = new[] { "Standard (linear least squares)", "Ranks (linear least squares on ranked data)" },
                                Help = "Regression method for quantitative case." },
                                new IntParam("KNN neighbors", 10){
                                Help = "Number of nearest neighbors to use for imputation of missing features values." }}),
                                new Parameters(new Parameter[] { new IntParam("Number of resamples", 20){
                                Help = "This number is used to construct test statistic." },
                                new IntParam("Number of resamples for permutations", 20){
                                Help = "Number of resamples used to construct test statistic for permutations." },
                                new BoolParam("Rounding up", true) {
                                Help = "Round up the main matrix. It will perform after multiplication."}})}
                                },
                            new SingleChoiceParam("Normalization"){
                                Values = new[] { "scale", "quantile", "none" },
                                Help = "The normalization methods." },
                            new IntParam("Number of permutations", 100){
                            Help = "Number of permutations used to estimate false discovery rates." },
                            new DoubleParam("Delta", 0.4){
                            Help = "Number of permutations used to estimate false discovery rates." },
                            new SingleChoiceParam("Data type"){
                            Values = new[] { "Two class unpaired", "Two class paired" },
                            Help = "The type of input data: Quantitative for a continuous parameter." },
                            new BoolParam("Log2", true){
                            Help = "Whether input data is log2 scaled." },
                            new IntParam("Seed", 100){
                            Help = "Initial seed for random number generator." }}),
                        new Parameters(new Parameter[]{ new IntParam("B", 1000){
                            Help = "The number of bootstrap and permutation resamplings." },
                            new StringParam("K", "Default"){
                                Help = "The largest top list size considered. Default: 1/4 of the features are used." },
                            new SingleChoiceParam("Norm. method"){
                                Values = new[] { "scale", "quantile", "none" },
                                Help = "The normalization methods." },
                            new BoolParam("Paired", false){
                                Help = "Whether a paired test is performed." },
                            new BoolParam("Log2", true){
                                Help = "Whether input data is log2 scaled." },
                            new IntParam("Seed", 100){
                            Help = "Initial seed for random number generator." }})
                    },
                    Help = "The program for doing differential expression analysis.",
                    ParamNameWidth = 150,
                    TotalWidth = 731
                }, new SingleChoiceWithSubParams("Group")
                {
                    Values = mdata.CategoryRowNames,
                    SubParams = subParams,
                    Help = "The categorical row that the analysis should be based on.",
                    ParamNameWidth = 70,
                    TotalWidth = 731
                }, new SingleChoiceWithSubParams("Log2 Fold Change") {
                    Help = "The Log2 Fold Change threshold of the significant features. " + 
                    "Number is for using exact values of Fold change as the cutoff. Percentage " + 
                    "is for using the user-defined percentage of maximum absolute fold change as the cutoff." ,
                    Values = new[] { "Percentage", "Number", "None" },
                    SubParams = new[]{
                        new Parameters(new DoubleParam("Up-regluation", 8) {
                        Help = "The cutoff of logFC for up-regulation. " +
                        "If the expressed values are located at higher percentile of the assigned value, " +
                        "it would be considered as up-regulated genes/proteins."},
                        new DoubleParam("Down-regluation", -8){
                        Help = "The cutoff of logFC for down-regulation. " +
                        "If the expressed values are located at lower percentile of the assigned value, " +
                        "it would be considered as down-regulated genes/proteins."}),
                        new Parameters(new DoubleParam("Up-regluation", 1) {
                        Help = "The minimum value of logFC for up-regulation. " +
                        "If the expressed values are higher than this value, " +
                        "it would be considered as up-regulated genes/proteins."},
                        new DoubleParam("Down-regluation", -1) {
                        Help = "The maximum value of logFC for down-regulation. " +
                        "If the expressed values are lower than this value, " +
                        "it would be considered as down-regulated genes/proteins."}),
                        new Parameters(),
                    },
                    ParamNameWidth = 90,
                    TotalWidth = 731
                }, new BoolWithSubParams("Adjusted p-value (FDR)") {
                    Help = "The Adjusted p-value (FDR) threshold of the significant features.",
                    Value = true,
                    SubParamsTrue = new Parameters(new DoubleParam("Max. Adjusted p-value (FDR)", 0.05)),
                    ParamNameWidth = 90,
                    TotalWidth = 731
                }, new BoolWithSubParams("P-value")
                {
                    Help = "The p-value threshold of the significant features.",
                    Value = false,
                    SubParamsTrue = new Parameters(new DoubleParam("Max. p-value", 0.05)),
                    ParamNameWidth = 90,
                    TotalWidth = 731
                }
            };
        }
    }
}