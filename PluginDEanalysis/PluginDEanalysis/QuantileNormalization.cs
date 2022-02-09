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
    public class QuantileNorm : PluginInterop.R.MatrixProcessing
    {
        public override string Heading => "Normalization";
        public override string Name => "Quantile normalization";
        public override string Description => "Run quantile normalization";
        public override bool IsActive => true;
        public override string Url => "http://coxdocs.org/doku.php?id=perseus:user:activities:matrixprocessing:normalization:QuantileNormalization";

        protected override bool TryGetCodeFile(Parameters param, out string codeFile)
        {
            byte[] code = (byte[])Resources.ResourceManager.GetObject("QuantileNorm");
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
            return new Parameter[] { };
        }
    }
}