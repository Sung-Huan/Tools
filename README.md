# PluginDEanalysis

The integration of the most widely used differential expression analysis packages

- DESeq2 (doi: 10.1186/s13059-014-0550-8)
- EdgeR (doi: 10.1093/bioinformatics/btp616)
- Limma (doi: 10.1093/nar/gkv007)
- DEqMS (doi: 10.1074/mcp.TIR119.001646)
- SAM (doi: 10.1073/pnas.091062498)
- ROTS (doi: 10.1371/journal.pcbi.1005562)

For the usage

1. Download and install Visual Studio for C# compile (https://visualstudio.microsoft.com/)
2. Install PerseusR (https://github.com/cox-labs/PerseusR). The protocols of generating Perseus plugin can be viewed at https://doi.org/10.1002/cpbi.105
3. Compile and build the scripts in Visual Studio
4. Copy PluginDEanalysis.dll and PluginDEanalysis.pdb from PluginDEanalysis\bin\Debug\netstandard2.0 to bin folder of your Perseus (ex: Perseus\bin)
