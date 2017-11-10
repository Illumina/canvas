using CanvasCommon;
using Isas.Framework.DataTypes.Maps;

namespace CanvasPedigreeCaller
{
    internal interface IVariantCaller
    {
        /// <summary>
        /// Identify variant with the highest likelihood at a given setPosition and assign relevant scores
        /// </summary>
        void CallVariant(ISampleMap<CanvasSegment> canvasSegments, ISampleMap<SampleMetrics> samplesInfo,
            ISampleMap<CopyNumberModel> copyNumberModel,
            PedigreeInfo pedigreeInfo);
    }
}