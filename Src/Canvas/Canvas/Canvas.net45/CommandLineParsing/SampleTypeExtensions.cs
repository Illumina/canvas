using System;

namespace Canvas.CommandLineParsing
{
    public static class SampleTypeExtensions
    {
        public static string GetOptionName(this SampleType type)
        {
            switch (type)
            {
                case SampleType.Proband:
                    return SmallPedigreeOptionsParser.ProbandOptionName;
                case SampleType.Mother:
                    return SmallPedigreeOptionsParser.MotherOptionName;
                case SampleType.Father:
                    return SmallPedigreeOptionsParser.FatherOptionName;
                default:
                    throw new ArgumentOutOfRangeException(nameof(type), type, null);
            }
        }
    }
}