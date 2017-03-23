using System.Collections.Generic;

namespace CanvasCommon.CommandLineParsing.CoreOptionTypes
{
    public class MultipleValueOptionInfo : ValueOptionInfo<List<string>>
    {
        public MultipleValueOptionInfo(bool required, string description, params string[] names) : base(required, description, names)
        {
        }

        public MultipleValueOptionInfo(bool required, IOptionInfo optionInfo) : base(required, optionInfo)
        {
        }
    }
}