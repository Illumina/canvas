using CanvasCommon.CommandLineParsing.CoreOptionTypes;

namespace Canvas.CommandLineParsing
{
    public abstract class ModeParser : Option<IModeRunner>
    {
        public string Name { get; }
        public string Description { get; }

        protected ModeParser(string name, string description)
        {
            Name = name;
            Description = description;
        }
    }
}