using System;
using System.Collections.Generic;
using System.Linq;
using Canvas.CommandLineParsing.OptionProcessing;
using ILMNcommon.Common;

namespace Canvas.CommandLineParsing.CoreOptionTypes
{
    public class OptionInfo<T> : Option<T>, IOptionInfo
    {
        public string RawDescription { get; }
        public virtual string Description => RawDescription;

        public IEnumerable<string> Names { get; }

        public virtual string GetPrototype()
        {
            return string.Join("|", Names);
        }

        public string Name { get; }

        public OptionInfo(string description, params string[] names)
        {
            if (!names.Any())
                throw new ArgumentException("An option requires at least one name");
            if (names.Any(string.IsNullOrEmpty))
                throw new ArgumentException("option names cannot be null or empty");
            Names = names;
            RawDescription = description;
            Name = Names.MaxBy(s => s.Length);
        }

        protected OptionInfo(IOptionInfo info) : this(info.Description, info.Names.ToArray())
        {
        }

        public override OptionCollection<T> GetOptions()
        {
            return new OptionCollection<T>();
        }

        public override ParsingResult<T> Parse(SuccessfulResultCollection parseInput)
        {
            return ParsingResult<T>.SuccesfulResult(parseInput.Get(this));
        }
    }

    public interface IOptionInfo : IOption
    {
        string RawDescription { get; }
        string Description { get; }
        string Name { get; }
        IEnumerable<string> Names { get; }
        string GetPrototype();
    }
}