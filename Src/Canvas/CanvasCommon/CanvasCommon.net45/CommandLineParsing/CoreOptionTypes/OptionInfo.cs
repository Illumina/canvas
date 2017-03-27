using System;
using System.Collections.Generic;
using System.Linq;
using CanvasCommon.CommandLineParsing.OptionProcessing;
using Illumina.Common;

namespace CanvasCommon.CommandLineParsing.CoreOptionTypes
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
            return ParsingResult<T>.SuccessfulResult(parseInput.Get(this));
        }

        public override bool Equals(object obj)
        {
            // If parameter is null return false.
            if (obj == null)
            {
                return false;
            }

            // If parameter cannot be cast to OptionInfo<T> return false.
            var p = obj as OptionInfo<T>;
            if (p == null)
            {
                return false;
            }

            // Return true if the fields match:
            return Equals(p);
        }

        public bool Equals(OptionInfo<T> p)
        {
            // If parameter is null return false:
            if (p == null)
            {
                return false;
            }

            return Name == p.Name;
        }

        public override int GetHashCode()
        {
            return Name.GetHashCode();
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