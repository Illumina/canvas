using System;
using System.Collections.Generic;
using System.Linq;

namespace Canvas.CommandLineParsing.CoreOptionTypes
{
    public class ValueOptionInfo<T> : OptionInfo<T>
    {
        private readonly bool _required;

        public ValueOptionInfo(bool required, string description, params string[] names) : base(description, names)
        {
            _required = required;
        }

        public ValueOptionInfo(bool required, IOptionInfo info) : this(required, info.RawDescription, info.Names.ToArray())
        {
            _required = required;
        }

        public override string Description
        {
            get
            {
                string description = base.Description;
                if (typeof(T) == typeof(List<string>))
                {
                    if (!description.EndsWith("."))
                        description += ".";
                    description += " Option can be specified multiple times.";
                }
                else if (typeof(T) != typeof(string))
                {
                    throw new ApplicationException(
                        $"Unexpected type for value option info: '{typeof(T)}' . Expected {typeof(string)} or {typeof(List<string>)}");
                }
                if (_required)
                {
                    description += " (required)";
                }
                return description;
            }
        }

        public override string GetPrototype()
        {
            return base.GetPrototype() + "=";
        }
    }
}