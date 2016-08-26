using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Illumina.Common;
using Ploeh.AutoFixture;
using Ploeh.AutoFixture.AutoNSubstitute;
using Ploeh.AutoFixture.Kernel;
using Ploeh.AutoFixture.Xunit2;
using Xunit;
using Xunit.Sdk;

namespace UnitTests
{
    /// <summary>
    /// This attribute allows for combining multiple lists of complex input parameters for an xUnit [Theory] while
    /// still supporting everything that the [InlineAutoNSubstituteData] attribute supports.
    /// For details about how lists can be used to provide complex input parameters for xUnit theories, read up on xUnit's [MemberData] attribute 
    /// 
    /// Adapted from http://stackoverflow.com/a/28767887/1999165
    /// </summary>
    public class InlineAutoNSubstituteMemberDataAttribute : DataAttribute
    {
        private readonly object[] _values;
        private readonly Func<IFixture> _createFixture = () => GetFixture();

        public InlineAutoNSubstituteMemberDataAttribute(params object[] values)
        {
            _values = values;
        }

        public static IFixture GetFixture()
        {
            var fixture = new Fixture();
            fixture.Customizations.Add(new PropertyOmitter());
            return fixture.Customize(
                new CompositeCustomization(
                    new AutoNSubstituteCustomization()));
        }

        public Type PropertyHost { get; set; }

        private IEnumerable<object[]> GetAllParameterObjects(MethodInfo methodUnderTest)
        {
            List<IEnumerable<object[]>> members = new List<IEnumerable<object[]>>();
            foreach (var value in _values)
            {
                var objects = GetAllParameterObjects(value as string, methodUnderTest);
                if (objects == null)
                    members.Add(value.ToSingleItemEnumerable().ToArray().ToSingleItemEnumerable());
                else
                    members.Add(objects);
            }

            // we want to produce all possible combinations for the parameters
            return Combinations(members);
        }

        private IEnumerable<object[]> Combinations(IEnumerable<IEnumerable<object[]>> members)
        {
            if (members.Count() == 1) return members.First();
            List<object[]> parameters = new List<object[]>();

            var firstMember = members.First();
            var remainingCombinations = Combinations(members.Skip(1));
            foreach (var remainingCombinationsRow in remainingCombinations)
            {
                foreach (var firstMemberRow in firstMember)
                {
                    parameters.Add(firstMemberRow.Concat(remainingCombinationsRow).ToArray());
                }
            }
            return parameters;
        }

        private IEnumerable<object[]> GetAllParameterObjects(string propertyName, MethodInfo methodUnderTest)
        {
            if (propertyName == null) return null;
            var type = PropertyHost ?? methodUnderTest.DeclaringType;
            var property = type.GetProperty(propertyName, BindingFlags.Static | BindingFlags.Public | BindingFlags.FlattenHierarchy);

            if (property == null)
            {
                return null;
            }
            var obj = property.GetValue(null, null);
            if (obj == null)
                throw new ArgumentException($"Property {propertyName} on {type.FullName} is null");

            var enumerable = obj as IEnumerable<object[]>;
            if (enumerable != null)
                return enumerable;

            var singleEnumerable = obj as IEnumerable<object>;
            if (singleEnumerable != null)
                return singleEnumerable.Select(x => new[] { x });

            throw new ArgumentException($"Property {propertyName} on {type.FullName} did not return IEnumerable<object[]> or IEnumerable<object>");
        }

        private object[] GetObjects(object[] parameterized, ParameterInfo[] parameters, IFixture fixture)
        {
            var result = new object[parameters.Length];

            for (int i = 0; i < parameters.Length; i++)
            {
                if (i < parameterized.Length)
                    result[i] = parameterized[i];
                else
                    result[i] = CustomizeAndCreate(fixture, parameters[i]);
            }

            return result;
        }

        private object CustomizeAndCreate(IFixture fixture, ParameterInfo p)
        {
            var customizations = p.GetCustomAttributes(typeof(CustomizeAttribute), false)
                .OfType<CustomizeAttribute>()
                .Select(attr => attr.GetCustomization(p));

            foreach (var c in customizations)
            {
                fixture.Customize(c);
            }

            var context = new SpecimenContext(fixture);
            return context.Resolve(p);
        }

        public override IEnumerable<object[]> GetData(MethodInfo methodUnderTest)
        {
            foreach (var values in GetAllParameterObjects(methodUnderTest))
            {
                yield return GetObjects(values, methodUnderTest.GetParameters(), _createFixture());
            }
        }
    }

    public class InlineAutoNSubstituteDataAttribute : CompositeDataAttribute
    {
        public InlineAutoNSubstituteDataAttribute(params object[] values)
            : base(
                new InlineDataAttribute(values),
                new AutoDataAttribute(GetFixture()))
        {
        }

        public static IFixture GetFixture()
        {
            var fixture = new Fixture();
            fixture.Customizations.Add(new PropertyOmitter());
            return fixture.Customize(
                new CompositeCustomization(
                    new AutoNSubstituteCustomization()));
        }
    }

    public class AutoNSubstituteDataAttribute : AutoDataAttribute
    {
        public AutoNSubstituteDataAttribute()
            : base(new Fixture().Customize(new AutoNSubstituteCustomization()))
        {
        }
    }

    public class PropertyOmitter : ISpecimenBuilder
    {
        public object Create(object request, ISpecimenContext context)
        {
            var propInfo = request as PropertyInfo;
            if (propInfo != null)
                return new OmitSpecimen();

            return new NoSpecimen(request);
        }
    }
}
