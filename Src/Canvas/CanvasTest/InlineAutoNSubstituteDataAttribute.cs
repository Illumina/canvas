using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using Ploeh.AutoFixture;
using Ploeh.AutoFixture.AutoNSubstitute;
using Ploeh.AutoFixture.Kernel;
using Ploeh.AutoFixture.Xunit2;
using Xunit;
using Xunit.Sdk;

namespace UnitTests
{
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
