using ProtoBuf;

namespace Isas.Shared
{
    // Note: 2015-10-23: (ImplicitFields = ImplicitFields.AllFields, AsReferenceDefault = true)
    //    causes "ProtoBuf.ProtoException: Possible recursion detected"
    // See TSAW-203

    // Note: 2015-10-23: (ImplicitFields = ImplicitFields.AllFields, AsReferenceDefault = true)
    //    causes "ProtoBuf.ProtoException: Possible recursion detected"
    // See TSAW-203

    /// <summary>
    ///     Store the information in the [Header] section
    /// </summary>
    [ProtoContract]
    public class HeaderSection
    {
        [ProtoMember(1)]
        public string Manifest; // Manifest name from DesignStudio
    }
}
