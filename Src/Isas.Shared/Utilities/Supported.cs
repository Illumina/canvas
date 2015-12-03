using System;

namespace Isas.Shared
{

    public static class SupportedNames
    {
        static public SupportedCNVCallers ParseCNVCallerName(string name)
        {
            string tempSetting = name.ToLowerInvariant();
            switch (tempSetting)
            {
                case "0":
                case "false":
                case "none":
                    return SupportedCNVCallers.None;
                case "1":
                case "true":
                case "canvas":
                    return SupportedCNVCallers.Canvas;
                default:
                    throw new Exception($"Invalid RunCNVDetection setting: {name}");
            }
        }
        static public SupportedSVCallers ParseSVCallerName(string name)
        {
            string tempSetting = name.ToLowerInvariant();
            switch (tempSetting)
            {
                case "0":
                case "false":
                case "none":
                    return SupportedSVCallers.None;
                case "true":
                case "1":
                case "manta":
                    return SupportedSVCallers.Manta;
                default:
                    throw new Exception($"Invalid SV caller name: {name}");
            }
        }
    }

	// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
	public enum SupportedAligners
    {
        bwa,
		bwamem,
        AmpliconSmithWaterman,
        Isaac,
        bowtie
		// ReSharper restore InconsistentNaming
	};

	// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
	public enum SupportedVariantCallers
    {
        None = 0,
        GATK,
        Somatic,
        Starling,
        Samtools,
        Strelka
		// ReSharper restore InconsistentNaming
	};

    public enum SupportedSVCallers
    {
        None = 0,
        Manta
    }

    public enum SupportedCNVCallers
    {
        None = 0,
        Canvas
    }

    public enum SupportedVariantAnnotators
    {
        None = 0,
        SVAnnotate,
        Nirvana
    }

	// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
	public enum IndelRealignment
    {
        None = 0,
        GATK,
        Isaac,
        Both,
        Isas,  //Realign only
        IsasThree, //Realign and shift 3'
        IsasFive, //Realign and shift 5'
        Three, //Shift 3' only
        Five //Shift 5' only
		// ReSharper restore InconsistentNaming
	}
}