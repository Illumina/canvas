using System;
using System.Collections.Generic;
using System.Diagnostics;
using SequencingFiles;

namespace Isas.Shared
{
    public enum GenotypeMethod
    {
        MaxGT,
        Poly
    }

	// ReSharper disable InconsistentNaming - prevents ReSharper from renaming serializeable members that are sensitive to being changed
	public class dbSNPEntry
    {
        public float AverageHeterozygosity;
        public int Bin;
        public string Chromosome;
        public int ChromosomeEnd;
        public int ChromosomeStart;
        public string Class;
        public string Function;
        public string LociType;
        public string MoleculeType;
        public string Observed;
        public string ReferenceSNPIdentifier;
        public string ReferenceUCSC;
        public string ReferencedbSNP;
        public int Score;
        public float StandardErrorHeterozygosity;
        public char Strand;
        public string Valid;
        public int Weight;
		// ReSharper restore InconsistentNaming
	}

    public class RefGeneEntry
    {
        public int Bin;
        public string Chromosome;
        public int CodingEnd;
        public string CodingEndStat;
        public int CodingStart;
        public string CodingStartStat;
        public int[] ExonEnd;
        public int[] ExonFrame;
        public int[] ExonStart;
        public string GeneID;
        public int NumberOfExons;
        public int Score;
        public char Strand;
        public string TranscriptID;
        public int TranscriptionEnd;
        public int TranscriptionStart;
    }
}