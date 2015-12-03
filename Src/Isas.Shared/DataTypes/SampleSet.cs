using Isas.Shared;
using System;
using System.Collections.Generic;
using System.Linq;
using Illumina.SecondaryAnalysis.Workflow;
using Isas.Shared.Utilities;

namespace Illumina.SecondaryAnalysis
{
    /// <summary>
    /// SampleSet is a container to hold data for a collection of uniquely identifiable samples
    /// Extends Dictionary<SampleInfo, T> 
    /// Value T is the data associated with the SampleInfo
    /// If Value T is a subclass of SampleInfo, this is essentially a HashSet<T>
    /// Wrappers/Workflows will typically want to accept a SampleSet as input and/or return a SampleSet as output
    /// Provides Join methods for creating a new SampleSet from multiple SampleSets
    /// </summary>
    public class SampleSet<T> : Dictionary<SampleInfo, T>
    {
        public SampleSet()
        {
        }

        public SampleSet(IDictionary<SampleInfo, T> dict)
            : base(dict)
        {
        }

        public SampleSet(IEnumerable<KeyValuePair<SampleInfo, T>> keyValuePairs)
            : base(keyValuePairs.ToDictionary(kvp => kvp.Key, kvp => kvp.Value))
        { }

        public IEnumerable<SampleInfo> SampleInfos
        {
            get { return Keys; }
        }

        public IEnumerable<T> SampleData
        {
            get { return Values; }
        }

        /// <summary>
        /// Create a new SampleSet by combining the contents of several individual SampleSets. 
        /// All SampleSets being joined should contain the same set of samples.
        /// E.g. SampleSet&lt;Vcf&gt; + SampleSet&lt;Bam&gt; =&gt; SampleSet&lt;Bam&amp;Vcf&gt;
        /// </summary>
        /// <typeparam name="X">Type of content of SampleSet to be joined</typeparam>
        /// <typeparam name="Y">Type of content of the resulting joined SampleSet</typeparam>
        /// <param name="set1">SampleSet to be joined with this SampleSet</param>
        /// <param name="c">Function creating a content type Y from a content type T and X</param>
        /// <returns>Combined SampleSet&lt;Y&gt; containing all samples in this SampleSet</returns>
        public SampleSet<Y> Join<X, Y>(SampleSet<X> set1, Func<T, X, Y> c)
        {
            return new SampleSet<Y>(this.Keys.ToDictionary(sample => sample,
                sample => c(this[sample], set1[sample])));
        }

        /// <summary>
        /// Create a new SampleSet by combining the contents of several individual SampleSets. 
        /// All SampleSets being joined should contain the same set of samples.
        /// </summary>
        public SampleSet<Y> Join<X, Y>(SampleSet<X> set1, Func<SampleInfo, T, X, Y> c)
        {
            return new SampleSet<Y>(this.Keys.ToDictionary(sample => sample,
                sample => c(sample, this[sample], set1[sample])));
        }

        /// <summary>
        /// Create a new SampleSet by combining the contents of several individual SampleSets. 
        /// All SampleSets being joined should contain the same set of samples.
        /// </summary>
        public SampleSet<Y> Join<X1, X2, Y>(SampleSet<X1> set1, SampleSet<X2> set2, Func<T, X1, X2, Y> c)
        {
            return new SampleSet<Y>(this.Keys.ToDictionary(sample => sample,
                sample => c(this[sample], set1[sample], set2[sample])));
        }

        /// <summary>
        /// Create a new SampleSet by combining the contents of several individual SampleSets. 
        /// All SampleSets being joined should contain the same set of samples.
        /// </summary>
        public SampleSet<Y> Join<X1, X2, Y>(SampleSet<X1> set1, SampleSet<X2> set2, Func<SampleInfo, T, X1, X2, Y> c)
        {
            return new SampleSet<Y>(this.Keys.ToDictionary(sample => sample,
                sample => c(sample, this[sample], set1[sample], set2[sample])));
        }

        /// <summary>
        /// Create a new SampleSet by combining the contents of several individual SampleSets. 
        /// All SampleSets being joined should contain the same set of samples.
        /// </summary>
        public SampleSet<Y> Join<X1, X2, X3, Y>(SampleSet<X1> set1, SampleSet<X2> set2, SampleSet<X3> set3, Func<T, X1, X2, X3, Y> c)
        {
            return new SampleSet<Y>(this.Keys.ToDictionary(sample => sample,
                sample => c(this[sample], set1[sample], set2[sample], set3[sample])));
        }

        /// <summary>
        /// Create a new SampleSet by combining the contents of several individual SampleSets. 
        /// All SampleSets being joined should contain the same set of samples.
        /// </summary>
        public SampleSet<Y> Join<X1, X2, X3, Y>(SampleSet<X1> set1, SampleSet<X2> set2, SampleSet<X3> set3, Func<SampleInfo, T, X1, X2, X3, Y> c)
        {
            return new SampleSet<Y>(this.Keys.ToDictionary(sample => sample,
                sample => c(sample, this[sample], set1[sample], set2[sample], set3[sample])));
        }

        /// <summary>
        /// Create a new SampleSet by combining the contents of several individual SampleSets. 
        /// All SampleSets being joined should contain the same set of samples.
        /// </summary>
        public SampleSet<Y> Join<X1, X2, X3, X4, Y>(SampleSet<X1> set1, SampleSet<X2> set2, SampleSet<X3> set3, SampleSet<X4> set4, Func<T, X1, X2, X3, X4, Y> c)
        {
            return new SampleSet<Y>(this.Keys.ToDictionary(sample => sample,
                sample => c(this[sample], set1[sample], set2[sample], set3[sample], set4[sample])));
        }

        /// <summary>
        /// Create a new SampleSet by combining the contents of several individual SampleSets. 
        /// All SampleSets being joined should contain the same set of samples.
        /// </summary>
        public SampleSet<Y> Join<X1, X2, X3, X4, Y>(SampleSet<X1> set1, SampleSet<X2> set2, SampleSet<X3> set3, SampleSet<X4> set4, Func<SampleInfo, T, X1, X2, X3, X4, Y> c)
        {
            return new SampleSet<Y>(this.Keys.ToDictionary(sample => sample,
                sample => c(sample, this[sample], set1[sample], set2[sample], set3[sample], set4[sample])));
        }

        /// <summary>
        /// Create a new SampleSet containing the samples from both this SampleSet and another SampleSet
        /// The SampleSets being combined must not have any samples in common.
        /// </summary>
        public SampleSet<T> Union(SampleSet<T> set)
        {
            SampleSet<T> result = new SampleSet<T>(this);
            foreach (var kvp in set)
            {
                if (ContainsKey(kvp.Key))
                    throw new ArgumentException("SampleSet already contains this SampleInfo");
                result.Add(kvp.Key, kvp.Value);
            }
            return result;
        }

        public SampleSet<T> WhereData(Func<T, bool> predicate)
        {
            SampleSet<T> result = new SampleSet<T>();
            foreach (var kvp in this)
            {
                if (predicate(kvp.Value))
                    result.Add(kvp.Key, kvp.Value);
            }
            return result;
        }

        public SampleSet<T> WhereData(Func<SampleInfo, T, bool> predicate)
        {
            SampleSet<T> result = new SampleSet<T>();
            foreach (var kvp in this)
            {
                if (predicate(kvp.Key, kvp.Value))
                    result.Add(kvp.Key, kvp.Value);
            }
            return result;
        }

        public SampleSet<Y> SelectData<Y>(Func<T, Y> select)
        {
            SampleSet<Y> result = new SampleSet<Y>();
            foreach (var kvp in this)
            {
                result.Add(kvp.Key, select(kvp.Value));
            }
            return result;
        }

        public SampleSet<Y> SelectData<Y>(Func<SampleInfo, T, Y> select)
        {
            SampleSet<Y> result = new SampleSet<Y>();
            foreach (var kvp in this)
            {
                result.Add(kvp.Key, select(kvp.Key, kvp.Value));
            }
            return result;
        }

        /// <summary>
        /// Return a new SampleSet containing only entries for the specified SampleInfos
        /// </summary>
        public SampleSet<T> Subset(IEnumerable<SampleInfo> infos)
        {
            return new SampleSet<T>(this.SelectKeys(infos));
        }

        /// <summary>
        /// Return a new SampleSet excluding entries for the specified SampleInfos
        /// </summary>
        public SampleSet<T> Exclude(IEnumerable<SampleInfo> infos)
        {
            return new SampleSet<T>(this.ExcludeKeys(infos));
        }
    }

    public static class SampleSetExtensions
    {
        /// <summary>
        /// Create a new SampleSet by combining the contents of several individual SampleSets. 
        /// All SampleSets being joined should contain the same set of samples.
        /// </summary>
        public static SampleSet<Y> Join<X1, X2, Y>(SampleSet<X1> set1, SampleSet<X2> set2, Func<X1, X2, Y> function)
        {
            return set1.Join(set2, function);
        }

        /// <summary>
        /// Create a new SampleSet by combining the contents of several individual SampleSets. 
        /// All SampleSets being joined should contain the same set of samples.
        /// </summary>
        public static SampleSet<Y> Join<X1, X2, Y>(SampleSet<X1> set1, SampleSet<X2> set2, Func<SampleInfo, X1, X2, Y> function)
        {
            return set1.Join(set2, function);
        }

        /// <summary>
        /// Create a new SampleSet by combining the contents of several individual SampleSets. 
        /// All SampleSets being joined should contain the same set of samples.
        /// </summary>
        public static SampleSet<Y> Join<X1, X2, X3, Y>(SampleSet<X1> set1, SampleSet<X2> set2, SampleSet<X3> set3, Func<X1, X2, X3, Y> function)
        {
            return set1.Join(set2, set3, function);
        }

        /// <summary>
        /// Create a new SampleSet by combining the contents of several individual SampleSets. 
        /// All SampleSets being joined must contain the same set of samples.
        /// </summary>
        public static SampleSet<Y> Join<X1, X2, X3, Y>(SampleSet<X1> set1, SampleSet<X2> set2, SampleSet<X3> set3, Func<SampleInfo, X1, X2, X3, Y> function)
        {
            return set1.Join(set2, set3, function);
        }

        /// <summary>
        /// Create a new SampleSet by combining the contents of several individual SampleSets. 
        /// All SampleSets being joined should contain the same set of samples.
        /// </summary>
        public static SampleSet<Y> Join<X1, X2, X3, X4, Y>(SampleSet<X1> set1, SampleSet<X2> set2, SampleSet<X3> set3, SampleSet<X4> set4, Func<X1, X2, X3, X4, Y> function)
        {
            return set1.Join(set2, set3, set4, function);
        }

        /// <summary>
        /// Create a new SampleSet by combining the contents of several individual SampleSets. 
        /// All SampleSets being joined should contain the same set of samples.
        /// </summary>
        public static SampleSet<Y> Join<X1, X2, X3, X4, Y>(SampleSet<X1> set1, SampleSet<X2> set2, SampleSet<X3> set3, SampleSet<X4> set4, Func<SampleInfo, X1, X2, X3, X4, Y> function)
        {
            return set1.Join(set2, set3, set4, function);
        }

        /// <summary>
        /// Create a new SampleSet by combining the contents of several individual SampleSets. 
        /// All SampleSets being joined should contain the same set of samples.
        /// </summary>
        public static SampleSet<Y> Join<X1, X2, X3, X4, X5, Y>(SampleSet<X1> set1, SampleSet<X2> set2, SampleSet<X3> set3, SampleSet<X4> set4, SampleSet<X5> set5, Func<X1, X2, X3, X4, X5, Y> function)
        {
            return set1.Join(set2, set3, set4, set5, function);
        }

        /// <summary>
        /// Create a new SampleSet by combining the contents of several individual SampleSets. 
        /// All SampleSets being joined should contain the same set of samples.
        /// </summary>
        public static SampleSet<Y> Join<X1, X2, X3, X4, X5, Y>(SampleSet<X1> set1, SampleSet<X2> set2, SampleSet<X3> set3, SampleSet<X4> set4, SampleSet<X5> set5, Func<SampleInfo, X1, X2, X3, X4, X5, Y> function)
        {
            return set1.Join(set2, set3, set4, set5, function);
        }
    }

    public class SampleSetLoadingConvention<TIn, T> : ILoadingConvention<SampleSet<TIn>, SampleSet<T>>
    {
        private readonly NamingConventionGenerator _loadingConventionGenerator;
        private readonly SampleStubNamingConvention _convention;

        public delegate ILoadingConvention<TIn, T> NamingConventionGenerator(IFileLocation fileNameStub);

        public SampleSetLoadingConvention(SampleStubNamingConvention sampleStubNamingConvention, NamingConventionGenerator loadingConventionGenerator)
        {
            _loadingConventionGenerator = loadingConventionGenerator;
            _convention = sampleStubNamingConvention;
        }

        public SampleSet<T> Load(SampleSet<TIn> input)
        {
            var set = new SampleSet<T>();
            foreach (KeyValuePair<SampleInfo, TIn> sample in input)
            {
                IFileLocation stub = _convention(sample.Key);
                ILoadingConvention<TIn, T> namingConvention = _loadingConventionGenerator(stub);
                T output = namingConvention.Load(sample.Value);
                set.Add(sample.Key, output);
            }
            return set;
        }

        public void Move(SampleSet<T> source, Action<IFileLocation, IFileLocation> move)
        {
            foreach (var sampleInfo in source.Keys)
            {
                IFileLocation stub = _convention(sampleInfo);
                ILoadingConvention<TIn, T> namingConvention = _loadingConventionGenerator(stub);
                namingConvention.Move(source[sampleInfo], move);
            }
        }
    }
}
