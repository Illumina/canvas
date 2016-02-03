using Isas.Shared;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Isas.Shared
{
    public interface IMoveableResult<T>
    {
        /// <summary>
        /// Move the relevant members (those that identify file system resources) of the object
        /// to match those in the <paramref name="destination"/> object.
        /// </summary>
        /// <param name="destination"></param>
        /// <returns>A new object where relevant fields have been updated to match the <paramref name="destination"/> object</returns>
        T Move(T destination);
    }

    public static class MoveableResultSampleSetExtensions 
    {
        public static SampleSet<T> Move<T>(this SampleSet<T> samples, SampleSet<T> destinations) where T : IMoveableResult<T>
        {
            SampleSet<T> results = new SampleSet<T>();
            foreach (SampleInfo sample in destinations.Keys)
            {
                results.Add(sample, samples[sample].Move(destinations[sample]));
            }
            return results;
        }

        public static SampleSet<IEnumerable<T>> Move<T>(this SampleSet<IEnumerable<T>> samples, SampleSet<IEnumerable<T>> destinations) where T : IMoveableResult<T>
        {
            SampleSet<IEnumerable<T>> results = new SampleSet<IEnumerable<T>>();
            foreach (SampleInfo sample in destinations.Keys)
            {
                results.Add(sample, samples[sample].Move(destinations[sample]));
            }
            return results;
        }

        public static IEnumerable<T> Move<T>(this IEnumerable<T> samples, IEnumerable<T> destinations) where T : IMoveableResult<T>
        {
            List<T> sampleList = samples.ToList();
            List<T> destList = destinations.ToList();
            List<T> results = new List<T>();
            if (sampleList.Count != destList.Count)
            {
                throw new Exception("Source and destination lists of IMoveableResult's must be the same length.");
            }
            for (int i = 0; i < sampleList.Count; i++)
            {
                results.Add(sampleList[i].Move(destList[i]));
            }
            return results;
        }

        public static SampleSet<IFileLocation> MoveAndLink(this SampleSet<IFileLocation> samples, SampleSet<IFileLocation> destinations)
        {
            return samples.SelectData((sample, file) => file.MoveAndLink(destinations[sample]));
        }
    }
}
