using System;

namespace Illumina.SecondaryAnalysis
{
    public class Benchmark
    {
        #region members

        private DateTime _startTime;

        #endregion

        // constructor
        public Benchmark()
        {
            Reset();
        }

        /// <summary>
        ///     returns the number of elapsed time since the last reset
        /// </summary>
        public string GetElapsedTime()
        {
            DateTime stopTime = DateTime.Now;
            TimeSpan span = new TimeSpan(stopTime.Ticks - _startTime.Ticks);

            if (span.Days > 0)
            {
                return string.Format("{0}:{1:D2}:{2:D2}:{3:D2}.{4:D1}", span.Days, span.Hours, span.Minutes, span.Seconds, span.Milliseconds / 100);
            }

            return string.Format("{0:D2}:{1:D2}:{2:D2}.{3:D1}", span.Hours, span.Minutes, span.Seconds, span.Milliseconds / 100);
        }

        /// <summary>
        ///     returns the number of elapsed time since the last reset
        /// </summary>
        public string GetElapsedIterationTime(int numUnits, string unitName)
        {
            DateTime stopTime = DateTime.Now;
            TimeSpan span = new TimeSpan(stopTime.Ticks - _startTime.Ticks);

            double unitsPerSecond = numUnits / span.TotalSeconds;

            return string.Format("{0:D2}:{1:D2}:{2:D2}.{3:D1} ({4:0.0} {5}/s)", span.Hours, span.Minutes, span.Seconds, span.Milliseconds / 100, unitsPerSecond, unitName);
        }

        /// <summary>
        ///     resets the benchmark start time
        /// </summary>
        public void Reset()
        {
            _startTime = DateTime.Now;
        }
    }
}