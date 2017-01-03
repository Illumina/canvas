using System.Collections.Generic;
using Isas.SequencingFiles;

namespace CanvasBin
{
    /// <summary>
    /// Hit counts by position - counting 0, 1, or several hits per unique 35mer.
    /// </summary>
    public class HitArray
    {
        #region Members
        public byte[] Data = null;
        #endregion

        public HitArray(int length)
        {
            this.Data = new byte[length];
        }
        public HitArray(byte[] data)
        {
            this.Data = data;
        }

        public int CountSetBits()
        {
            int count = 0;
            for (int index = 0; index < Data.Length; index++)
            {
                if (Data[index] > 0) count++;
            }
            return count;
        }

        /// <summary>
        /// Count set bits in the targeted regions
        /// </summary>
        /// <param name="regions">sorted list of manifest regions</param>
        /// <returns></returns>
        public int CountSetBits(List<NexteraManifest.ManifestRegion> regions)
        {
            if (regions == null) { return CountSetBits(); }

            int count = 0;
            int index = -1;
            foreach (var region in regions) 
            {
                if (index < region.Start) // avoid overlapping targeted regions
                {
                    index = region.Start -1; // index is 0-based; manifest coordinates are 1-based.
                }
                for (; index < Data.Length && index < region.End; index++) 
                {
                    if (Data[index] > 0) count++;
                }
            }

            return count;
        }


        public void Set(int i)
        {
            Data[i] = (Data[i] == 255 ? (byte)255 : (byte)(Data[i] + 1));
        }
    }
}
