using System.Collections.Generic;
using System.Linq;
using ILMNcommon.Common;

namespace SequencingFiles.Vcf
{
    public interface IOverlappable
    {
        bool VariantBasesOverlap(ReferenceInterval interval);
    }

    public class VariantIntervalOverlapper : IOverlappable
    {
        private readonly ReferenceInterval _referenceInterval;

        public VariantIntervalOverlapper(ReferenceInterval referenceInterval)
        {
            _referenceInterval = referenceInterval;
        }

        public bool VariantBasesOverlap(ReferenceInterval interval)
        {
            return interval.Overlaps(_referenceInterval);
        }
    }

    public class NonVariantIntervalOverlapper : IOverlappable
    {
        public bool VariantBasesOverlap(ReferenceInterval interval)
        {
            return false;
        }
    }

    public class InsertionOverlapper : IOverlappable
    {
        private readonly ReferencePosition _positionPrecedingInsertion;
        private readonly string _followingReferenceSequence;
        private readonly string _insertedSequence;

        public InsertionOverlapper(ReferencePosition positionPrecedingInsertion, string followingReferenceSequence, string insertedSequence)
        {
            _positionPrecedingInsertion = positionPrecedingInsertion;
            _followingReferenceSequence = followingReferenceSequence;
            _insertedSequence = insertedSequence;
        }

        public bool VariantBasesOverlap(ReferenceInterval interval)
        {
            return PossibleInsertionPositions.Any(position => VariantBasesOverlap(position, interval));
        }

        private IEnumerable<ReferencePosition> PossibleInsertionPositions
        {
            get
            {
                yield return _positionPrecedingInsertion;

                var shortestRepeat = _insertedSequence.ShortestRepeatingSubstring();
                var followingReferenceSequence = _followingReferenceSequence;
                int insertIndex = 0;
                while (followingReferenceSequence.StartsWith(shortestRepeat))
                {
                    insertIndex += shortestRepeat.Length;
                    yield return _positionPrecedingInsertion.Shift(insertIndex);
                    followingReferenceSequence = followingReferenceSequence.Substring(shortestRepeat.Length);
                }
            }
        }

        private bool VariantBasesOverlap(ReferencePosition positionPrecedingInsertion, ReferenceInterval interval)
        {
            var positionAfterInsertion = positionPrecedingInsertion.Next();
            return
                positionPrecedingInsertion.Overlaps(interval) &&
                positionAfterInsertion.Overlaps(interval);
        }
    }
}