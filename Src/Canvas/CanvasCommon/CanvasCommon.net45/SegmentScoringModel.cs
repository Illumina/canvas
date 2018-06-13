using System;
using System.Collections.Generic;
using System.Linq;

namespace CanvasCommon
{
    public partial class CanvasSegment
    {
        /// <summary>
        /// Computes QScore using one of the available methods
        /// </summary>
        public enum QScoreMethod
        {
            BinCountLinearFit,
            GeneralizedLinearFit,
            Logistic,
            LogisticGermline
        };

        public int ComputeQScore(QScoreMethod qscoreMethod, QualityScoreParameters qscoreParameters)
        {
            double score;
            int qscore;
            switch (qscoreMethod)
            {
                case QScoreMethod.LogisticGermline:
                    // Logistic model using a new selection of features.  Gives ROC curve area 0.921
                    score = qscoreParameters.LogisticGermlineIntercept;
                    score += GetQScorePredictor(QScorePredictor.LogBinCount) *
                             qscoreParameters.LogisticGermlineLogBinCount;
                    score += GetQScorePredictor(QScorePredictor.ModelDistance) *
                             qscoreParameters.LogisticGermlineModelDistance;
                    score += GetQScorePredictor(QScorePredictor.DistanceRatio) *
                             qscoreParameters.LogisticGermlineDistanceRatio;
                    score = Math.Exp(score);
                    score = score / (score + 1);
                    // Transform probability into a q-score:
                    qscore = (int) (Math.Round(-10 * Math.Log10(1 - score)));
                    qscore = Math.Min(40, qscore);
                    qscore = Math.Max(2, qscore);
                    return qscore;
                case QScoreMethod.Logistic:
                    // Logistic model using a new selection of features.  Gives ROC curve area 0.8289
                    score = qscoreParameters.LogisticIntercept;
                    score += GetQScorePredictor(QScorePredictor.LogBinCount) * qscoreParameters.LogisticLogBinCount;
                    score += GetQScorePredictor(QScorePredictor.ModelDistance) * qscoreParameters.LogisticModelDistance;
                    score += GetQScorePredictor(QScorePredictor.DistanceRatio) * qscoreParameters.LogisticDistanceRatio;
                    score += GetQScorePredictor(QScorePredictor.BinCountAmpDistance);
                    double coreScore = score;
                    score = Math.Exp(score);
                    score = score / (score + 1);
                    // Transform probability into a q-score:
                    qscore = (int) Math.Round(-10 * Math.Log10(1 - score));
                    qscore = Math.Min(60, qscore);
                    qscore = Math.Max(2, qscore);
                    //if (CopyNumber>20)
                    //{
                    //    Console.WriteLine($"HiCN: {CopyNumber} from {this.MedianCount} distance {ModelDistance} next {RunnerUpModelDistance} bins {BinCount}");
                    //    Console.WriteLine($"      Prelogit {coreScore} = intercept {qscoreParameters.LogisticIntercept} + bins { GetQScorePredictor(QScorePredictor.LogBinCount) * qscoreParameters.LogisticLogBinCount} + dist {GetQScorePredictor(QScorePredictor.ModelDistance) * qscoreParameters.LogisticModelDistance} + ratio { GetQScorePredictor(QScorePredictor.DistanceRatio) * qscoreParameters.LogisticDistanceRatio} + amp {GetQScorePredictor(QScorePredictor.BinCountAmpDistance)}");
                    //    Console.WriteLine($"      Logit {score} --> init qscore {(int)Math.Round(-10 * Math.Log10(1 - score))} --> qscore");
                    //}
                    return qscore;
                case QScoreMethod.BinCountLinearFit:
                    if (this.BinCount >= 100)
                        return 61;
                    else
                        return
                            (int)
                            Math.Round(-10 * Math.Log10(1 - 1 / (1 + Math.Exp(0.5532 - this.BinCount * 0.147))), 0,
                                MidpointRounding.AwayFromZero);
                case QScoreMethod.GeneralizedLinearFit: // Generalized linear fit with linear transformation to QScore
                    double linearFit = qscoreParameters.GeneralizedLinearFitIntercept;
                    linearFit += qscoreParameters.GeneralizedLinearFitLogBinCount *
                                 GetQScorePredictor(QScorePredictor.LogBinCount);
                    linearFit += qscoreParameters.GeneralizedLinearFitModelDistance *
                                 GetQScorePredictor(QScorePredictor.ModelDistance);
                    linearFit += qscoreParameters.GeneralizedLinearFitMajorChromosomeCount *
                                 GetQScorePredictor(QScorePredictor.MajorChromosomeCount);
                    linearFit += qscoreParameters.GeneralizedLinearFitMafMean *
                                 GetQScorePredictor(QScorePredictor.MafMean);
                    linearFit += qscoreParameters.GeneralizedLinearFitLogMafCv *
                                 GetQScorePredictor(QScorePredictor.LogMafCv);
                    linearFit += GetQScorePredictor(QScorePredictor.BinCountAmpDistance);
                    score = -11.9 - 11.4 * linearFit; // Scaling to achieve 2 <= qscore <= 61
                    score = Math.Max(2, score);
                    score = Math.Min(61, score);
                    return (int) Math.Round(score, 0, MidpointRounding.AwayFromZero);
                default:
                    throw new Exception("Unhandled qscore method");
            }
        }

        /// <summary>
        /// Computes QScore predictor
        /// </summary>
        public enum QScorePredictor
        {
            BinCount,
            LogBinCount,
            BinCountAmpDistance,
            BinMean,
            BinCv,
            MafCount,
            MafMean,
            MafCv,
            LogMafCv,
            ModelDistance,
            RunnerUpModelDistance,
            DistanceRatio,
            CopyNumber,
            MajorChromosomeCount
        };

        public double GetQScorePredictor(QScorePredictor predictorId)
        {
            switch (predictorId)
            {
                case QScorePredictor.BinCount:
                    return (double) this.BinCount;

                case QScorePredictor.LogBinCount:
                    return Math.Log10(1 + this.BinCount);

                case QScorePredictor.BinCountAmpDistance:
                    return this.CopyNumber >= 15 ? Math.Log10(1 + this.BinCount) : 0.0;

                case QScorePredictor.BinMean:
                    if (this.Counts.Count == 0) return 0;
                    return this.Counts.Average();

                case QScorePredictor.BinCv:
                    if (this.Counts.Count == 0) return 0;
                    if (this.Counts.Average() == 0) return 0;
                    return Utilities.CoefficientOfVariation(this.Counts);

                case QScorePredictor.MafCount:
                    return Balleles.Frequencies.Count;

                case QScorePredictor.MafMean:
                    if (Balleles.Frequencies.Count == 0) return 0;
                    return Balleles.Frequencies.Average();

                case QScorePredictor.MafCv:
                    if (Balleles.Frequencies.Count == 0) return 0;
                    if (Balleles.Frequencies.Average() == 0) return 0;
                    return Utilities.CoefficientOfVariation(Balleles.Frequencies);

                case QScorePredictor.LogMafCv:
                    return Math.Log10(1 + GetQScorePredictor(QScorePredictor.MafCv));

                case QScorePredictor.ModelDistance:
                    // HACK: high copy number needs some help: higher variance and general lack of interest in exact copy number  
                    return this.ModelDistance / Math.Max(1.0,this.CopyNumber-4.0);

                case QScorePredictor.RunnerUpModelDistance:
                    return this.RunnerUpModelDistance;

                case QScorePredictor.DistanceRatio:
                    if (this.RunnerUpModelDistance == 0) return 0;
                    return this.ModelDistance / this.RunnerUpModelDistance;

                case QScorePredictor.CopyNumber:
                    return (double) this.CopyNumber;

                case QScorePredictor.MajorChromosomeCount:
                    // Force a double:
                    if (!this.MajorChromosomeCount.HasValue) return Math.Ceiling(this.CopyNumber / 2f);
                    return (double) this.MajorChromosomeCount;
            }
            return 0;
        }
    }
}